"""
 Defines an equation of state (eos) class which reads in a stellarcollapse.org
 EOS .h5 file.  Then any quantity in the table can be queried.
 Note: You may specify the independent variables rho and temp
  for physical state in terms of either logspace or non-logspace;
  Code which does table computation should always call validatePointDict
  to convert rho and temp to logspace if they are specified in non-logspace.

Jeff Kaplan  Feb, 2013  <jeffkaplan@caltech.edu>
"""
import ast
import h5py
import math
import numpy
import consts
from consts import CGS_C, CGS_H, CGS_EV, N_AVO
from utils import multidimInterp, linInterp, solveRootBisect, \
    relativeError, lookupIndexBisect
import scipy.optimize as scipyOptimize

#Plank's constant * speed of light in units of MeV * cm
hc_mevcm = CGS_H / (CGS_EV * 1.e6) * CGS_C

class eosDriver(object):

    valuesDict = None
    # Independent variables; assumes these are the independent variables in the
    # table.  Routines check input to see if variables are non-log10'd
    # independent variables are REORDERED in __init__ so that they
    # correspond to the correct ordering of the dependent variable table axes.
    indVars = ('ye', 'logtemp', 'logrho')
    # Physical state contains the  values of the independent variables in the
    # order determined in __init__
    physicalState = (None, None, None)

    #independent variables stored as logvars in table
    indLogVars =('rho', 'temp')
    logVars = ('energy', 'press', 'rho', 'temp')

    h5file = None

    #dependent variables table shape; set as shape of 'logpress'
    #refactor to remove tableShape and only have tableShapeDict!
    tableShape = None
    tableShapeDict = None

    #energy shift is constant value added to all values of 'energy' in the table
    #so that 'logenergy' is always positive
    energy_shift = None

    def __init__(self, tableFilename):
        """
        eosDriver class constructor takes a h5 file filename for a
        stellarcollapse.org EOS file
        """
        self.h5file = h5py.File(tableFilename, 'r')
        self.tableShape = numpy.shape(self.h5file['logpress'])
        self.energy_shift = self.h5file['energy_shift'][0]
        #print self.energy_shift
        #Determine the ordering of independent variable axes by identifying with
        # the number of points for that indVar axis
        newOrdering = [None for _ in self.indVars]
        for indVar in self.indVars:
            key = 'points' + indVar.split('log')[-1]
            points = self.h5file[key][0]
            for ithAxis, ithAxesPoints in enumerate(self.tableShape):
                if ithAxesPoints == points:
                    newOrdering[ithAxis] = indVar
                    break
        self.indVars = tuple(newOrdering)
        self.tableShapeDict = dict([(indVar, self.tableShape[i])
                                    for i, indVar in enumerate(self.indVars)])

    def __del__(self):
        """
        Without this, shit went very bad with multithreading in the thermalSupport module
        """
        self.h5file.close()

    def validatePointDict(self, pointDict):
        """
        Checks if points in pointDict are valid indVars.
        Converts non-log10'd vars to logvars
        !! MUTATES pointDict !!
        """
        assert isinstance(pointDict, dict)
        for key in pointDict.keys():
            assert key in self.indVars or key in self.indLogVars, \
                "'%s' of your pointDict is not a valid independent variable!" \
                "pointDict: %s \n" \
                "self.indVars: %s \n" \
                "self.indLogVars: %s" % (key, pointDict, self.indVars, self.indLogVars)

            assert isinstance(pointDict[key], float), \
                "Wow, it seems you've tried to specify %s as something " \
                "other than a float: %s" % (key, str(pointDict[key]))
            if key in self.indLogVars:
                pointDict.update({'log' + key: numpy.log10(pointDict[key])})
                del pointDict[key]

    def tempOfLog10RhobFuncFromPrescription(self, tempPrescription, ye):
        """
        Three temperature prescriptions available:
        1) Fix a quantity to determine T
         {'quantity':  a dependent variable in the EOS table,
          'target':    the target value you wish to fix 'quantity to'}
          E.g.: {'quantity': 'entropy', 'target': 1.0}

        2) Isothermal with a temperature roll-off in log10(rho_b CGS) space
         {'T': Isothermal temperature (MeV) at high densities (eosTmax),
          'rollMid': Midpoint of roll-off in log10( rho_b CGS) space,
          'rollScale': e-folding falloff of T to min value log10(rho_b CGS),
          'eosTmin': minimum temperature (MeV) in roll-off prescription }
          E.g: {'T': 30.0, 'rollMid': 14.0, 'rollScale': 0.5, 'eosTmin': 0.5}

        3) ManualTofLogRho
           Sets temperature using function specified.
           {'funcTofLogRho': string specifying python function to use
              must be valid python function that exists in eosDriver.py
              and is added to databaseOfManualFunctions! }

        Please only set one or the other of these; doing otherwise may result
        in UNDEFINED BEHAVIOR
        """
        databaseOfManualFunctions = dict([(f.__name__, f) for f in
                                          (kentaDataTofLogRhoFit1,
                                           kentaDataTofLogRhoFit2
                                          )
                                         ])
        isothermalKeys = ('T', 'rollMid', 'rollScale', 'eosTmin')
        fixedQuantityKeys = ('quantity', 'target')
        manualTofLogRhoKeys = ('funcTofLogRho',)
        assert isinstance(tempPrescription, dict)
        isothermalPrescription = all([key in tempPrescription.keys()
                                      for key in isothermalKeys])
        fixedQuantityPrescription = all([key in tempPrescription.keys()
                                         for key in fixedQuantityKeys])
        manualTofLogRhoPrescription = all([key in tempPrescription.keys()
                                           for key in manualTofLogRhoKeys])
        assert not(isothermalPrescription and fixedQuantityPrescription), "See docstring!"
        assert not(manualTofLogRhoPrescription and fixedQuantityPrescription), "See docstring!"
        assert not(manualTofLogRhoPrescription and isothermalPrescription), "See docstring!"

        #defines 1D root solver to use in routine
        solveRoot = scipyOptimize.brentq  # solveRootBisect
        tol = 1.0e-6

        tempOfLog10Rhob = lambda lr: None

        if isothermalPrescription:
            print "Identified isothermal prescription in tempOfLog10RhobFuncFromPrescription"
            Tmax = tempPrescription['T']
            Tmin = tempPrescription['eosTmin']
            mid = tempPrescription['rollMid']
            scale = tempPrescription['rollScale']
            tempOfLog10Rhob = getTRollFunc(Tmax, Tmin, mid, scale)
        elif manualTofLogRhoPrescription:
            funcName = tempPrescription['funcTofLogRho']
            assert funcName in databaseOfManualFunctions, \
                "funcName for manualTofLogRhoPrescription not recognized! If you have written it " \
                "be sure to add it to databaseOfManualFunctions at the start of writeRotNSeosfile"
            tempOfLog10Rhob = databaseOfManualFunctions[funcName]()
            print "Using manual prescription %s in tempOfLog10RhobFuncFromPrescription" % funcName
        elif fixedQuantityPrescription:
            print "Using fixed quantity prescription in tempOfLog10RhobFuncFromPrescription"
            quantity = tempPrescription['quantity']
            target = tempPrescription['target']
            if ye is not None:
                def getTemp(t, lr):
                    answer = multidimInterp((ye, t, lr),
                                            [self.h5file['ye'][:],
                                             self.h5file['logtemp'],
                                             self.h5file['logrho']],
                                            self.h5file[quantity][...],
                                            linInterp, 2) - target
                    #print t, lr, answer
                    return answer

                def solveTemp(lr):
                    try:
                        answer = solveRoot(lambda T: getTemp(T, lr),
                                           self.h5file['logtemp'][0],
                                           self.h5file['logtemp'][-1],
                                           (), tol)
                    except ValueError as err:
                        print "Root for log10(T) not bracketed on entire table! " + str(err)
                        answer = self.h5file['logtemp'][0]
                        print "Recovering with lowest table value, answer: %s" % answer
                    return numpy.power(10.0, answer)
                tempOfLog10Rhob = solveTemp
            else:
                # Otherwise tempOfLog10Rhob will return None,
                # this will trigger use of setConstQuantityAndBetaEq
                # in the output loop
                pass
        return tempOfLog10Rhob

    def writeRotNSeosfile(self, filename, tempPrescription, ye=None, addNeutrinoTerms=False):
        """
        Not setting ye will give results for neutrinoless beta equilibrium.
        Setting ye to NuFull will give results for nu-full beta equlibrium.
        Nu-full beta eq not supported with set const quantity!

        Add neutrino terms will add neutrino pressure according to eq (3) of
        NSNSThermal support paper via queryTrappedNuPress
        It will also add the contribution of the neutrinos to the energy density
        via eps_nu = 3 P_nu / rho
        """
        assert not(ye == 'NuFull' and 'quantity' not in tempPrescription), \
            "NuFull Beta-Eq not supported for fixedQuantityPrescription"

        tempOfLog10Rhob = self.tempOfLog10RhobFuncFromPrescription(tempPrescription, ye)

        log10numberdensityMin = 2.67801536139756E+01
        log10numberdensityMax = 3.97601536139756E+01  # = 1e16 g/cm^3

        npoints = 600

        #dlogn = (log10numberdensityMax - log10numberdensityMin) / (npoints - 1.0)

        logns = numpy.linspace(log10numberdensityMin, log10numberdensityMax, npoints)

        outfile = open(filename, 'w')

        for logn in logns:
            numberdensityCGS = numpy.power(10.0, logn)
            rho_b_CGS = numberdensityCGS * consts.CGS_AMU
            logrho_b_CGS = numpy.log10(rho_b_CGS)

            temp = tempOfLog10Rhob(logrho_b_CGS)

            if ye is None and temp is not None:
                self.setBetaEqState({'rho':rho_b_CGS, 'temp': temp})
            elif ye is None and temp is None:
                quantity = tempPrescription['quantity']
                target = tempPrescription['target']
                self.setConstQuantityAndBetaEqState({'rho': rho_b_CGS},
                                                    quantity,
                                                    target)
            elif ye == 'NuFull':
                self.setNuFullBetaEqState({'rho': rho_b_CGS, 'temp': temp})
            else:
                self.setState({'rho': rho_b_CGS, 'temp': temp, 'ye': ye})

            logpress, logeps, munu = self.query(['logpress', 'logenergy', 'munu'])
            #print logpress, rho_b_CGS
            P_nu = 0.0
            if addNeutrinoTerms:
                P_nu = P_nu_of(rho_b_CGS, temp, munu)
            logpress = numpy.log10(numpy.power(10.0, logpress) + P_nu)
            eps_nu = 3.0 * P_nu / rho_b_CGS

            #Total energy density (1.0 + eps) is in AGEO units
            eps = (numpy.power(10.0, logeps) - self.energy_shift + eps_nu)\
                  * consts.AGEO_ERG / consts.AGEO_GRAM
            #print logpress
            #print
            totalEnergyDensity = rho_b_CGS * (1.0 + eps)
            logTotalEnergyDensity = numpy.log10(totalEnergyDensity)
            #print logrho_b_CGS ,logeps,  numpy.power(10.0, logeps), eps, self.energy_shift
            #print logrho_b_CGS, eps, logpress
            outfile.write("{:24.14e}{:24.14e}{:24.14e}\n".format(logn,
                                                                 logTotalEnergyDensity,
                                                                 logpress))

    #todo: add option for picking different recovery methods
    def solveForQuantity(self, pointDict, quantity, target, bounds=None,
                         function=(lambda x, q: q),
                         pointAsFunctionOfSolveVar=lambda x: None,
                         tol=1.e-6):
        """
        Solve for independent variable left out of pointDict so that
        quantity=target.
        If bounds for root solve not supplied, will try the table max and min.
        Function lets the user specify an arbitrary function of the independent
        variable, x, and the quantity q as what we are solving for.
        Defaults to simply the quantity.
        pointAsFunctionOfSolveVar lets you specify logtemp as a function of
        solveVar.  Right now hardcoded so only uses logtemp :(
        """
        assert isinstance(pointDict, dict)
        self.validatePointDict(pointDict)
        assert len(pointDict) < 3, "Can't solve anything if you've specified more than 2 indVars!"
        assert len(pointDict) > 1, "Solve is under-determined with less than 2 indVars!"

        solveRoot = scipyOptimize.brentq
        #solveRoot = solveRootBisect
        solveVar = [indVar for indVar in self.indVars if indVar not in pointDict][0]

        #NOTE POINTASFUNCTIONOFSOLVERDICT MUST BE IN SAME FORMAT AS POINTDICT
        #TODO FIX THIS HARD CODING FUCK FUKC FUCK
        if pointAsFunctionOfSolveVar(14.0) is None:
            val = pointDict['logtemp']
            pointAsFunctionOfSolveVar = lambda x: val
        #todo: add some good asserts for bounds
        #NOTE BOUNDS MUST BE IN LOGVAR!!!
        if bounds is not None:
            boundMin = bounds[0]
            boundMax = bounds[1]
        else:
            boundMin = self.h5file[solveVar][0]
            boundMax = self.h5file[solveVar][-1]

        indVarsTable = self.getIndVarsTable()

        def quantityOfSolveVar(x):
            #Here we construct the point to interpolate at, but we
            # must do it carefully since we don't know apriori what
            # solveVar is
            point = []
            #todo factor this for out of quantityOfSolveVar
            for indVar in self.indVars:
                if indVar not in pointDict:
                    #print "NOT", indVar
                    value = x
                else:
                    value = pointDict[indVar]
                    if indVar == 'logtemp':
                        value = pointAsFunctionOfSolveVar(x)
                #print indVar, value
                point.append(value)
            point = tuple(point)
            #print point
            #print indVarsTable
            answer = function(x, multidimInterp(point, indVarsTable,
                                                self.h5file[quantity][...],
                                                linInterp, 2)
                              ) - target
            return answer

        try:
            answer = solveRoot(quantityOfSolveVar, boundMin, boundMax, (), tol)
        except ValueError as err:
            #todo: note this is slightly incorrect if pointAsFunctionOfSolveVar is specified
            print "Error in root solver solving for %s: " % solveVar, str(err)
            answer = self.findIndVarOfMinAbsQuantity(solveVar,
                                                     self.pointFromDict(pointDict),
                                                     quantity,
                                                     function,
                                                     target)
            print "Recovering with findIndVarOfMinAbsQuantity, answer: %s" % answer
        return answer

    def tableIndexer(self, indVar, i):
        """
        Given an indVar and an index for that indVar returns
        a numpy array slice tuple which indexes the slice for
        the dependent variable table with indVar fixed to index i
        """
        assert indVar in self.indVars, "Input indVar %s is not a valid indVar!" % indVar
        assert i >= 0, "Index in tableIndexer, %s, is less than zero. Dummy." % i
        assert i < self.tableShapeDict[indVar], \
            "Index in tableIndexer, %s, is greater than table extends in indVar direction!" % i

        result = []
        for var in self.indVars:
            if indVar == var:
                result.append(i)
            else:
                result.append(slice(None))
        return tuple(result)

    def getIndVarsTable(self, omitTheseIndVars=()):
        """
        Returns a table where the each entry is the list of table point
        values for each indVar.  Ordered list.
        """
        assert all([indVar in self.indVars for indVar in omitTheseIndVars]), \
            "Can't omit an indVar that is not a proper indVar!"

        result = []
        for indVar in self.indVars:
            if indVar not in omitTheseIndVars:
                result.append(self.h5file[indVar])
        #returning a tuple prevents inadvertent mutating of the result
        return tuple(result)

    def pointFromDict(self, pointDict):
        """
        Returns tuple point of indVars in their CORRECT ORDERING
        given a pointDict.  Also DELOGS logvars
        """
        assert isinstance(pointDict, dict)
        assert all([key in self.indVars for key in pointDict.keys()])

        result = []
        for indVar in self.indVars:
            if indVar in pointDict:
                if indVar in self.indLogVars:
                    result.append(numpy.log10(pointDict[indVar]))
                else:
                    result.append(pointDict[indVar])

        return tuple(result)

    def findIndVarOfMinAbsQuantity(self, indVar, point, quantity,
                                   function=(lambda x, q: q), target=0.0):
        """
        Given an independent variable indVar,
         The values of the other independent variables as point,
         and a dependent variable quantity,
        Find for what value of indVar's table grid-points is
        quantity closest to zero.  Uses simple sequential search.
        Designed for use in conjunction with solveForQuantity.
        If function is specified, it searches for
        function(indVar, quantity closest to zero)
        If target specified finds closest value to target
        """
        assert indVar in self.indVars, "Input indVar %s is not a valid indVar!" % indVar

        indVarsTable = self.getIndVarsTable(omitTheseIndVars=(indVar,))
        closestIndVar = None
        closestQuantity = 1.0e300

        for i, var in enumerate(self.h5file[indVar][:]):
            index = self.tableIndexer(indVar, i)
            thisQuantity = function(var,
                                    multidimInterp(point, indVarsTable,
                                                   self.h5file[quantity][index],
                                                   linInterp, 2)
                                    ) - target
            if abs(thisQuantity) < closestQuantity:
                closestIndVar = var
                closestQuantity = abs(thisQuantity)
        return closestIndVar

    def setNuFullBetaEqState(self, pointDict, Ylep=None):
        """
        Set's 'neutrino-full' beta equilibrium for a rho/temp point dict
        and a lepton fraction.  If no lepton fraction defined, it will
        be calculated as the value at cold neutrino-less beta eq.
        """
        assert isinstance(pointDict, dict)
        assert 'ye' not in pointDict, "Can't set ye since we're solving for it!"
        self.validatePointDict(pointDict)

        logtempMin = self.h5file['logtemp'][0]
        if Ylep is None:
            coldBetaEqDict = pointDict.copy()
            coldBetaEqDict['logtemp'] = logtempMin
            Ylep = self.setBetaEqState(coldBetaEqDict)
            #We don't actually want to use this state, so clear it to prevent
            # unintended use
            self.clearState()

        temp = numpy.power(10.0, pointDict['logtemp'])
        rho = numpy.power(10.0, pointDict['logrho'])

        def YePlusYnu(ye, munu):

            eta = munu / temp 
            n_nu = 4 * numpy.pi * (temp / hc_mevcm) ** 3 \
                * 1.0 / 3.0 * eta * (eta ** 2 + numpy.pi ** 2)
            Ynu = n_nu / (rho * N_AVO) 

            return ye + Ynu

        currentYe = self.solveForQuantity(pointDict, 'munu', Ylep, function=YePlusYnu)
        newDict = pointDict.copy()
        newDict['ye'] = currentYe
        self.setState(newDict)
        return currentYe

    def rhobFromEnergyDensity(self, ed, pointDict):
        assert isinstance(pointDict, dict)
        assert 'rho' not in pointDict and 'logrho' not in pointDict, \
            "Can't set rho/logrho since we're solving for it!"
        self.validatePointDict(pointDict)

        edFunc = lambda x, q: numpy.power(10.0,x) \
                              * (1.0 + (numpy.power(10.0, q) - self.energy_shift)/ CGS_C**2)
        result = self.solveForQuantity(pointDict, 'logenergy', ed,
                                       bounds=(4., 16.),
                                       function=edFunc)
        return result


    def rhobFromEnergyDensityWithTofRho(self, ed, ye, TofRhoFunc):

        edFunc = lambda x, q: numpy.power(10.0,x) \
                              * (1.0 + (numpy.power(10.0, q) - self.energy_shift)/ CGS_C**2)
        tempFunc = lambda x: numpy.log10(TofRhoFunc(x))
        result = self.solveForQuantity({'logtemp': 0.0, 'ye': ye}, 'logenergy', ed,
                                       bounds=(4., 16.),
                                       function=edFunc,
                                       pointAsFunctionOfSolveVar=tempFunc)
        return result

    def setState(self, pointDict):
        """
        Takes dictionary defining a physical state, aka the EOS table's
        independent variables, and sets the physical state.
        Modifies self.physicalState
        """
        self.validatePointDict(pointDict)
        state = []
        for indVar in self.indVars:
            assert indVar in pointDict,\
                "You have not specified an required independent variable in pointDict!"
            state.append( pointDict[indVar] )
        self.physicalState = tuple(state)

    def clearState(self):
        """
        Resets the physicalState member variable to Nones.
        Should prevent the query member from being called.
        Modifies self.physicalState
        """
        self.physicalState = (None for unused in self.indVars)


    def setConstQuantityState(self, pointDict, quantity, target):
        assert len(pointDict) < 3, "State overdetermined for more than 2 indVars!"
        self.validatePointDict(pointDict)
        assert False, "setConstQuantityState not implemented yet!"

    def getTemperatureFromQuantityTYe(self, pointDict, quantity, target):
        # assign vars
        self.validatePointDict(pointDict)
        ye = pointDict['ye']
        xt = pointDict['logtemp']
        lr = pointDict['logrho']

        #defines 1D root solver to use in routine
        solveRoot = scipyOptimize.brentq  # solveRootBisect
        tol = 1.0e-12

        def getTemp(t, lr):
            answer = multidimInterp((ye, t, lr),
                                    [self.h5file['ye'][:],
                                     self.h5file['logtemp'],
                                     self.h5file['logrho']],
                                    self.h5file[quantity][...],
                                    linInterp, 2) - target
            return answer

        def solveTemp(lr):
            try:
                answer = solveRoot(lambda T: getTemp(T, lr),
                                   self.h5file['logtemp'][0],
                                   self.h5file['logtemp'][-1],
                                   (), tol)
            except ValueError as err:
                print "Root for log10(T) not bracketed on entire table! " + str(err)
                # see if lower or upper temperature bound best
                logtemp = self.h5file['logtemp']
                answer1 = multidimInterp((ye, logtemp[0], lr),
                                         [self.h5file['ye'][:],
                                          self.h5file['logtemp'],
                                          self.h5file['logrho']],
                                         self.h5file[quantity][...],
                                         linInterp, 2) - target
                answer2 = multidimInterp((ye, logtemp[-1], lr),
                                         [self.h5file['ye'][:],
                                          self.h5file['logtemp'],
                                          self.h5file['logrho']],
                                         self.h5file[quantity][...],
                                         linInterp, 2) - target

                if (abs(answer1) < abs(answer2)):
                    answer = self.h5file['logtemp'][0]
                    print "Recovering with lowest table value, answer: %s" % answer
                else:
                    answer = self.h5file['logtemp'][-1]
                    print "Recovering with highest value, answer: %s" % answer

            return numpy.power(10.0, answer)

        return solveTemp(lr)
        
    #todo: RIGHT NOW HARD CODED TO BE GIVEN RHO AND FIND T!! FIX
    #     solveVar is T and otherVar is rho!
    def setConstQuantityAndBetaEqState(self, pointDict, quantity, target):
        """
        Does inefficient 2D root solve to set state at neutrino-less
        beta equilibrium and quantity = target. Sets the physicalState
         and then returns ye, temp
        Modifies self.physicalState
        """
        print "setConstQuantityAndBetaEqState: ", pointDict
        assert 'ye' not in pointDict, "You can't SPECIFY a Ye if you're " \
                                      "setting neutrinoless beta equlibrium!"
        self.validatePointDict(pointDict)
        assert len(pointDict) < 2, "State overdetermined for more than 1 indVars!"
        #todo: check quantity is valid 3D table

        #defines 1D root solver to use in routine
        solveRoot = scipyOptimize.brentq  # solveRootBisect

        solveVarName = 'logtemp'
        currentSolveVar =  0.0
        currentYe = 0.25
        #previous variables used to measure convergence of solve
        # so set them to something significantly different than starting values
        previousSolveVar =  100.0
        previousYe = 100.0
        yeError = relativeError(currentYe, previousYe)
        solveVarError = relativeError(currentSolveVar, previousSolveVar)
        otherVarName = pointDict.keys()[0]
        otherVar = pointDict.values()[0]

        maxIters = 5
        tol = 1e-3

        iteration = 0
        while iteration < maxIters and yeError + solveVarError > tol/2.0:
            previousSolveVar = currentSolveVar
            previousYe = currentYe
            getSolveVar = lambda x: multidimInterp((currentYe, x, otherVar),
                                                   [self.h5file['ye'][:],
                                                    self.h5file[solveVarName],
                                                    self.h5file[otherVarName]],
                                                   self.h5file[quantity][...],
                                                   linInterp, 2) - target
            try:
                currentSolveVar = solveRoot(getSolveVar,
                                            self.h5file[solveVarName][0],
                                            self.h5file[solveVarName][-1],
                                            (),tol)
            except ValueError as err:
                print "Root for log10(T) not bracketed on entire table: " \
                      + str(err)
                # see if lower or upper temperature bound best
                logtemp = self.h5file['logtemp']
                answer1 = multidimInterp((currentYe, logtemp[0], otherVar),
                                         [self.h5file['ye'][:],
                                          self.h5file['logtemp'],
                                          self.h5file['logrho']],
                                         self.h5file[quantity][...],
                                         linInterp, 2) - target
                answer2 = multidimInterp((currentYe, logtemp[-1], otherVar),
                                         [self.h5file['ye'][:],
                                          self.h5file['logtemp'],
                                          self.h5file['logrho']],
                                         self.h5file[quantity][...],
                                         linInterp, 2) - target

                if (abs(answer1) < abs(answer2)):
                    currentSolveVar = self.h5file['logtemp'][0]
                    print "Recovering with lowest table value, answer: %s" % currentSolveVar
                else:
                    currentSolveVar = self.h5file['logtemp'][-1]
                    print "Recovering with highest value, answer: %s" % currentSolveVar

            getYe = lambda x : multidimInterp((x, currentSolveVar, otherVar),
                                              [self.h5file['ye'][:],
                                               self.h5file[solveVarName],
                                               self.h5file[otherVarName]],
                                              self.h5file['munu'][...],
                                              linInterp, 2)
            #check for bracketing error in root solve for ye
            try:
                currentYe = solveRoot(getYe,
                                      self.h5file['ye'][0],
                                      self.h5file['ye'][-1], (), tol)
            except ValueError as err:
                print "Error in scipy root solver solving for ye: ", str(err)
                currentYe = self.findYeOfMinAbsMunu((currentSolveVar, otherVar))
                print "Recovering with findYeOfMinAbsMunu, answer: %s" % currentYe
            #print "currentYe: ", currentYe, "\tcurrentT: ", currentSolveVar

            yeError = relativeError(currentYe, previousYe)
            solveVarError = relativeError(currentSolveVar, previousSolveVar)
            iteration += 1
            #print "errs: ", yeError, solveVarError

        newDict = pointDict.copy()
        newDict['ye'] = currentYe
        temp = numpy.power(10.0,currentSolveVar)  # TODO TEMP HARD CODE
        newDict['temp'] = temp
        self.setState(newDict)
        return currentYe, temp # TODO TEMP HARD CODE

    def findYeOfMinAbsMunu(self, point):
        """
        Given a point in T, rho, calculate for what value of the Ye
        table grid-points is munu the closest to zero.
        """
        closestYeToMunusZero = None
        closestMunuToZero = 1.0e300

        for i, ye in enumerate(self.h5file['ye'][:]):
            munu = multidimInterp(point, [self.h5file['logtemp'],
                                          self.h5file['logrho']],
                                  self.h5file['munu'][i, ...],
                                  linInterp, 2)
            if abs(munu) < closestMunuToZero:
                closestYeToMunusZero = ye
                closestMunuToZero = abs(munu)
        return closestYeToMunusZero

    # neutrino-less beta equilibrium occurs when mu_n = mu_e + mu_p
    # munu = mu_p - mu_n + mu_e, so when munu = 0, we have beta-eq!
    def setBetaEqState(self, pointDict, useThisYeIfSolveFails=None):
        """
        Takes dictionary for physical state values EXCEPT Ye, then sets Ye
        via solving for neutrino-less beta equilibrium.
        Modifies self.physicalState
        """
        assert isinstance(pointDict, dict)
        assert 'ye' not in pointDict, "You can't SPECIFY a Ye if you're " \
                                      "setting neutrinoless beta equlibrium!"
        self.validatePointDict(pointDict)
        assert len(pointDict) < 3, "State overdetermined for more than 2 indVars!"

        #defines 1D root solver to use in routine
        solveRoot = scipyOptimize.brentq  # solveRootBisect


        #ASSUME 2 INDEPENENT VARIABLES ARE rho & temp
        logtemp = pointDict['logtemp']
        logrho = pointDict['logrho']

        tol = 1.e-6
        getYe = lambda x : multidimInterp((x, logtemp, logrho),
                                          [self.h5file['ye'][:],
                                           self.h5file['logtemp'],
                                           self.h5file['logrho']],
                                          self.h5file['munu'][...],
                                          linInterp, 2)
        #check for bracketing error in root solve for ye
        try:
            currentYe = solveRoot(getYe,
                                  self.h5file['ye'][0],
                                  self.h5file['ye'][-1], (), tol)
        except ValueError as err:
            print "Error in scipy root solver solving for ye: ", str(err)
            if useThisYeIfSolveFails is None:
                currentYe = self.findYeOfMinAbsMunu((logtemp, logrho))
                print "Recovering with findYeOfMinAbsMunu, answer: %s" % currentYe
            else:
                currentYe = useThisYeIfSolveFails
                print "Setting Ye to useThisYeIfSolveFails, answer: %s" % currentYe

        newDict = pointDict.copy()
        newDict['ye'] = currentYe
        self.setState(newDict)
        return currentYe

    def queryTrappedNuPress(self, rho_trap=10. ** 12.5):
        """
        Queries the EOS table for the pressure due to trapped neutrinos.
        See eq (2) & (3) of NSNSthermal paper. rho_trap should be in g/cm^3
        Note: calls self.query!
        """
        # This is a way to recover the independent variables from the class
        # given that they are ordered properly in setState
        pointDict = dict()
        for i, var in enumerate(self.indVars):
            pointDict[var] = self.physicalState[i]
        temp_mev = numpy.power(10.0, pointDict['logtemp'])
        rho_cgs = numpy.power(10.0, pointDict['logrho'])

        munu_mev = self.query('munu')

        return P_nu_of(rho_cgs, temp_mev, munu_mev, rho_trap)

    #TODO: query should check to make sure quantity is a valid quantity in the h5file
    def query(self, quantities, deLog10Result=False):
        """
        Queries the EOS table looking for 'quantity' at set physical state
        Note: query clears physical state after quantity is determined!
        """
        assert all([var is not None for var in self.physicalState]), \
            "One or more independent variables " \
            "not set for this EOS's physical state!"

        tableIndexes = []
        for i, indVar in enumerate(self.indVars):
            value = self.physicalState[i]
            tableIndexes.append(lookupIndexBisect(value, self.h5file[indVar][:]))

        answers = self.interpolateTable(tableIndexes, quantities)
        self.clearState()
        if deLog10Result:
            answers = numpy.power(10.0, answers)
        return answers

    #just does trilinear interpolation; does NOT try and account/correct
    # for log spacing in independent or dependent variable
    # got lazy; am not generalizing, copying from
    # http://en.wikipedia.org/wiki/Trilinear_interpolation
    # Note: index lookup is separated from interpolation routine
    #       so that different interpolators may be written and plugged in
    def interpolateTable(self, tableIndex, quantities):
        """
        Given table indexes preceding location of independent variables
        (physical state), does trilinear interpolation of 'quantity' to the
         current physical state.
        """
        tableIndex = tuple(tableIndex)
        assert all([tableIndex[i] +1 < self.tableShape[i] for i in range(len(tableIndex))]), \
            "Table index + 1  specified to interpolateTable is out of tableShape range!"
        if isinstance(quantities,str):
            quantities = [quantities]
        ys = []
        #todo: maybe implement checking if quantity is logvar
        for quantity in quantities:
            ys.append(self.h5file[quantity])

        xs = []  # xs is vector x, y, z from wikipedia
        x0 = []  # x0 is vector x0, y0, z0 from wikipedia
        x1 = []  # x0 is vector x1, y1, z1 from wikipedia
        for i, indVar in enumerate(self.indVars):
            value = self.physicalState[i]
            xs.append( value  )
            x0.append( self.h5file[indVar][tableIndex[i]] )
            x1.append( self.h5file[indVar][tableIndex[i]+1] )
        xs = numpy.array(xs)
        x0 = numpy.array(x0)
        x1 = numpy.array(x1)

        dxs = (xs - x0) / (x1 - x0)  # dxs is xd, yd, zd from wikipedia

        answers=[]
        for y in ys:
            #matrix is c[i,j] in wikipedia
            matrix = numpy.zeros([2,2])
            for i in (0,1):
                for j in (0,1):
                    lowerInX = (tableIndex[0],tableIndex[1]+i,tableIndex[2]+j)
                    upperInX = (tableIndex[0] +1,tableIndex[1]+i,tableIndex[2]+j)
                    #print lowerInX, upperInX, y[lowerInX], y[upperInX]
                    matrix[i,j] = y[lowerInX] * (1.-dxs[0]) + y[upperInX]* dxs[0]
            #print matrix

            vector = numpy.zeros([2])
            for i in (0,1):
                vector[i] = matrix[0,i] * (1.-dxs[1]) + matrix[1,i] * dxs[1]
            #print vector
            answers.append(vector[0] * (1. - dxs[2]) + vector[1] * dxs[2])

        if len(answers) == 1:
            return answers[0]
        else:
            return answers


def P_nu_of(rho_cgs, temp_mev, munu_mev, rho_trap=10 ** 12.5):
    """
    Returns trapped neutrino pressure as function of variables it
    depends on.  Default trapping density is fiducial value of
    10^12.5 g/cm^3.
    """
    eta_nu = munu_mev / temp_mev

    # First term
    P_nu = 4.0 * numpy.pi / 3.0 * temp_mev ** 4.0 / hc_mevcm ** 3.0
    # Fermi integral term
    P_nu *= 21.0 / 60.0 * numpy.pi ** 4.0 \
            + 0.5 * eta_nu ** 2 * (numpy.pi ** 2 + 0.5 * eta_nu ** 2)
    # Decoupling at lower densities
    P_nu *= numpy.exp(-rho_trap / rho_cgs)

    # convert from MeV to erg and return
    return P_nu * (CGS_EV * 1.0e6)

def getTRollFunc(Tmax, Tmin, mid, scale):
    """
    Returns a tanh step function with a minimum of Tmin,
     a maximum of Tmax, a halfway point at mid and
     the transition e-folding length of scale
    """
    tempOfLog10Rhob = lambda lr: Tmin + (Tmax - Tmin) / 2.0 \
                                         * (numpy.tanh((lr - mid)/scale) + 1.0)
    return tempOfLog10Rhob

def kentaDataTofLogRhoFit1():
    """
    Returns a function that is parametrized fit for T(logrho)
     of Kenta's Shen135135 simulation data.
    """
    func = lambda lr: getTRollFunc(20.0, 0.0, 14.25 - 0.07, .25)(lr) \
                    + getTRollFunc(10.0, 0.01, 11.5, .25)(lr)
    return func

def kentaDataTofLogRhoFit2():
    """
    Returns a function that is parametrized fit for T(logrho)
     of Kenta's Shen135135 simulation data.
     EXCEPT that the temperature plateau has been cooled to 5 MeV
    """
    func = lambda lr: getTRollFunc(25.0, 0.0, 14.1875 - 0.07, .3125)(lr) \
                    + getTRollFunc(5.0, 0.01, 11.5, .25)(lr)
    return func
