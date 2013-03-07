"""
 Defines an equation of state (eos) class which reads in a stellarcollapse.org
 EOS .h5 file.  Then any quantity in the table can be queried.
 Note: specify the independent variables for physical state in non-logspace;
       the code checks to see if variable is in logspace and adjusts

Jeff Kaplan  Feb, 2013  <jeffkaplan@caltech.edu>
"""
import h5py
import math
import numpy
import consts
from utils import multidimInterp, linInterp, solveRootBisect,\
    BracketingError, relativeError, lookupIndexBisect
import scipy.optimize as scipyOptimize


class eosDriver(object):

    valuesDict = None
    # Independent variables; assumes these are the independent variables in the
    # table.  Routines check to see if they are stored as 'log' + var
    # independent variables are REORDERED in __init__ so that they
    # correspond to the correct ordering of the dependent variable table axes.
    indVars = ('rho', 'ye', 'temp')
    # Physical state contains the  values of the independent variables in the
    # order determined in __init__
    physicalState = (None, None, None)

    logVars = ('energy', 'press', 'rho', 'temp')

    h5file = None

    #dependent variables table shape; set as shape of 'logpress'
    tableShape = None

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
        newOrdering = [None for unused in self.indVars]
        for indVar in self.indVars:
            key = 'points' + indVar
            points = self.h5file[key][0]
            for ithAxis, ithAxesPoints in enumerate(self.tableShape):
                if ithAxesPoints == points:
                    newOrdering[ithAxis] = indVar
                    break
        self.indVars = tuple(newOrdering)

    def writeRotNSeosfile(self, filename, tempPrescription, ye=None):
        """
        Two temperature prescriptions available:
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

        Please only set one or the other of these; doing otherwise may result
        in UNDEFINED BEHAVIOR

        Not setting ye will give results for neutrinoless beta equilibrium.
        NOT IMPLEMENTED YET
        """

        isothermalKeys = ('T', 'rollMid', 'rollScale', 'eosTmin')
        fixedQuantityKeys = ('quantity', 'target')
        assert isinstance(tempPrescription, dict)
        isothermalPrescription = all([key in tempPrescription.keys()
                                      for key in isothermalKeys])
        fixedQuantityPrescription = all([key in tempPrescription.keys()
                                         for key in fixedQuantityKeys])
        assert not(isothermalPrescription and fixedQuantityPrescription), "See docstring!"

        #defines 1D root solver to use in routine
        solveRoot = scipyOptimize.brentq  # solveRootBisect
        tol = 1.0e-6

        tempOfLog10Rhob = lambda lr: None
        yeOfLog10Rhob = lambda ye: ye

        if isothermalPrescription:
            print "Using isothermal prescription in writeRotNSeosfile"
            Tmax = tempPrescription['T']
            Tmin = tempPrescription['eosTmin']
            mid = tempPrescription['rollMid']
            scale = tempPrescription['rollScale']
            tempOfLog10Rhob = lambda lr: Tmin + (Tmax - Tmin) / 2.0 \
                                         * (numpy.tanh((lr - mid)/scale) + 1.0)
        elif fixedQuantityPrescription:
            print "Using fixed quantity prescription in writeRotNSeosfile"
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
                self.setConstQuantityAndBetaEqState({'rho': rho_b_CGS},
                                                    quantity,
                                                    target)
            else:
                self.setState({'rho': rho_b_CGS, 'temp': temp, 'ye': ye})

            logpress, logeps = self.query(['logpress', 'logenergy'])

            #Total energy density (1.0 + eps) is in AGEO units

            eps = (numpy.power(10.0, logeps) - self.energy_shift)\
                  * consts.AGEO_ERG / consts.AGEO_GRAM

            totalEnergyDensity = rho_b_CGS * (1.0 + eps)
            logTotalEnergyDensity = numpy.log10(totalEnergyDensity)
            #print logrho_b_CGS ,logeps,  numpy.power(10.0, logeps), eps, self.energy_shift
            #print logrho_b_CGS, eps, logpress
            outfile.write("{:24.14e}{:24.14e}{:24.14e}\n".format(logn,
                                                                 logTotalEnergyDensity,
                                                                 logpress))


    def setState(self, pointDict):
        """
        Takes dictionary defining a physical state, aka the EOS table's
        independent variables, and sets the physical state.
        Modifies self.physicalState
        """
        state = []
        for indVar in self.indVars:
            assert indVar in pointDict,\
                "You have not specified an required independent variable in pointDict!"
            state.append( pointDict[indVar] )
        self.physicalState = list(state)

    def clearState(self):
        """
        Resets the physicalState member variable to Nones.
        Should prevent the query member from being called.
        Modifies self.physicalState
        """
        self.physicalState = (None for unused in self.indVars)


    def setConstQuantityState(self, pointDict, quantity, target):
        assert len(pointDict) < 3, "State overdetermined for more than 2 indVars!"
        assert False, "setConstQuantityState not implemented yet!"

    def getTemperatureFromQuantityTYe(self, pointDict, quantity, target):
        # assign vars
        ye = pointDict['ye']
        xt = numpy.log10(pointDict['temp'])
        lr = numpy.log10(pointDict['rho'])

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
        assert all([key in self.indVars for key in pointDict.keys()])
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
        if otherVarName in self.logVars:
                otherVar = math.log10(otherVar)
                otherVarName = 'log' + otherVarName

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
            except BracketingError as err:
                print "Root for ye not bracketed on entire table!" + str(err)
                currentYe = self.findYeOfMinAbsMunu((currentSolveVar, otherVar))
                print "Recovering with findYeOfMinAbsMunu, answer: %s" % currentYe
            #ValueError is thrown by scipy's brentq
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
        newDict['temp'] = numpy.power(10.0,currentSolveVar)  # TODO TEMP HARD CODE
        self.setState(newDict)
        return currentYe, newDict['temp'] # TODO TEMP HARD CODE

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
    def setBetaEqState(self, pointDict):
        """
        Takes dictionary for physical state values EXCEPT Ye, then sets Ye
        via solving for neutrino-less beta equilibrium.
        Modifies self.physicalState
        """
        assert isinstance(pointDict, dict)
        assert 'ye' not in pointDict, "You can't SPECIFY a Ye if you're " \
                                      "setting neutrinoless beta equlibrium!"
        assert all([key in self.indVars for key in pointDict.keys()])
        assert len(pointDict) < 3, "State overdetermined for more than 2 indVars!"

        #defines 1D root solver to use in routine
        solveRoot = scipyOptimize.brentq  # solveRootBisect


        for key, value in pointDict.items():
            if key in self.logVars:
                pointDict['log' + key] = numpy.log10(value)

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
        except BracketingError as err:
            print "Root for ye not bracketed on entire table!" + str(err)
            currentYe = self.findYeOfMinAbsMunu((logtemp, logrho))
            print "Recovering with findYeOfMinAbsMunu, answer: %s" % currentYe
        #ValueError is thrown by scipy's brentq
        except ValueError as err:
            print "Error in scipy root solver solving for ye: ", str(err)
            currentYe = self.findYeOfMinAbsMunu((logtemp, logrho))
            print "Recovering with findYeOfMinAbsMunu, answer: %s" % currentYe

        newDict = pointDict.copy()
        newDict['ye'] = currentYe
        self.setState(newDict)
        return currentYe

    #TODO: query should check to make sure quantity is a valid quantity in the h5file
    def query(self, quantities, deLog10Result=False):
        """
        Query's the EOS table looking for 'quantity' at set physical state
        Note: query clears physical state after quantity is determined!
        """
        assert all(self.physicalState), "One or more independent variables " \
                                        "not set for this EOS's physical state!"

        tableIndexes = []
        for i, indVar in enumerate(self.indVars):
            value = self.physicalState[i]
            if indVar in self.logVars:
                value = math.log10(value)
                indVar = 'log' + indVar
            tableIndexes.append(lookupIndexBisect(value, self.h5file[indVar][:]))

        answers = self.interpolateTable(tableIndexes, quantities)
        self.clearState()
        if deLog10Result:
            answers = numpy.power(10.0,answers)
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
            if indVar in self.logVars:
                value = math.log10(value)
                indVar = 'log' + indVar
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

        # print xs
        # print x0
        # print x1
        # print numpy.shape(self.h5file[quantity])
        # print self.h5file[quantity][tableIndex]

