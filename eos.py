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


class eos(object):

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

    def __init__(self, tableFilename):
        """
        eos class constructor takes a h5 file filename for a
        stellarcollapse.org EOS file
        """
        self.h5file = h5py.File(tableFilename, 'r')
        #Determine the ordering of independent variable axes by identifying with
        # the number of points for that indVar axis
        newOrdering = [None for unused in self.indVars]
        for indVar in self.indVars:
            key = 'points' + indVar
            points = self.h5file[key][0]
            for ithAxis, ithAxesPoints in enumerate(numpy.shape(self.h5file['logpress'])):
                if ithAxesPoints == points:
                    newOrdering[ithAxis] = indVar
                    break
        self.indVars = tuple(newOrdering)

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

    def setBetaEqState(self, pointDict):
        """
        Takes dictionary for physical state values EXCEPT Ye, then sets Ye
        via solving for neutrino-less beta equilibrium.
        Modifies self.physicalState
        """
        assert 'ye' not in pointDict, "You can't SPECIFY a Ye if you're " \
                                      "setting neutrinoless beta equlibrium!"
        assert all([key in self.indVars for key in pointDict.keys()])

        tableIndexes=[]
        partialNewState=[]
        for indVar in self.indVars:
            if indVar in pointDict:
                tableIndexes.append(self.lookupIndex(indVar, pointDict[indVar]))
                partialNewState.append(pointDict[indVar])
            elif indVar == 'ye':
                tableIndexes.append(None)
                partialNewState.append(None)
            else:
                assert False, "indVar %s is not ye or in pointDict!" % indVar
        # physical state for non-ye independent variables is required in
        # getYeBetaEqFromTable for interpolation
        newDict = pointDict.copy()
        newDict['ye'] = self.getYeBetaEqFromTable(tableIndexes, partialNewState)
        self.setState(newDict)

    #TODO: query should check to make sure quantity is a valid quantity in the h5file
    def query(self, quantity):
        """
        Query's the EOS table looking for 'quantity' at set physical state
        Note: query clears physical state after quantity is determined!
        """
        assert all(self.physicalState), "One or more independent variables " \
                                        "not set for this EOS's physical state!"

        tableIndexes = []
        for i, indVar in enumerate(self.indVars):
            value = self.physicalState[i]
            tableIndexes.append(self.lookupIndex(indVar, value))

        answer = self.interpolateTable(tableIndexes, quantity)
        self.clearState()
        return answer

    # neutrino-less beta equilibrium occurs when mu_n = mu_e + mu_p
    # munu = mu_p - mu_n + mu_e, so when munu = 0, we have beta-eq!
    # TODO: The logic of this routine combined with setBetaYe is not great; refactor it
    def getYeBetaEqFromTable(self, tableIndex, partialNewState):
        """
        Work routine to solve for Ye in neutrino-less beta equilibrium.
        Expects 'None' in partialNewState list where Ye needs filling in
        Modifies self.physicalState
        """
        munu = self.h5file['munu']
        ye = self.h5file['ye']

        previousPoint = tuple( 0 if val is None else val for val in tableIndex )
        currentMunu = previousMunu = munu[previousPoint]
        gotZero = False
        i = 0
        for i in range(len(munu)):
            thisPoint = []
            for j in tableIndex:
                if j is None:
                    thisPoint.append(i)
                else:
                    thisPoint.append(j)
            thisPoint = tuple(thisPoint)
            #adjust physical state to current state
            self.physicalState = tuple( ye[i] if val is None else val for val in partialNewState )
            currentMunu = self.interpolateTable(thisPoint, 'munu')
            #print thisPoint, currentMunu , self.physicalState
            if currentMunu * previousMunu < 0.0:
                gotZero = True
                break
            previousPoint = thisPoint
            previousMunu = currentMunu
        assert gotZero, \
            "Did not find zero of munu for all ye at non-ye parameters: %s" % partialNewState

        index = i - 1
        deltaMunu = -previousMunu / (currentMunu - previousMunu)
        return ye[index] * (1. - deltaMunu) + ye[index + 1] * deltaMunu

    def lookupIndex(self, indVar, value):
        """
        Returns the index directly preceding value
        Uses a dumb sequential search to find index.
        """
        assert indVar in self.indVars

        if indVar in self.logVars:
            value = math.log10(value)
            indVar = 'log' + indVar

        previousValue = -1.0e300
        thisVal = None
        for i, thisVal in enumerate(self.h5file[indVar][:]):
            assert thisVal > previousValue, \
                "Lookup index assumes independent variable table is increasing!"
            #print thisVal
            if value < thisVal:
                assert i > 0, "Uh oh, value %s for variable '%s' is below the table "\
                              "minimum: %s" % (value, indVar, thisVal)
                return i - 1
            previousValue = thisVal
        assert thisVal is not None, "Looks like table for %s is empty!" % indVar
        assert False, "Uh oh, value %s for variable '%s' is above the table "\
                      "maximum: %s" % (value, indVar, thisVal)

    #just does trilinear interpolation; does NOT try and account/correct
    # for log spacing in independent or dependent variable
    # got lazy; am not generalizing, copying from
    # http://en.wikipedia.org/wiki/Trilinear_interpolation
    # Note: index lookup is separated from interpolation routine
    #       so that different interpolators may be written and plugged in
    def interpolateTable(self, tableIndex, quantity):
        """
        Given table indexes preceding location of independent variables
        (physical state), does trilinear interpolation of 'quantity' to the
         current physical state.
        """
        tableIndex = tuple(tableIndex)
        y = self.h5file[quantity]

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

        return vector[0] * (1. - dxs[2]) + vector[1] * dxs[2]

        # print xs
        # print x0
        # print x1
        # print numpy.shape(self.h5file[quantity])
        # print self.h5file[quantity][tableIndex]

