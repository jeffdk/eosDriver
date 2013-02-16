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

        self.h5file = h5py.File(tableFilename, 'r')
        print self.h5file.keys()
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
        """
        self.physicalState = (None for unused in self.indVars)

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

    def lookupIndex(self, indVar, value):
        """
        Returns the index directly preceding value
        """
        assert indVar in self.indVars

        if indVar in self.logVars:
            value = math.log10(value)
            indVar = 'log' + indVar

        lastValue = -1.0e300
        thisVal = None
        for i, thisVal in enumerate(self.h5file[indVar][:]):
            assert thisVal > lastValue, \
                "Lookup index assumes independent variable table is increasing!"
            #print thisVal
            if value < thisVal:
                assert i > 0, "Uh oh, value %s for variable '%s' is below the table "\
                              "minimum: %s" % (value, indVar, thisVal)
                return i - 1
            lastValue = thisVal
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
        #print tableIndex, quantity
        tableIndex = tuple(tableIndex)
        y = self.h5file[quantity]


        xs = [] # xs is vector x, y, z from wikipeida
        x0 = [] # x0 is vector x0, y0, z0 from wikipedia
        x1 = [] # x0 is vector x1, y1, z1 from wikipedia
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

        dxs = (xs - x0)/(x1-x0)
        #print dxs

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

