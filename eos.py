import h5py
import math
import numpy


class eos(object):

    valuesDict = None
    indVars = ('rho', 'ye', 'temp')  #indepenednt variables
    logVars = ('energy', 'press', 'rho', 'temp')

    h5file = None

    def __init__(self, tableFilename):

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


    def query(self, pointDict, quantity):
        answer = 0

        tableIndexes = []
        for indVar in self.indVars:
            assert indVar in pointDict, \
                "You have not specified an required independent variable in pointDict!"
            #print indVar
            value = pointDict[indVar]

            #print indVar, value
            tableIndexes.append(self.lookupIndex(indVar, value))


        return self.interpolateTable(pointDict, tableIndexes, quantity)

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
    def interpolateTable(self, pointDict, tableIndex, quantity):
        #print tableIndex, quantity
        tableIndex = tuple(tableIndex)
        y = self.h5file[quantity]

        xs = []
        x0 = []
        x1 = []
        for i, indVar in enumerate(self.indVars):
            value = pointDict[indVar]
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

