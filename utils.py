"""
Defines interpolation and other routines.
Interpolators should take:
 x,  xs,  ys
 Where
"""
import numpy

def relativeError(a,b):
    """
    Calculates relative error between a and b as shown.
    Result is always positive.
    """
    assert not a == b == 0.0, "Relative error between 0 and 0 is undefined!"
    return abs( (a - b) / (abs(a)+abs(b)) )

#Ideally, interpolator should be an object that knows its
# passing interpolation order lets us recursively call
# multidimInterp without interpolating every point in the table
def multidimInterp(point, tablePoints, tableData, interpolator,  interpolationOrder=2):
    assert len(point) == tableData.ndim
    assert len(tablePoints) == tableData.ndim

    dim = tableData.ndim
    #
    offset = int(interpolationOrder / 2) - 1
    indexesStart = [lookupIndexBisect(point[i], tablePoints[i]) - offset
                    for i in range(dim)]
    indexesEnd = [i + interpolationOrder -1 for i in indexesStart]

    localTablePoints = [None for _ in tablePoints]
    for i in range(dim):
        indexesThisAxis = [j for j in range(indexesStart[i], indexesEnd[i] + 1 )]
        #print indexesThisAxis
        tableData = numpy.take(tableData, indexesThisAxis, axis=i)
        localTablePoints[i] = numpy.take(tablePoints[i], indexesThisAxis)

    answer = multidimInterpWork(point, localTablePoints, tableData, interpolator)


    return answer

# multidimInterp which handles reduction  of the table to something
#  which is no larger than needed by the interpolation stencil
# the real interpolation takes place in the function below
def multidimInterpWork(point, tablePoints, tableData, interpolator):
    assert isinstance(tableData, numpy.ndarray), "Input table must be a numpy array"
    assert len(point) == tableData.ndim, \
        "Point interpolating to must have same dimension as data table \n " \
        "PointDim: %s \t tableDataDim: %s " % (len(point),tableData.ndim)
    # print "Point: ", point
    # print "tablePoints: ", tablePoints
    # print "tableData shpe: ", numpy.size(tableData)
    # print "tableData: ", tableData
    interpolationOrder = numpy.shape(tableData)[0]

    if tableData.ndim == 1:
        return interpolator(point, tablePoints[0], tableData)[0]
    else:
        #print point[1:], tablePoints[1:], tableData[1:]
        npoints = numpy.size(tableData[1:])
        reducedTableData = numpy.zeros(npoints)
        #thisPoint = []
        for j in range(npoints):
            thisPoint = []
            for i in range(interpolationOrder):
                thisPoint.append(tableData[i, ...].flatten()[j] )
            reducedTableData[j] = interpolator(point[0],tablePoints[0],thisPoint)
        reducedTableData = reducedTableData.reshape(numpy.shape(tableData[1, ...]))
        return multidimInterpWork(point[1:], tablePoints[1:], reducedTableData, interpolator )

    pass

#Scipy root solvers throw value error when bracketing problem,
# so derive BracketingError from Value error
class BracketingError(ValueError):
    def __init__(self,x0,x1,fx0,fx1,func):
        self.message = " \n Root not bracketed between x0 = %s  x1 = %s" % (x0,x1)
        self.message += "\n f(x0) = %s   f(x1) = %s" % (fx0,fx1)
        self.message += "\n for function %s " % func

    def __str__(self):
        return self.message

def solveRootBisect(func, x0, x1, optionalParams=(), relTol=1.0e-8, maxIterations=40):
    """
    Finds a root of func between x0 and x1
    to relative tolerance relTol
    Optional parameters is not implemented yet but included so the function
    arguments match that of scipy's root solvers.
    """
    if func(x0) * func(x1) > 0.:
        raise BracketingError(x0,x1,func(x0),func(x1),func)
    if func(x0) == 0.:
        return x0
    if func(x1) == 0.:
        return x1

    upper = x1
    lower = x0

    foundRoot = False

    iteration = 0

    while not foundRoot:
        nextX= (upper - lower)/2. + lower
        fOfNextX = func(nextX)
        if func(upper) * fOfNextX < 0:
            lower = nextX
        elif func(lower) * fOfNextX  < 0:
            upper = nextX
        elif 0. == fOfNextX:
            return nextX



        if abs((upper - lower)/(lower + upper)) < relTol:
            foundRoot = True
            return lower
        iteration += 1
        assert iteration < maxIterations

#xs as in multiple values of x "x is one of the exes"
def linInterp(x, xs, ys):
    """
    Linear interpolation to point x from 1D table with independent
    variable 'xs' and dependent variable 'ys'.
    """
    assert len(xs) == len(ys), \
        "Length xs: %s \t Length ys: %s" % (len(xs),len(ys))

    tableXmin = xs[0]
    tableXmax = xs[-1]

    assert x >= tableXmin, \
        "x: %s \t is less than table minimum, %s" % (x, tableXmin)
    assert x <= tableXmax, \
        "x: %s \t is greater than table maximum, %s" % (x, tableXmax)
    i = lookupIndexBisect(x, xs)

    dx = (x - xs[i])/(xs[i+1] - xs[i])

    return ys[i+1] * dx + (1.- dx) * ys[i]

def lookupIndexBisect(x, xs):
    """
    Assuming a monotonically increasing list of values xs, returns
    the index directly preceding x.

    Gives error for values less than list minimum but happily
    returns maximum index for any value above the list maximum
    """
    assert x >= xs[0], "there is no index preceding x, %s, since it " \
                       "is below the table minimum: %s" % (x,xs[0])
    upper = len(xs) -1
    lower = 0

    maxIterations = upper
    foundIndex = False
    iteration = 0

    while not foundIndex:
        nextIndex = int((upper - lower)/2) + lower

        if x > xs[nextIndex]:
            lower = nextIndex
        elif x < xs[nextIndex]:
            upper = nextIndex
        elif x == xs[nextIndex]:
            return nextIndex-1

        if abs(upper - lower) == 1:
            #print "diff is 1"
            foundIndex = True
            return lower
        iteration += 1
        assert iteration < maxIterations