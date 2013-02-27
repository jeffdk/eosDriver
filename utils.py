"""
Defines interpolation and other routines.
Interpolators should take:
 x,  xs,  ys
 Where
"""

def solveRootBisect(func, x0, x1, relTol=1.0e-8, maxIterations=30):
    """
    Finds a root of func between x0 and x1
    to relative tolerance relTol
    """
    assert func(x0) * func(x1) <= 0., "Root not bracketed by x0 and x1!"
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
            #print "diff is 1"
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
        "x: %s \t is less than table minimum, %s" % (x,tableXmin)
    assert x <= tableXmax, \
        "x: %s \t is greater than table minimum, %s" % (x,tableXmin)
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