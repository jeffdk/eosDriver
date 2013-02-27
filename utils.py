"""
Defines interpolation and other routines.
Interpolators should take:
 x,  xs,  ys
 Where
"""


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

    assert x > tableXmin, \
        "x: %s \t is less than table minimum" % (x,tableXmin)
    assert x < tableXmax, \
        "x: %s \t is greater than table minimum" % (x,tableXmin)


def lookupIndexBisect(x, xs):
    """
    Assuming a monotonically increasing list of values xs, returns
    the index directly preceding x.
    """
    #todo: add out of bounds checks
    upper = len(xs)
    lower = 0

    maxIters = upper
    foundIndex = False
    iter = 0

    while not foundIndex:
        nextIndex = int((upper - lower)/2) + lower

        if x > xs[nextIndex]:
            lower = nextIndex
        elif x < xs[nextIndex]:
            upper = nextIndex

        if abs(upper - lower) == 1:
            print "diff is 1"
            foundIndex = True
            return lower
        iter += 1
        assert iter < maxIters