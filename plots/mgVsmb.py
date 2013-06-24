import matplotlib
from matplotlib import pyplot as plt
import numpy
import plot_defaults

def vline(x, **kwargs):
    x1, x2, y1, y2 = plt.axis()
    print y1,y2
    plt.plot([x, x], [y1, y2], **kwargs)

cols = {'index': 0,
        'rhob': 1,
        'mg': 2,
        'mb': 3,
        'r': 4}

eoss = ["HShen", "LS220"]

xVars = ['rhob']
#xLab = "$M_\mathrm{b}$"
xLab = r"$\rho_\mathrm{b,c}$"
yVars = ['mg', 'mb']
#yLab = "$M_\mathrm{g}$"
yLab = "$M_\mathrm{g}$/$M_\mathrm{b}$"
xFunc = lambda x: x
yFunc = lambda x, y:  x / y

for eos in eoss:

    models = numpy.loadtxt("/home/jeff/work/tov_sequence_%s_eostable_BetaEq_T=00.500.dat" % eos)
    models = models.transpose()
    arg_mbMax = numpy.argmax(models[cols['mb']])
    arg_mgMax = numpy.argmax(models[cols['mg']])
    mbMax = models[cols['rhob']][arg_mbMax]
    mgMax = models[cols['rhob']][arg_mgMax]
    xData = xFunc(*(models[cols[xVar]] for xVar in xVars))
    yData = yData = yFunc(*(models[cols[yVar]] for yVar in yVars))
    plt.plot(xData, yData, label=eos)
    #vline(mbMax, ls='--', dashes=plot_defaults.longDashes, label="$M_\mathrm{b,max}$")
    #vline(mgMax, ls='--', dashes=plot_defaults.goodDashDot, label="$M_\mathrm{g,max}$")
plt.xlabel(xLab)
plt.ylabel(yLab)
plt.legend(loc=1)

plt.show()

