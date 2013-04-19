import datetime
import matplotlib
import plot_defaults
from eosDriver import eosDriver, getTRollFunc, kentaDataTofLogRhoFit1, kentaDataTofLogRhoFit2
import numpy
import matplotlib.pyplot as plt

startTime = datetime.datetime.now()
ls220 = eosDriver('/home/jeff/work/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5')
shen = eosDriver('/home/jeff/work/HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5')

theEos = shen
tableName = "HShen"
# matplotlib.rc('legend', fontsize=20)
# mpl.rcParams['figure.subplot.bottom'] = 0.16
# mpl.minorticks_on()
# mpl.rcParams['figure.subplot.right'] = 0.85

ye = 0.1
tcold = 0.01

##
#Grid for plotting points
logrhos = numpy.arange(11.0, 15.5, 0.05)
rhos = numpy.power(10.0, logrhos)

# Colors for each tempPrescription as is consistent with their coloring
colors = ['r', 'g', 'b', 'k', 'c']  # , 'm', 'y', 'orange', 'violet']
legends = ["c40p0", "c30p0", "c20p0", "c30p10", "c30p5"]
cXXp0_params = [(40.0, 14.18, 0.5,), (30.0, 14.055, 0.375),  (20.0, 13.93, 0.25)]
tempFuncs = [getTRollFunc(params[0], tcold, params[1], params[2]) for params in cXXp0_params]
tempFuncs.append(kentaDataTofLogRhoFit1())
tempFuncs.append(kentaDataTofLogRhoFit2())

###
# First calculate all the pressures we want
pBetaColds = []
pHotNuFull = [[] for unused in colors]
pHotNuLess = [[] for unused in colors]
pHotYeFixed = [[] for unused in colors]
pNuNuFull = [[] for unused in colors]
pNuNuLess = [[] for unused in colors]
pNuYeFixed = [[] for unused in colors]
for lr in logrhos:
    coldYe = theEos.setBetaEqState({'logrho': lr, 'temp': tcold})
    pBetaColds.append(theEos.query('logpress', deLog10Result=True))
    for i, tRollFunc in enumerate(tempFuncs):
        thisT = tRollFunc(lr)
        print thisT,

        nuFullYe = theEos.setNuFullBetaEqState({'logrho': lr, 'temp': thisT}, coldYe)
        pHotNuFull[i].append(theEos.query('logpress', deLog10Result=True))
        theEos.setState({'logrho': lr, 'temp': thisT, 'ye': nuFullYe})
        pNuNuFull[i].append(theEos.queryTrappedNuPress())
        pHotNuFull[i][-1] += pNuNuFull[i][-1] #Add neutrino pressure to hotNuFull

        nuLessYe = theEos.setBetaEqState({'logrho': lr, 'temp': thisT})
        pHotNuLess[i].append(theEos.query('logpress', deLog10Result=True))
        theEos.setState({'logrho': lr, 'temp': thisT, 'ye': nuLessYe})
        pNuNuLess[i].append(theEos.queryTrappedNuPress())

        theEos.setState({'logrho': lr, 'temp': thisT, 'ye': ye})
        pHotYeFixed[i].append(theEos.query('logpress', deLog10Result=True))
        theEos.setState({'logrho': lr, 'temp': thisT, 'ye': ye})
        pNuYeFixed[i].append(theEos.queryTrappedNuPress())

    print
pBetaColds = numpy.array(pBetaColds)
pHotNuFull = numpy.array(pHotNuFull)
pHotNuLess = numpy.array(pHotNuLess)
pHotYeFixed = numpy.array(pHotYeFixed)
pNuNuFull = numpy.array(pNuNuFull)
pNuNuLess = numpy.array(pNuNuLess)
pNuYeFixed = numpy.array(pNuYeFixed)
print pBetaColds

###
# Setup plot environment
#basics
myfig = plt.figure(figsize=(12, 22))
myfig.subplots_adjust(left=0.14)
myfig.subplots_adjust(bottom=0.08)
myfig.subplots_adjust(top=0.98)
myfig.subplots_adjust(right=0.97)
myfig.subplots_adjust(hspace=1.0e-10)


###
# First plot: P/P(cold_NuLess) - 1
plt.subplot(211)
plt.xlim([11.0, 15.3])

fracDiff = lambda a, b: a / b - 1.0
for i, color in enumerate(colors):
    plt.semilogy(logrhos, fracDiff(pHotNuFull[i], pBetaColds), c=color, label=legends[i])
    plt.semilogy(logrhos, fracDiff(pHotNuLess[i], pBetaColds), c=color, ls='-.')
    plt.semilogy(logrhos, fracDiff(pHotYeFixed[i], pBetaColds), c=color, ls='--')

plt.text(14.5, 1, tableName, fontsize=24)
plt.legend(loc=4)
plt.ylabel(r"$P_{`hot'}/P_{\nu\mathrm{-less, cold}} - 1$", labelpad=12)
plt.ylim([2e-4, 2e1])
# hide x-axis labeling of upper panels
ax = plt.gca()
for label in ax.get_xticklabels():
    label.set_visible(False)

###
# Second plot: P/P(cold_NuLess) - 1
plt.subplot(212)
plt.xlim([11.0, 15.3])

for i, color in enumerate(colors):
    plt.semilogy(logrhos, pNuNuFull[i]/pHotNuFull[i], c=color)
    plt.semilogy(logrhos, pNuNuLess[i]/pHotNuLess[i], c=color, ls='-.')
    plt.semilogy(logrhos, pNuYeFixed[i]/pHotYeFixed[i], c=color, ls='--')

plt.ylim([2e-5, 5e-1])
plt.ylabel(r"$P_{\nu}/P_\mathrm{total}$", labelpad=12)

plttxt = "$\\nu$-full $\\beta$-equil. \&  $P_\\nu$: --,\,\,\,\, $\\nu$-less: -.,\,\,\,\, $Y_e=0.1$: -- --"
plt.text(0.05, -0.2, plttxt, fontsize=26, horizontalalignment="left", transform=ax.transAxes)

plt.xlabel(r"$\mathrm{log10}(\rho_b$ [g/cm$^3$])", labelpad=12)

print "TIME DIFFERENCE: ", datetime.datetime.now() - startTime
plt.show()
exit()


ylocs, ylabels = plt.yticks()
plt.yticks(ylocs, map(lambda x: "%.1f" % x, numpy.log10(ylocs)))
