import datetime
import plot_defaults
from eosDriver import eosDriver, getTRollFunc
import numpy
from utils import lookupIndexBisect, linInterp, solveRootBisect, multidimInterp
import matplotlib.pyplot as mpl

startTime = datetime.datetime.now()
ls220 = eosDriver('/home/jeff/work/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5')
shen = eosDriver('/home/jeff/work/HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5')

theEos = ls220
tableName = "ls220"
mpl.rcParams['figure.subplot.bottom'] = 0.16
mpl.minorticks_on()
theEos.setState({'rho': 1.0e15, 'ye': 0.15, 'temp': 30.0})
canonicalEntropy = theEos.query('entropy')
print canonicalEntropy
print theEos.setBetaEqState({'rho': 1.0e15,  'temp': 30.0})
canonicalEntropy = theEos.query('entropy')
print canonicalEntropy
canonicalEntropy = 0.9
fixedYe = 0.15


print
print "---------------"
print
ye = 0.15
tcold = 0.01

logrhos = numpy.arange(12.0, 15.5, 0.05)
rhos = numpy.power(10.0,logrhos)
#print logrhos

def entropyOfT(logtemp, logrho, target):


    return multidimInterp((ye, logtemp, logrho), [theEos.h5file['ye'][:],
                                                  theEos.h5file['logtemp'][:],
                                                  theEos.h5file['logrho'][:]],
                                                theEos.h5file['entropy'][:, :, :],
                                                linInterp, 2) - target




esses=[0.5, 1.0, 2.0, 3.0]
tRollFunc = getTRollFunc(40.0, tcold, 14.0, 0.5)
escolors=['r', 'g', 'b', 'k', 'c', 'm', 'y', 'orange', 'violet']
proll = []
prollBeta = []
legends = []
lrsList = []
pBetaList = []
for i, s in enumerate(esses):
    ts=[]
    tsBeta=[]
    lps =[]
    lpsBeta = []
    for lr in logrhos:
        rho = numpy.power(10.0, lr)
        theEos.setState({'rho': rho, 'ye': ye, 'temp': tcold})
        pCold = theEos.query('logpress')
        theEos.setBetaEqState({'rho': rho, 'temp': tcold})
        pBetaCold = theEos.query('logpress')

        sOfT = lambda T: entropyOfT(T, lr, s)
        getT = solveRootBisect(sOfT, -2.0, 2.4)
        thisT = numpy.power(10.0, getT)
        ts.append(thisT)
        unused, tbeta = theEos.setConstQuantityAndBetaEqState({'rho': rho}, 'entropy', s)
        tsBeta.append(tbeta)
        lpsBeta.append(theEos.query('logpress') / pBetaCold - 1.)
        theEos.setState({'rho': rho, 'ye': ye, 'temp': thisT})
        #print thisT, tbeta, tcold
        lps.append(theEos.query('logpress') / pCold - 1.)

    pBetaList.append(lpsBeta)
    lrsList.append(logrhos)
    legends.append("s = " + str(s))
    mpl.semilogy(logrhos,lps, c=escolors[i])
    #legends.append("s = " + str(s) + "  BetaEq")
    #legends.append(None)
    #mpl.plot(logrhos, lpsBeta, ls='-.', c=escolors[i])

tempRolls = [(40.0, 14.0, 0.5,), (20.0, 14.0, 0.5),  (40.0, 13.5, 0.5), (40.0, 14.0, 0.25)]
for i, troll in enumerate(tempRolls):
    lpsBeta = []
    lps = []
    tRollFunc = getTRollFunc(troll[0], tcold, troll[1], troll[2])
    for lr in logrhos:
        rho = numpy.power(10.0, lr)
        theEos.setState({'rho': rho, 'ye': ye, 'temp': tcold})
        pCold = theEos.query('logpress')
        theEos.setBetaEqState({'rho': rho, 'temp': tcold})
        pBetaCold = theEos.query('logpress')

        thisT = tRollFunc(lr)
        ts.append(thisT)
        theEos.setBetaEqState({'rho': rho, 'temp': thisT})
        lpsBeta.append(theEos.query('logpress') / pBetaCold - 1.)

        theEos.setState({'rho': rho, 'ye': ye, 'temp': thisT})
        lps.append(theEos.query('logpress') / pCold - 1.)

    pBetaList.append(lpsBeta)
    lrsList.append(logrhos)
    legends.append(str(troll).replace(" ", ""))
    mpl.semilogy(logrhos, lps, c=escolors[i + len(esses)])
mpl.rcParams['legend.labelspacing'] = -.10
mpl.rcParams['legend.labelspacing'] = -.10
lg = mpl.legend(legends, loc=4)
lg.draw_frame(False)
mpl.xlabel(r"$\mathrm{log10}(\rho_b$ CGS)", labelpad=12)
#mpl.ylabel(r"$P^{hot}/P^{cold}_{T=" + str(tcold) + "MeV} - 1$")
mpl.ylabel(r"$\mathrm{log10}(P_\mathrm{hot}/P_\mathrm{cold} - 1)$", labelpad=12)
lg = mpl.legend(legends, loc=4)
lg.draw_frame(False)
#mpl.title("LS220")
for i in range(len(pBetaList)):
    mpl.semilogy(lrsList[i], pBetaList[i], ls='--', c=escolors[i])
mpl.axes().annotate(tableName, xy=(12.2, 4.0), fontsize=20 )
ylocs, ylabels = mpl.yticks()
mpl.yticks(ylocs, map(lambda x: "%.1f" % x, numpy.log10(ylocs)))
print "TIME DIFFERENCE: ", datetime.datetime.now() - startTime
mpl.ylim([10**-6, 16])
mpl.xlim([12.0, 15.5])
mpl.show()