import datetime
import plot_defaults
from eosDriver import eosDriver
import numpy
import matplotlib.pyplot as plt
from utils import lookupIndexBisect, linInterp, solveRootBisect, multidimInterp
import matplotlib.pyplot as mpl

startTime = datetime.datetime.now()
ls220 = eosDriver('/home/jeff/work/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5')
shen = eosDriver('/home/jeff/work/HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5')

theEos = ls220

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
getLogRhoIndex = lambda lr: lookupIndexBisect(lr, theEos.h5file['logrho'][:])
def entropyOfT(logtemp, logrho, target):

    lrindex = getLogRhoIndex(logrho)
    return multidimInterp((logtemp, logrho), [theEos.h5file['logtemp'][:],
                                              theEos.h5file['logrho'][:]],
                          theEos.h5file['entropy'][11, :, :],
                          linInterp, 2) - target


logtemp = numpy.log10(30.)

ye = 0.15
tcold = 0.5

logrhos = numpy.arange(13.0,15.5,0.1)
rhos = numpy.power(10.0,logrhos)
#print logrhos

tindex = lookupIndexBisect(logtemp, theEos.h5file['logtemp'][:])
yeindex = lookupIndexBisect(ye,  theEos.h5file['ye'][:])
#print "ye index:  ", yeindex
#print ls220.h5file['ye'][yeindex+1]

scheck=[]
rollTs =[]
esses=[1.0,2.0,3.0]
escolors=['r','g','b']
proll = []
prollBeta = []
legends = []
for i, s in enumerate(esses):
    ts=[]
    tsBeta=[]
    lps =[]
    lpsBeta = []
    for lr in logrhos:
        rollT = (30/2.0*(1+numpy.tanh((lr-14.0)/0.5)) +tcold)
        theEos.setState({'rho': numpy.power(10.0,lr), 'ye': 0.15,
                         'temp': tcold})
        pCold = theEos.query('logpress')
        theEos.setBetaEqState({'rho': numpy.power(10.0,lr),
                               'temp': tcold})
        pBetaCold = theEos.query('logpress')
        if s ==1.0:
            rollTs.append(rollT)
            theEos.setState({'rho': numpy.power(10.0,lr), 'ye': 0.15,
                             'temp': rollT})
            proll.append(theEos.query('logpress')/pCold)
            theEos.setBetaEqState({'rho': numpy.power(10.0,lr),
                                   'temp': rollT})
            prollBeta.append(theEos.query('logpress')/pBetaCold)
        sOfT = lambda T: entropyOfT(T, lr, s)
        getT = solveRootBisect(sOfT,-2.0,2.4)
        ts.append(numpy.power(10.0,getT))
        unused, tbeta = theEos.setConstQuantityAndBetaEqState({'rho': numpy.power(10.0,lr)}, 'entropy', s)
        tsBeta.append(tbeta)
        lpsBeta.append(theEos.query('logpress')/pBetaCold)
        theEos.setState({'rho': numpy.power(10.0,lr), 'ye': 0.15, 'temp': numpy.power(10.0,getT)})
        lps.append(theEos.query('logpress')/pCold)

    legends.append("s = " + str(s) + "  ye = " + str(ye))
    print ts
    print logrhos
    mpl.plot(logrhos,lps, c=escolors[i])
    legends.append("s = " + str(s) + "  BetaEq")
    #legends.append(None)
    mpl.plot(logrhos, lpsBeta, ls='-.', c=escolors[i])
#print scheck
#print rollT
legends.append("IsoT + roll  ye = " + str(ye))
legends.append("IsoT + roll  ye = BetaEq")
mpl.plot(logrhos,proll, ls='-', c='m')
mpl.plot(logrhos,prollBeta,ls='-.',c='m')
lg = mpl.legend(legends, loc=1)
lg.draw_frame(False)
mpl.xlabel(r"$Log10(\rho_b$ CGS)")
mpl.ylabel(r"$P^{hot}/P^{cold}_{T=0.5MeV}$")
#mpl.title("LS220")

print "TIME DIFFERENCE: ", datetime.datetime.now() - startTime
mpl.show()
