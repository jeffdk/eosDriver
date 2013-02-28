
import plot_defaults
from eos import eos
import numpy
import matplotlib.pyplot as plt
from utils import lookupIndexBisect, linInterp, solveRootBisect
import matplotlib.pyplot as mpl

ls220 = eos('/home/jeff/work/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5')
shen = eos('/home/jeff/work/HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5')

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
getLogRhoIndex = lambda lr: lookupIndexBisect(lr, ls220.h5file['logrho'][:])
def entropyOfT(logtemp, logrho, target):

    lrindex = getLogRhoIndex(logrho)
    return linInterp(logtemp,ls220.h5file['logtemp'][:],
                     ls220.h5file['entropy'][11, :, lrindex]) - target


logtemp = numpy.log10(30.)

ye = 0.15

logrhos = numpy.arange(14.0,15.5,0.05)

#print logrhos

tindex = lookupIndexBisect(logtemp, ls220.h5file['logtemp'][:])
yeindex = lookupIndexBisect(ye,  ls220.h5file['ye'][:])
#print "ye index:  ", yeindex
#print ls220.h5file['ye'][yeindex+1]

scheck=[]
rollT =[]
esses=[1.0,1.46,2.0,3.0]
legends = []
for s in esses:
    ts=[]
    for lr in logrhos:
        sOfT = lambda T: entropyOfT(T, lr, s)
        getT = solveRootBisect(sOfT,-2.0,2.4)
        ts.append(numpy.power(10.0,getT))
        ls220.setState({'rho': numpy.power(10.0,lr), 'ye': 0.15, 'temp': numpy.power(10.0,getT)})
        scheck.append(ls220.query('entropy'))
        if s ==1.0:
            rollT.append(30/2.0*(1+numpy.tanh((lr-14.0)/0.5)) +0.5)
    legends.append("s = " + str(s))
    #print ts
    #print logrhos
    mpl.plot(numpy.power(10.0,logrhos),ts)
#print scheck
#print rollT
lg = mpl.legend(legends, loc=2)
lg.draw_frame(False)
mpl.plot(numpy.power(10.0,logrhos),rollT,ls='--')
mpl.xlabel(r"$\rho_b$ CGS")
mpl.ylabel(r"T($\rho_b$) ")
#mpl.title("LS220")
mpl.show()




print canonicalEntropy
