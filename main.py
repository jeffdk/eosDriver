from eos import eos
import numpy
import matplotlib.pyplot as mpl
from utils import lookupIndexBisect, linInterp, solveRootBisect
import plot_defaults

shen = eos('/home/jeff/work/HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5')

ls220 = eos('/home/jeff/work/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5')

#print ls220.lookupIndex('rho',1e4)
#ls220.lookupIndex('temp',10.0)
ls220.lookupIndex('rho',10000.0)
ls220.setState({'rho': 1e14, 'ye': .1, 'temp': 0.5})
print ls220.query('logpress')

#print ls220.h5file['logtemp'][:]
interpVal = 3.1
i =  lookupIndexBisect(interpVal, ls220.h5file['logrho'][:])
print  i, ls220.h5file['logrho'][i],  ls220.h5file['logrho'][i + 1]

print i, ls220.h5file['logpress'][10,10,i], ls220.h5file['logpress'][10,10,i+1]
print linInterp(interpVal,ls220.h5file['logrho'][:],
                ls220.h5file['logpress'][10,10,:] )


print
print "---------------"
print
getLogRhoIndex = lambda lr: lookupIndexBisect(lr, ls220.h5file['logrho'][:])
def entropyOfT(logtemp, logrho):

    lrindex = getLogRhoIndex(logrho)
    return linInterp(logtemp,ls220.h5file['logtemp'][:],
                     ls220.h5file['entropy'][11, :, lrindex]) - 1.46


logtemp = numpy.log10(30.)

ye = 0.15

logrhos = numpy.arange(12.0,16,0.1)

print logrhos

tindex = lookupIndexBisect(logtemp, ls220.h5file['logtemp'][:])
yeindex = lookupIndexBisect(ye,  ls220.h5file['ye'][:])
print "ye index:  ", yeindex
print ls220.h5file['ye'][yeindex+1]
ts=[]
scheck=[]
rollT =[]
for lr in logrhos:
    sOfT = lambda T: entropyOfT(T, lr)
    getT = solveRootBisect(sOfT,-2.0,2.4)
    ts.append(numpy.power(10.0,getT))
    ls220.setState({'rho': numpy.power(10.0,lr), 'ye': 0.15, 'temp': numpy.power(10.0,getT)})
    scheck.append(ls220.query('entropy'))
    rollT.append(30/2.0*(1+numpy.tanh((lr-14.0)/0.5)) +0.5)

print scheck
print rollT
mpl.plot(logrhos,ts)
mpl.plot(logrhos,rollT,ls='--')
mpl.xlabel(r"Log10($\rho_b$ CGS)")
mpl.ylabel(r"T($\rho_b$)  $\,\,s = 1.46\, \frac{k_B}{baryon}$ ")
#mpl.title("LS220")
mpl.show()



exit()

def func(rho):

    global ls220
    global shen
    global ye
    global temp
    theDict = {'rho': rho, 'ye': ye, 'temp': temp}
    if ye > 0:
        shen.setState(theDict)
    else:
        theDict = {'rho': rho, 'temp': temp}
        shen.setBetaEqState(theDict)
    return shen.query('logpress')

vfunc = numpy.frompyfunc(func,1,1)


logrhos = numpy.arange(13.0,15.5,0.1)
rhos = numpy.power(10.0,logrhos)

print logrhos
print rhos
legend = ["Ye = BetaEq"]

ye =  -1.0 #beta equilibrum
temp = 0.5

xs = logrhos
yscold = numpy.power(10.0, vfunc(rhos))
print '--'
temp = 20.0
yshot = numpy.power(10.0, vfunc(rhos))
mpl.plot(xs,yshot/yscold)
yes = [0.05, 0.1, 0.2, 0.3]

for thisYe in yes:
    global ye
    ye = thisYe
    temp = 0.5
    yscold = numpy.power(10.0, vfunc(rhos))
    temp = 20.0
    yshot = numpy.power(10.0, vfunc(rhos))
    legend.append("Ye = " + str(ye))
    mpl.plot(xs,yshot/yscold)
mpl.grid()
mpl.legend(legend)
mpl.xlabel("Log10(rho_b CGS)")
mpl.ylabel("P_20MeV / P_0.5MeV")
mpl.title("Hot vs. Cold Pressure HShen")
mpl.show()