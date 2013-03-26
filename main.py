from consts import CGS_C
from eosDriver import eosDriver, getTRollFunc, kentaDataTofLogRhoFit1
import numpy
import matplotlib.pyplot as mpl
from utils import lookupIndexBisect, linInterp, solveRootBisect, multidimInterp
import plot_defaults

shen = eosDriver('/home/jeff/work/HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5')

ls220 = eosDriver('/home/jeff/work/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5')

point = {'rho': 1e14, 'ye': .1, 'temp': 0.5}


ls220.setState(point)
print ls220.query(['logpress','entropy','logenergy'])

max = 30.0
min = 0.01

logrhos = numpy.arange(10.0,16.0,0.05)
rhos = numpy.power(10.0,logrhos)

midsAndScales=[(14.0, 0.5,), (13.5, 0.5), (14.0, 0.25)]
#ls220.writeRotNSeosfile("test2.eos", {'funcTofLogRho': 'kentaDataTofLogRhoFit2'}, 0.15)
labels = []

for mid, scale in midsAndScales:
    mpl.plot( logrhos,  getTRollFunc(max,min, mid, scale)(logrhos),
              logrhos,  kentaDataTofLogRhoFit1()(logrhos))
    labels.append("mid=" + str(mid) + " scale=" + str(scale))
mpl.legend(labels, loc =2)
mpl.ylabel("T (MeV)")
mpl.xlabel(r"log10($\rho_b$ CGS)")
#mpl.show()

#print ls220.solveForQuantity({'rho': 1e7, 'temp': 1.0}, 'munu', 0., bounds=None)
#print ls220.solveForQuantity({'rho': 1e15, 'ye': 0.1}, 'entropy', 1., bounds=None)

tempFunc = lambda x: numpy.log10(kentaDataTofLogRhoFit1()(x))
print ls220.solveForQuantity({'logtemp': 1., 'ye': 0.1}, 'logenergy',
                             numpy.log10(.1 * CGS_C**2 - ls220.energy_shift),
                             bounds=(11.,16.),
                             pointAsFunctionOfSolveVar={'logtemp': tempFunc,
                                                        'ye': None})

exit()
#print ls220.findIndVarOfMinAbsQuantity('ye', (0.0, 7), 'munu')
#exit()

print
print "---------------"
print

#ls220.setConstQuantityAndBetaEqState({'rho': 1e14}, 'entropy', 1.0)

#ls220.writeRotNSeosfile('EOSisoentropic.dat', {'quantity': 'entropy', 'target': 1.0}, ye=0.15)

print ls220.setConstQuantityAndBetaEqState({'rho': 1e8}, 'entropy', 25.0 )
print ls220.query(['entropy','logpress','logenergy'])
ls220.writeRotNSeosfile('EOSisothermalBetaEq.dat', {'T': 30.0,
                                              'rollMid': 14.0,
                                              'rollScale': 0.5,
                                              'eosTmin': 0.5}, ye=None)
print "--------Dur----------"
print ls220.setBetaEqState({'rho': 1e7, 'temp': 1.0})

exit()

def func(rho):

    global ls220
    global shen
    #global ye
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
    #global ye
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