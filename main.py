from eos import eos
import numpy
import matplotlib.pyplot as mpl

shen = eos('/home/jeff/work/HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5')

ls220 = eos('/home/jeff/work/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5')

#print ls220.lookupIndex('rho',1e4)
#ls220.lookupIndex('temp',10.0)
ls220.lookupIndex('rho',10000.0)
ls220.setState({'rho': 1e14, 'ye': .1, 'temp': 0.5})
print ls220.query('logpress')

def func(rho):

    global ls220
    global shen
    global ye
    global temp
    theDict = {'rho': rho, 'ye': ye, 'temp': temp}
    print theDict
    shen.setState(theDict)
    return shen.query('logpress')

vfunc = numpy.frompyfunc(func,1,1)


logrhos = numpy.arange(13.0,15.5,0.1)
rhos = numpy.power(10.0,logrhos)

print logrhos
print rhos

ye = 0.1
temp = 0.5

xs = logrhos
yscold = numpy.power(10.0, vfunc(rhos))

temp = 20.0
yshot = numpy.power(10.0, vfunc(rhos))

yes = [0.05, 0.1, 0.2, 0.3]
legend = []
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