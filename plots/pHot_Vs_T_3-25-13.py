
import plot_defaults
from eosDriver import eosDriver
import numpy
from matplotlib import rcParams
import matplotlib.pyplot as plt
rcParams['legend.fontsize'] = 16
ls220 = eosDriver('/home/jeff/work/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5')
shen = eosDriver('/home/jeff/work/HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5')

theEos = ls220

print ls220.h5file['ye'][:]
#print ls220.h5file['logrho'][:]
#print ls220.h5file['logtemp'][:]
#exit()

logTmin = theEos.h5file['logtemp'][0]
Tmin = numpy.power(10.0,logTmin)
logtemps = numpy.arange(logTmin, 2.5, .01)
#logtemps = ls220.h5file['logtemp'][:]
#logtemps = logtemps[:-1]
temps = numpy.power(10.0,logtemps)
print temps

def getP(rho,temp,ye):
    #print rho, temp, ye
    if ye > 0:
        theEos.setState({'rho': rho, 'ye': ye, 'temp': temp})
    else:
        theEos.setBetaEqState({'rho': rho, 'temp': temp})
    answer = numpy.power(10.0,theEos.query('logpress'))
#    if rho == 5.0e13 and temp < 4.0:
#        print rho,  temp, ye, answer
    return answer

print getP(1e14,1.0,0.2)


rhos = [1.0e13, 5.0e13, 1.0e14, 5.0e14, 1.0e15]

ye = 0.1
legend=[]
for rho in rhos:
    pOfTmin = getP(rho, Tmin, ye)
    oneDFunc = lambda T : getP(rho,T,ye)
    vectorFunc = numpy.frompyfunc(oneDFunc,1,1)
    legend.append(r"$\rho_b$ = " + str(rho))
    plt.loglog(temps, vectorFunc(temps))
plt.loglog(temps, 10e29*(10+temps*temps), '--')
legend.append(r"$10^{29} (10 + T^2)$")
#plt.loglog(temps, 10e28*(10+numpy.power(temps, 3)), '--')
#legend.append(r"$10^{28} (10 + T^3)$")
plt.loglog(temps, 10e26*(10+numpy.power(temps, 4)), '--')
legend.append(r"$10^{26} (10 + T^4)$")
lg = plt.legend(legend, loc=4)
lg.draw_frame(False)
plt.xlabel("T (MeV)")
plt.ylabel(r"${P(T)}$")
if ye < 0:
    ye="~BetaEq"
plt.title("LS220, Ye=" + str(ye))
plt.xlim([1,100])
plt.show()