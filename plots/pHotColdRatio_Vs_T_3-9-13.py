
import plot_defaults
from eosDriver import eosDriver
import numpy
import matplotlib.pyplot as plt
import matplotlib.ticker as mplTicker

ls220 = eosDriver('/home/jeff/work/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5')
shen = eosDriver('/home/jeff/work/HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5')

theEos = ls220


logTmin = theEos.h5file['logtemp'][0]
Tmin = numpy.power(10.0,logTmin)
logtemps = numpy.arange(logTmin, 1.9, .01)
#logtemps = ls220.h5file['logtemp'][:]
#logtemps = logtemps[:-1]
temps = numpy.power(10.0,logtemps)
#print temps

def getP(rho,temp, ye):
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


rhos = [1.0e13, 1.0e14, 5.0e14, 1.0e15]
colors = ['r', 'b', 'g', 'm']
ye = 0.15
legend=[]
#minorLocator = mplTicker.MultipleLocator(2)
#ax = plt.axes()
#ax.xaxis.set_minor_locator(minorLocator)
plt.minorticks_on()
for i, rho in enumerate(rhos):
    pOfTmin = getP(rho, Tmin, ye)
    oneDFunc = lambda T : getP(rho,T,ye) / pOfTmin - 1.0
    vectorFunc = numpy.frompyfunc(oneDFunc,1,1)
    legend.append(r"$\rho_b$ = " + plot_defaults.fixScientifcNotation(rho))
    plt.semilogy(temps, vectorFunc(temps), c=colors[i])
lg = plt.legend(legend, loc=4)
#lg.draw_frame(False)

#labelpad is the spacing in points between the label and the x-axis
plt.xlabel("T (MeV)", labelpad=14)
plt.ylabel(r"$\mathrm{log10}({P(T)}/{P(T_\mathrm{min})}$ - 1)", labelpad=12)
if ye < 0:
    ye="~BetaEq"
#plt.title("Fractional increase in pressure, Shen, Ye=" + str(ye))

ye = 0
for i, rho in enumerate(rhos):
    pOfTmin = getP(rho, Tmin, ye)
    oneDFunc = lambda T : getP(rho,T,ye) / pOfTmin - 1.0
    vectorFunc = numpy.frompyfunc(oneDFunc,1,1)
    legend.append(r"$\rho_b$ = " + plot_defaults.fixScientifcNotation(rho))
    plt.semilogy(temps, vectorFunc(temps), c=colors[i], ls='--')

ylocs, ylabels = plt.yticks()
plt.yticks(ylocs, map(lambda x: "%.1f" % x, numpy.log10(ylocs)))
plt.ylim([10**-5, 10**3])
plt.xlim([0, 70])
plt.show()