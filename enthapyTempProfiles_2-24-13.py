
import plot_defaults
from eos import eos
import numpy
import matplotlib.pyplot as plt


ls220 = eos('/home/jeff/work/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5')
shen = eos('/home/jeff/work/HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5')

theEos = ls220

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

def getS(rho,temp,ye):
    #print rho, temp, ye
    if ye > 0:
        theEos.setState({'rho': rho, 'ye': ye, 'temp': temp})
    else:
        theEos.setBetaEqState({'rho': rho, 'temp': temp})
    answer = theEos.query('entropy')
#    if rho == 5.0e13 and temp < 4.0:
#        print rho,  temp, ye, answer
    return answer

#pick initial entropy

theEos.setState({'rho': 1.0e15, 'ye': 0.15, 'temp': 30.0})
canonicalEntropy = theEos.query('entropy')
print canonicalEntropy
print theEos.setBetaEqState({'rho': 1.0e15,  'temp': 30.0})
canonicalEntropy = theEos.query('entropy')
print canonicalEntropy
canonicalEntropy = 0.9
fixedYe = 0.15




print canonicalEntropy
