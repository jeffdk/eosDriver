#!/opt/local/bin/python

import sys
from units import *
from numpy import linspace, zeros, log10, pi, sqrt, exp
from eosDriver import eosDriver
import makeeostable


# global constants

myeos = eosDriver('LS220_234r_136t_50y_analmu_20091212_SVNr26.h5')

kb_erg = 1.380658e-16
avo = 6.0221367e23
hc_mevcm = 1.97326966e-11*2.0*pi
pi = 3.14159265358979e0
mev_to_erg = 1.60217733e-6

def get_ynu(eta,rho,temp):
    
    nnu = 4*pi*(temp/hc_mevcm)**3 * 1.0/3.0 * eta * (eta**2 + pi**2)
    xynu = nnu/(rho*avo)

    return xynu

def get_pnu(eta,temp,rho,rhotrap):
    pnu =  mev_to_erg * 4.0*pi/3.0 * (temp**4/hc_mevcm**3) * \
        (21.0*pi**4 / 60.0 + 0.5*eta**2 * \
             (pi**2 + 0.5*eta**2)) * exp(-rhotrap/rho)
    return pnu

rho = 1.0e14
ylep = 0.1

temp = 0.1
ylep = myeos.setBetaEqState({'rho': rho,
                           'temp': temp})
print "initial Y_lep", ylep
ye = ylep
temp = 30.0

count = 0
prec = 1.0e-6
F = 1.0

minye = min(myeos.h5file['ye'])*1.01
maxye = max(myeos.h5file['ye'])*0.99

# get two initial guesses for ynu
ye1 = ye
ye2 = ye*1.01
myeos.setState({'rho': rho,
                'ye': ye1,
                'temp': temp})
munu1 = myeos.query('munu')
eta1 = munu1 / temp
ynu1 = get_ynu(eta1,rho,temp)

myeos.setState({'rho': rho,
                'ye': ye2,
                'temp': temp})
munu2 = myeos.query('munu')
eta2 = munu2 / temp
ynu2 = get_ynu(eta2,rho,temp)

dynudye = (ynu2-ynu1)/(ye2-ye1)
tiny = 1.0e-10
fac = 1.0
while count < 1000 and abs(F) > prec:
    #if(count>1000):
    #    fac = 0.1
    dF = -dynudye - 1.0
    ye2 = min(max(ye1 - F/dF*fac,minye),maxye)
    myeos.setState({'rho': rho,
                    'ye': ye2,
                    'temp': temp})
    munu2 = myeos.query('munu')
    eta = munu2/temp
    ynu2 = get_ynu(eta,rho,temp)
    count += 1
    F = ylep - ye2 - ynu2
    dynudye = (ynu2-ynu1)/(ye2-ye1+tiny)
    ye1 = ye2
    ynu1 = ynu2
    print "%5d %15.6E %15.6E %15.6E %15.6E %15.6E" % (count,munu2,ylep,
                                                      ynu1,ye1,abs(F))

ye = ye1

print count,ye,F

myeos.setState({'rho': rho,
                'ye': ye,
                'temp': temp})

press = 10.0e0**myeos.query('logpress')

print eta, press, get_pnu(eta,temp,rho,10.0**12.5)



sys.exit()

myeos.setState({'rho': rho,
                'ye': ye,
                'temp': temp})

ent = myeos.query('entropy')
print ent


ent = 3.0
print myeos.getTemperatureFromQuantityTYe({'rho': rho,
                                     'ye': ye,
                                     'temp': temp}, 'entropy', ent)




sys.exit()


rhomin = 2.0e3
rhomax = 1.0e12
nrhos = 100

logrhos = linspace(log10(rhomin),log10(rhomax),nrhos)
eostable = zeros((nrhos,2))
energy_shift = 0.0



rho = 1.000139999E+15
temp = 10.0
ye = 0.10

myeos.setState({'rho': rho,
                'ye': ye,
                'temp': temp})

press = myeos.query('logpress')

print 10.0**press

sys.exit()
    
epsdata = zeros(nrhos)

for i in range(nrhos):
    myeos.setState({'rho': 10.0**logrhos[i],
                    'ye': ye,
                    'temp': temp})
    
    press = myeos.query('logpress')
    
    myeos.setState({'rho': 10.0**logrhos[i],
                    'ye': ye,
                    'temp': temp})
    
    epsdata[i] = myeos.query('logenergy')
            

outfile = open("eosout.dat","w")
for i in range(nrhos):
    bstring = "%15.6E %15.6E\n" % (logrhos[i],epsdata[i])
    outfile.write(bstring)

outfile.close()
