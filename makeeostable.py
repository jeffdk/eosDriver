import sys
from units import *
from scipy import *
from eosDriver import eosDriver
from utils import lookupIndexBisect, linInterp, solveRootBisect

myeos = eosDriver('LS220_234r_136t_50y_analmu_20091212_SVNr26.h5')

def makeeostable(nrhos,rhomin,rhomax,myeos,mytype):

    logrhos = linspace(log10(rhomin),log10(rhomax),nrhos)
    eostable = zeros((nrhos,2))
    energy_shift = 0.0

    if(mytype == 'fixed_ye_temp'):

        ye = 0.15
        temp = 0.5
        
        energy_shift = myeos.h5file['energy_shift'][0]
        
        for i in range(nrhos):
            myeos.setState({'rho': 10.0**logrhos[i],\
                                'ye': ye,\
                                'temp': temp})
            
            press = myeos.query('logpress')

            myeos.setState({'rho': 10.0**logrhos[i],\
                                'ye': ye,\
                                'temp': temp})

            eps = myeos.query('logenergy')
            
            # convert units
            eostable[i,0] = log10(10.0**press * press_gf)
            eostable[i,1] = log10(10.0**eps * eps_gf)

        energy_shift = energy_shift*eps_gf

    elif(mytype == 'poly_G2_K100'):
        eostable[:,0] = log10(100.0*(10.0**logrhos[:]*rho_gf)**2.0)
        eostable[:,1] = log10(10.0**eostable[:,0]/(10.0**logrhos[:]*rho_gf) )

    else:
        print "This kind of table can't be done yet: ",mytype
    

    logrhos = log10(10.0**logrhos * rho_gf)
    dlrho = logrhos[1]-logrhos[0]
    return (eostable,energy_shift,dlrho,logrhos)

