#!/opt/local/bin/python

import sys
from units import *
from numpy import linspace, zeros, log10, pi, sqrt
from eosDriver import eosDriver
import makeeostable
from tov import *

myeos = eosDriver('LS220_234r_136t_50y_analmu_20091212_SVNr26.h5')

# make fixed temperature, Y_e sequence
tovinfo = tovinfoclass()
tovinfo.polyK = 100.0
tovinfo.polyG = 2.0
tovinfo.nzones = 80000
tovinfo.rmax = 50.0
tovinfo.eostype = 3

mytype = "fixed_ye_entropy"
par1 = 0.5
par2 = 0.1
rhomin = 1.0e7
rhomax = 8.0e15
tovinfo.eoslrhomin = log10(rhomin*rho_gf)
tovinfo.eoslrhomax = log10(rhomax*rho_gf)
(tovinfo.eostable,tovinfo.eosepsshift,dlrho,tovinfo.logrhos) \
    = makeeostable.makeeostable(\
    tovinfo.nrhos,rhomin,rhomax,myeos,mytype,par1,par2)

tovinfo.eosdlrho = dlrho
tovinfo.eosdlrhoi = 1.0/dlrho
(tovinfo.minpress,bogus) = tabeos_press_eps(rho_gf*rhomin,tovinfo)


rho_c = 4.717563e-04
rho_c = 2.0e14*rho_gf
rho_c = 8.615461e-04
rho_c = 4.0e14*rho_gf
rho_c = 6.864997E-04
rho_c = 7.705537E-3
rho_c = 1.258022215E+15*rho_gf
print rho_c*inv_rho_gf
(tovout,isurf,rad,dr) = tov_integrate(rho_c,tovinfo)


print "%4d %15.6E %15.6E %15.6E %15.6E" % \
    (isurf,rho_c,tovout[isurf,3],tovout[isurf,4],rad[isurf]*inv_length_gf)

print rho_c*inv_rho_gf, tovout[0,0]*inv_rho_gf
outfile = open("pytov.dat","w")
for i in range(isurf):
    rline = "%18.9E %18.9E %18.9E %18.9E\n" % (rad[i]*inv_length_gf,\
                                                   tovout[i,0]*inv_rho_gf,\
                                                   tovout[i,1]*inv_press_gf,\
                                                   tovout[i,3])
    outfile.write(rline)

outfile.close()
