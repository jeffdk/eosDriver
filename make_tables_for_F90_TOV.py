#!/opt/local/bin/python
import sys
from units import *
from numpy import linspace, zeros, log10, pi, sqrt
from eosDriver import eosDriver
import makeeostable
from tov import *

#define EOSs

EOSlist = [ ["ls220","LS220_234r_136t_50y_analmu_20091212_SVNr26.h5"], 
            ["shen", "HShenEOS_rho234_temp136_ye50_version2.0_20120706.h5"]   
          ]

yes = [0.1,0.15]
tmin = 0.5
tmax = 50.0
dtemp = 0.5
ntemp = int((tmax-tmin)/dtemp)+1
temps = zeros(ntemp)
for i in range(ntemp):
	temps[i] = 0.5 + dtemp*i

rhomin = 1.0e6
rhomax = 8.0e15
nrhos = 600

for ieos in range(len(EOSlist)):
    myeos = eosDriver(EOSlist[ieos][1])
    for ii in range(len(temps)):
        for jj in range(len(yes)):
            mytype = "fixed_ye_temp"
            par1 = temps[ii]
            par2 = yes[jj]
            print "T = %5.2f, Y_e = %5.2f" % (par1,par2)
            print "Preparing EOS table: ",mytype
            makeeostable.makeeostable(\
                nrhos,rhomin,rhomax,myeos,EOSlist[ieos][0],\
                    mytype,par1,par2)

    del myeos
    
    
    

