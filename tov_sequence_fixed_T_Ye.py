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
tovinfo.rmax = 100.0
tovinfo.eostype = 3

yes = [0.1,0.15,0.2,0.25,0.3]

tmin = 0.5
tmax = 50.0
dtemp = 0.5
ntemp = int((tmax-tmin)/dtemp)+1
temps = zeros(ntemp)
for i in range(ntemp):
	temps[i] = 0.5 + dtemp*i

for ii in range(len(temps)):
    for jj in range(len(yes)):
        mytype = "fixed_ye_temp"
        par1 = temps[ii]
        par2 = yes[jj]
        print "T = %5.2f, Y_e = %5.2f" % (par1,par2)
        print "Preparing EOS table: ",mytype
        rhomin = 1.0e6
        rhomax = 8.0e15
        tovinfo.eoslrhomin = log10(rhomin*rho_gf)
        tovinfo.eoslrhomax = log10(rhomax*rho_gf)
        (tovinfo.eostable,tovinfo.eosepsshift,dlrho,tovinfo.logrhos) \
            = makeeostable.makeeostable(\
            tovinfo.nrhos,rhomin,rhomax,myeos,mytype,par1,par2)
    
        tovinfo.eosdlrho = dlrho
        tovinfo.eosdlrhoi = 1.0/dlrho
        (tovinfo.minpress,bogus) = tabeos_press_eps(rho_gf*rhomin,tovinfo)

        outdata = tov_sequence(3.0e14,7.0e15,50,tovinfo)
    
        filename = mytype+"_T%06.3f_Ye%06.3f.dat" % (par1, par2)
        outfile=open(filename,"w")
        for i in range(len(outdata[:,0])):
            line = "%15.6E %15.6E\n" % (outdata[i,3]*inv_length_gf,outdata[i,1])
            outfile.write(line)
        outfile.close()

        # Look where baryonic mass is maximal, not gravitational
        # mass, because that could just be heat and not actual mass.
        # This has the advantage of excluding puffed-up low-density
        # configurations.
	# As a second condition, we also demand that the radius be
	# smaller than 20 Msun
        imax = outdata[:,2].argmax() 
	if(outdata[imax,3] > 20.0):
		imax = 0
		i = 0
		mmax = 0.0
		while(i < len(outdata[:,0])):
			if(outdata[i,2]>mmax and outdata[i,3]<=20.0):
				mmax = outdata[imax,2]
				imax = i
			i +=1

        print "T = %5.2f, Y_e = %5.2f" % (par1,par2)
        print "Maximum mass: M_grav = %15.6E   M_bary = %15.6E " % \
            (outdata[imax,1],outdata[imax,2])
        print "Radius at maximum mass: r_max = %15.6E cm" % (outdata[imax,3]*inv_length_gf)
        print "Maximum gravitational mass at: rho_c = %15.6E" % outdata[imax,0]
        print "                         rho_c CGS   = %15.6E" % (outdata[imax,0] * inv_rho_gf)

	sfile = open("summary_fixed_T_Ye.dat","aw") 
	outstring = "%15.6E %15.6E %15.6E %15.6E %15.6E %15.6E\n"  % \
            (par1,par2,outdata[imax,1],outdata[imax,2],outdata[imax,0]*inv_rho_gf,\
                 outdata[imax,3]*inv_length_gf)

        sfile.write(outstring)
	sfile.close()


