#!/opt/local/bin/python
import sys
from units import *
from numpy import linspace, zeros, log10, pi, sqrt, exp
from eosDriver import eosDriver, getTRollFunc, kentaDataTofLogRhoFit1, kentaDataTofLogRhoFit2
import makeeostable
from tov import *

#define EOSs

EOSlist = [ ["LS220","LS220_234r_136t_50y_analmu_20091212_SVNr26.h5"]] 
#EOSlist = [ ["HShen", "HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5"]]
#EOSlist = [ ["HShen", "HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5"], ["LS220","LS220_234r_136t_50y_analmu_20091212_SVNr26.h5"]]

#EOSlist = [             ["GShenNL3",
#             "GShen_NL3EOS_rho280_temp180_ye52_version_1.1_20120817.h5"]]

#EOSlist = [             ["GShenFSU2.1",
#             "GShenFSU_2.1EOS_rho280_temp180_ye52_version_1.1_20120824.h5"] ]

#EOSlist = [            ["Hempel_DD2",
#             "Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5"]]

#EOSlist = [             ["Hempel_SFHo",
#             "Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5"]]

#EOSlist = [ ["Hempel_SFHx",
#             "Hempel_SFHxEOS_rho234_temp180_ye60_version_1.1_20120817.h5"]]

#EOSlist = [ ["LS375","LS375_234r_136t_50y_analmu_20091212_SVNr26.h5"]]

nothing = [ ["HShen", "HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5"],
            ["GShenFSU2.1",
             "GShenFSU_2.1EOS_rho280_temp180_ye52_version_1.1_20120824.h5"],
            ["GShenNL3",
             "GShen_NL3EOS_rho280_temp180_ye52_version_1.1_20120817.h5"],
            ["Hempel_DD2",
             "Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5"],
            ["Hempel_SFHo",
             "Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5"],
            ["Hempel_SFHx",
             "Hempel_SFHxEOS_rho234_temp180_ye60_version_1.1_20120817.h5"],
            ["LS375","LS375_234r_136t_50y_analmu_20091212_SVNr26.h5"],
          ]

yes = [0.1,0.15]
tmin = 0.5
tmax = 0.5
dtemp = 0.5
ntemp = int((tmax-tmin)/dtemp)+1
temps = zeros(ntemp)
for i in range(ntemp):
	temps[i] = tmin + dtemp*i


rhomin = 1.0e6
rhomax = 0.0
nrhos = 400

def fixed_ye_temp(EOSlist,temps,yes):
    for ieos in range(len(EOSlist)):
        myeos = eosDriver(EOSlist[ieos][1])
        rhomax = (10.0e0**max(myeos.h5file['logrho']))*0.995
        for ii in range(len(temps)):
            for jj in range(len(yes)):
                mytype = "fixed_ye_temp"
                par1 = temps[ii]
                par2 = yes[jj]
                print "T = %5.2f, Y_e = %5.2f" % (par1,par2)
                print "Preparing EOS table: ",EOSlist[ieos][0], mytype
                makeeostable.makeeostable(\
                    nrhos,rhomin,rhomax,myeos,EOSlist[ieos][0],\
                        mytype,par1,par2)

        del myeos
  
def fixed_temp_betaeq(EOSlist,temps):
	for ieos in range(len(EOSlist)):
		myeos = eosDriver(EOSlist[ieos][1])
		rhomax = (10.0e0**max(myeos.h5file['logrho']))*0.995
		for ii in range(len(temps)):
			mytype = "fixed_temp_betaeq"
			par1 = temps[ii]
			par2 = 0.0
			print "T = %5.2f, betaeq" % (par1)
			print "Preparing EOS table: ",EOSlist[ieos][0], mytype
			makeeostable.makeeostable(\
				nrhos,rhomin,rhomax,myeos,EOSlist[ieos][0],\
					mytype,par1,par2)
		del myeos


def special_fixed_Ye(EOSlist,temps,yes,mytype):
	for ieos in range(len(EOSlist)):
		myeos = eosDriver(EOSlist[ieos][1])
                rhomax = (10.0e0**max(myeos.h5file['logrho']))*0.995
		for jj in range(len(yes)):
			par1 = 0.0
			par2 = yes[jj]
			print "Y_e: %15.6E" % (par2)
			print "Preparing EOS table: ",EOSlist[ieos][0], mytype
			makeeostable.makeeostable(\
				nrhos,rhomin,rhomax,myeos,EOSlist[ieos][0],\
					mytype,par1,par2)
		del myeos

def get_pnu(eta,temp,rho,rhotrap):
	mev_to_erg = 1.60217733e-6
	pi = 3.14159265358979e0
	hc_mevcm = 1.97326966e-11*2.0*pi
	pnu =  mev_to_erg * 4.0*pi/3.0 * (temp**4/hc_mevcm**3) * \
	    (21.0*pi**4 / 60.0 + 0.5*eta**2 * \
		     (pi**2 + 0.5*eta**2)) * exp(-rhotrap/rho)
	return pnu

def fixed_temp_nfbetaeq_pnu(EOSlist,temps):
	for ieos in range(len(EOSlist)):
		myeos = eosDriver(EOSlist[ieos][1])
		rhomax = (10.0e0**max(myeos.h5file['logrho']))*0.995
		logrhos = linspace(log10(rhomin),log10(rhomax),nrhos)
		eostable = zeros((nrhos,2))
		energy_shift = myeos.h5file['energy_shift'][0]
		energy_shift = energy_shift*eps_gf
		rhotrap = 10.0e0**12.5
		for ii in range(len(temps)):
			eostablename = EOSlist[ieos][0]+\
			    "_eostable_NFBetaEq_pnu_T=%06.3f.dat" % \
			    (temps[ii])
			for i in range(nrhos):
				temp_cold = 0.5e0
				temp = temps[ii]
				rho = 10.0**logrhos[i]
				ylep = myeos.setBetaEqState({'rho': rho,
							     'temp': temp_cold})

				ye2 = myeos.setBetaEqState({'rho': rho,
							   'temp': temp})

				ye1 = myeos.setNuFullBetaEqState({'rho': rho,
								  'temp': temp},ylep)

				ye = ye2*(1-exp(-rhotrap/rho)) + ye1*(exp(-rhotrap/rho))
				
				
				myeos.setState({'rho': rho,
						'temp': temp, 
						'ye' : ye})
				
				print "Making EOS: %15.6E %15.6E %15.6E %15.6E %15.6E %15.6E" %\
				    (10.0**logrhos[i],temp,ye,ye1,ye2,ylep)

				(press,eps,munu) = myeos.query(['logpress','logenergy','munu'])
				eta = munu/temp 
				pnu = get_pnu(eta,temp,rho,rhotrap)
				epsnu = 3.0e0*pnu / rho
				# convert units
				eostable[i,0] = log10((10.0**press+pnu) * press_gf)
				eostable[i,1] = log10((10.0**eps+epsnu) * eps_gf)
				if (i>0 and eostable[i,0] < eostable[i-1,0]):
					eostable[i,0] = eostable[i-1,0]
					eostable[i,0] = eostable[i-1,0]
				

			# write EOS table
			eosfile=open(eostablename,"w")
			esstring = "%18.9E\n" % (energy_shift/eps_gf)
			eosfile.write(esstring)
			for i in range(len(eostable[:,0])):
				sline = "%15.6E %15.6E %15.6E\n" % \
				    (logrhos[i],eostable[i,0],eostable[i,1])
				eosfile.write(sline)
			eosfile.close()
		del myeos


def special_NFBetaEq_pnu(EOSlist,temps,mytype):
	for ieos in range(len(EOSlist)):
		myeos = eosDriver(EOSlist[ieos][1])
                rhomax = (10.0e0**max(myeos.h5file['logrho']))*0.995
		eostablename = EOSlist[ieos][0]+\
		    "_eostable_NFBetaEq_pnu_"+mytype+".dat"
		energy_shift = myeos.h5file['energy_shift'][0]
		logrhos = linspace(log10(rhomin),log10(rhomax),nrhos)
		eostable = zeros((nrhos,2))
		energy_shift = energy_shift*eps_gf
		rhotrap = 10.0e0**12.5
		if(mytype == "c30p5"):
			tempfunc = kentaDataTofLogRhoFit2()
		elif(mytype == "c30p10"):
			tempfunc = kentaDataTofLogRhoFit1()
		elif(mytype == "c30p0"):
			tempfunc = getTRollFunc(30.0,0.01,14.055,0.375)
		elif(mytype == "c20p0"):
			tempfunc = getTRollFunc(20.0,0.01,14.0-0.07,0.125)
		elif(mytype == "c40p0"):
			tempfunc = getTRollFunc(40.0,0.01,14.25-0.07,0.5)

		for i in range(nrhos):
			temp_cold = 0.5e0
			temp = tempfunc(logrhos[i])
			rho = 10.0**logrhos[i]

			ylep = myeos.setBetaEqState({'rho': rho,
						     'temp': temp_cold})

			ye2 = myeos.setBetaEqState({'rho': rho,
						    'temp': temp})

			ye1 = myeos.setNuFullBetaEqState({'rho': rho,
							  'temp': temp},ylep)

			ye = ye2*(1-exp(-rhotrap/rho)) + ye1*(exp(-rhotrap/rho))

			myeos.setState({'rho': rho,
					'temp': temp, 
					'ye' : ye})
				
			print "Making EOS: %15.6E %15.6E %15.6E %15.6E %15.6E %15.6E" %\
			    (10.0**logrhos[i],temp,ye,ye1,ye2,ylep)

			(press,eps,munu) = myeos.query(['logpress','logenergy','munu'])
			eta = munu/temp 
			pnu = get_pnu(eta,temp,rho,rhotrap)
			epsnu = 3.0e0*pnu / rho
			# convert units
			eostable[i,0] = log10((10.0**press+pnu) * press_gf)
			eostable[i,1] = log10((10.0**eps+epsnu) * eps_gf)
			if (i>0 and eostable[i,0] < eostable[i-1,0]):
				eostable[i,0] = eostable[i-1,0]
				eostable[i,0] = eostable[i-1,0]

		# write EOS table
		eosfile=open(eostablename,"w")
		esstring = "%18.9E\n" % (energy_shift/eps_gf)
		eosfile.write(esstring)
		for i in range(len(eostable[:,0])):
			sline = "%15.6E %15.6E %15.6E\n" % \
			    (logrhos[i],eostable[i,0],eostable[i,1])
			eosfile.write(sline)
		eosfile.close()
			

		del myeos




def special_BetaEq(EOSlist,temps,yes,mytype):
	for ieos in range(len(EOSlist)):
		myeos = eosDriver(EOSlist[ieos][1])
                rhomax = (10.0e0**max(myeos.h5file['logrho']))*0.995
		eostablename = EOSlist[ieos][0]+\
		    "_eostable_BetaEq_"+mytype+".dat"
		energy_shift = myeos.h5file['energy_shift'][0]
		if(mytype == "c30p5"):
			tempfunc = kentaDataTofLogRhoFit2()
		elif(mytype == "c30p10"):
			tempfunc = kentaDataTofLogRhoFit1()
		elif(mytype == "c30p0"):
			tempfunc = getTRollFunc(30.0,0.01,14.055,0.375)
		elif(mytype == "c20p0"):
			tempfunc = getTRollFunc(20.0,0.01,14.0-0.07,0.125)
		elif(mytype == "c40p0"):
			tempfunc = getTRollFunc(40.0,0.01,14.25-0.07,0.5)

		logrhos = linspace(log10(rhomin),log10(rhomax),nrhos)
		eostable = zeros((nrhos,2))
		for i in range(nrhos):
			temp = tempfunc(logrhos[i])
			ye = myeos.setBetaEqState({'rho': 10.0**logrhos[i],
						   'temp': temp})
			(press,eps) = myeos.query(['logpress','logenergy'])
			# convert units
			eostable[i,0] = log10(10.0**press * press_gf)
			eostable[i,1] = log10(10.0**eps * eps_gf)
			if (i>0 and eostable[i,0] < eostable[i-1,0]):
				eostable[i,0] = eostable[i-1,0]

			
			print "Making EOS: %15.6E %15.6E %15.6E" % (10.0**logrhos[i],temp,ye)
		energy_shift = energy_shift*eps_gf

		# write EOS table
		eosfile=open(eostablename,"w")
		esstring = "%18.9E\n" % (energy_shift/eps_gf)
		eosfile.write(esstring)
		for i in range(len(eostable[:,0])):
			sline = "%15.6E %15.6E %15.6E\n" % \
			    (logrhos[i],eostable[i,0],eostable[i,1])
			eosfile.write(sline)
		eosfile.close()
			

		del myeos


def poly(gamma,K):
	def press(rho,gamma,K):
		return K*rho**gamma
	def eps(rho,gamma,K):
		return K*rho**gamma / (gamma - 1.0) / rho
        rhomax = 8.0e15
	rhomin = 1.0e6
	logrhos = linspace(log10(rhomin*rho_gf),log10(rhomax*rho_gf),nrhos) 
	eostable = zeros((nrhos,2))
	eostable[:,0] = log10(press(10.0**logrhos[:],gamma,K))
	eostable[:,1] = log10(eps(10.0**logrhos[:],gamma,K))
	# write EOS table
	sgamma = "%05.4f" % (gamma)
	sK = "%15.8E" % (K)
	eostablename = "poly"+"_eostable_G"+sgamma+"_K"+sK+".dat"
	eosfile=open(eostablename,"w")
	esstring = "%18.7E\n" % (0.0)
	eosfile.write(esstring)
	for i in range(len(eostable[:,0])):
		sline = "%15.6E %15.6E %15.6E\n" % \
		    (log10(10.0**logrhos[i]/rho_gf),eostable[i,0],eostable[i,1])
		eosfile.write(sline)
	eosfile.close()
	

#poly(2.75,30000.0)

#fixed_ye_temp(EOSlist,temps,yes)

fixed_temp_betaeq(EOSlist,temps)

#fixed_temp_nfbetaeq_pnu(EOSlist,temps)

#mytypes = ["c30p5_fixed_Ye","c30p10_fixed_Ye","c30p0_fixed_Ye",
#	   "c20p0_fixed_Ye", "c40p0_fixed_Ye"]

#for t in mytypes:
#    special_fixed_Ye(EOSlist,temps,yes,t)

#mytypes2 = ["c20p0","c30p5","c30p10","c30p0","c40p0"]
#for t in mytypes2:
#    special_NFBetaEq_pnu(EOSlist,temps,t)
