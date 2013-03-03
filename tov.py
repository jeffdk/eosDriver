#!/opt/local/bin/python

import sys
from units import *
from numpy import linspace, zeros, log10, pi, sqrt
from eosDriver import eosDriver
import makeeostable



class tovinfoclass(object):

    def __init__(self):
        self.eostype = 1
        self.polyG = 0.0
        self.polyK = 0.0
        self.nzones = 0
        self.rmax = 0.0
        self.nrhos = 600
        self.eostable = zeros((self.nrhos,2))
        self.eoslrhomin = 0.0
        self.eoslrhomax = 0.0
        self.eosepsshift = 0.0
        self.eosdlrho = 0.0
        self.eosdlrhoi = 0.0
        self.minpress = 0.0

# global constants
ggrav = 1.0
clite = 1.0
msun = 1.0
precision = 1.0e-10

def tabeos_press_eps(rho,tovinfo):
    # get pressure and eps for given density
    lrho = log10(rho)
    irho = 1 + int( (lrho - tovinfo.eoslrhomin) * \
                           tovinfo.eosdlrhoi) 

#    print tovinfo.logrhos[irho],lrho,tovinfo.logrhos[irho-1]

    lp = (tovinfo.eostable[irho,0] - tovinfo.eostable[irho-1,0]) \
        * tovinfo.eosdlrhoi * (lrho - tovinfo.logrhos[irho-1]) + \
        tovinfo.eostable[irho-1,0]
    
    press = 10.0**lp

    leps = (tovinfo.eostable[irho,1] - tovinfo.eostable[irho-1,1]) \
        * tovinfo.eosdlrhoi * (lrho - tovinfo.logrhos[irho-1]) + \
        tovinfo.eostable[irho-1,1]
    
    eps = 10.0**leps - tovinfo.eosepsshift

    return(press,eps)                  

def get_rho_eps(press,rho_old,tovinfo):
    polyK = tovinfo.polyK
    polyG = tovinfo.polyG
    if tovinfo.eostype == 1:
        # analytic polytropic case
        rho = (press / polyK)**(1.0/polyG)
        eps = press / (polyG - 1.0) / rho
    else:
        # tabulated EOS
        rho_guess = rho_old
        cont = True
        counter = 0
        fac = 1.0
        while(cont):
            #print counter, rho_guess, rho_old, 10.0**tovinfo.eoslrhomin
            counter += 1
            rho_guess2 = rho_guess * 1.0001
            (xprs2,xeps2) = tabeos_press_eps(rho_guess2,tovinfo)
            (xprs,xeps) = tabeos_press_eps(rho_guess,tovinfo)
            mydpdrho = (xprs2-xprs)/(rho_guess2-rho_guess)
            if (abs(1.0-xprs/press) < precision) :
                cont = False
                rho = rho_guess
                eps = xeps
            else:
                if(fac*(press-xprs)/mydpdrho/rho_guess > 0.1):
                    rho_guess = 0.99*rho_guess
                else:
                    rho_guess = rho_guess + fac*(press-xprs)/mydpdrho

            if (counter > 100):
                fac = 0.01

            if (rho_guess <= 10.0**tovinfo.eoslrhomin):
                rho_guess = 10.0**tovinfo.eoslrhomin
                press = tovinfo.minpress
                eps = 10.0**tovinfo.eostable[0,1]-tovinfo.eosepsshift
                rho = rho_guess*0.99
                return(rho,eps)

            if (counter > 10000):
                print "error in rho(press) iteration"
                sys.exit()

    return (rho,eps)



def set_grid(rmax,nzones):
    rad = linspace(0.0,rmax,nzones)
    dr = rad[1]-rad[0]
    return (rad,dr)

def tov_RHS(r,data,tovinfo,rho_old):

    rhs = zeros(9)

    mass = data[1]
    press = max(data[0],tovinfo.minpress)

    (rho,eps) = get_rho_eps(press,rho_old,tovinfo)
#    print r, press, rho_old, rho, eps

    u = rho * (1.0 + eps/clite**2)

    if(r < 1.0e-5 and mass < 1.0e-5*msun):
        rhs[0] = 0.0
#        rhs[1] = 4.0*pi*r**2 * u
#        rhs[2] = 4.0*pi*r**2 * rho
        rhs[1] = 0.0
        rhs[2] = 0.0
        rhs[3] = 0.0
        rhs[4] = 0.0
        rhs[5] = 0.0
        rhs[6] = 0.0
    else:
        # RHS inside the star -- solver for
        # hydrostatic equilibrium
        if(press > tovinfo.minpress):
            rhs[0] = - ggrav / r**2 * \
                (rho + rho*eps/clite**2 + press/clite**2) * \
                (mass + 4.0*pi*r**3 * press/clite**2) * \
                1.0/(1.0 - 2.0*ggrav*mass/r/clite**2)
        else:
            rhs[0] = 0.0

        # gravitational mass
        rhs[1] = 4.0*pi*r**2 * u

        # baryonic mass
        rhs[2] = 4.0*pi*r**2 * rho * \
            1.0/sqrt(1.0-2.0*ggrav*mass/r/clite**2)

        # for other equations (not implemented)
        rhs[3] = 0.0
        rhs[4] = 0.0
        rhs[5] = 0.0
        rhs[6] = 0.0

    return rhs

def tov_RK2(old_data,r,dr,tovinfo,rho_old):

    ktemp = zeros((9,2))
    ktemp[:,0] = dr*tov_RHS(r,old_data,tovinfo,rho_old)
    ktemp[:,1] = dr*\
        tov_RHS(r+0.5*dr,old_data+0.5*ktemp[:,0],tovinfo,rho_old)

    return old_data + ktemp[:,1]
    
def tov_RK3(old_data,r,dr,tovinfo,rho_old):
    ktemp = zeros((9,3))

    # first step
    ktemp[:,0] = dr * tov_RHS(r,old_data,tovinfo,rho_old)

    # second_step
    ktemp[:,1] = dr * tov_RHS(r + 0.5*dr, \
                                  old_data + 0.5*ktemp[:,0],tovinfo,rho_old)

    # third step
    ktemp[:,2] = dr *\
        tov_RHS(r + dr, old_data - ktemp[:,0] + 2.0*ktemp[:,1],\
                    tovinfo,rho_old)

    return old_data + \
        1.0/6.0 * (ktemp[:,0] + 4.0*ktemp[:,1] + ktemp[:,2])

def tov_RK4(old_data,r,dr,tovinfo,rho_old):
    ktemp = zeros((9,4))

    # first step
    ktemp[:,0] = dr * tov_RHS(r,old_data,tovinfo,rho_old)

    # second_step
    ktemp[:,1] = dr * \
        tov_RHS(r + 0.5*dr, old_data + 0.5*ktemp[:,0],tovinfo,rho_old)

    # third_step
    ktemp[:,2] = dr * \
        tov_RHS(r + 0.5*dr, old_data + 0.5*ktemp[:,1],tovinfo,rho_old)

    # fourth step
    ktemp[:,3] = dr * \
        tov_RHS(r + dr, old_data + ktemp[:,2],tovinfo,rho_old)

    return old_data + \
        1.0/6.0 * (ktemp[:,0] + 2*ktemp[:,1] + 2*ktemp[:,2] +\
                   ktemp[:,3])



def tov_integrate(rho_c,tovinfo):

    # set up grid
    (rad,dr) = set_grid(tovinfo.rmax,tovinfo.nzones)

    nzones = tovinfo.nzones
    # initialize some variables
    tovdata = zeros((nzones,9))
    rhos = zeros(nzones)
    # 0 -- press
    # 1 -- mgrav
    # 2 -- mbary 


    tovout = zeros((nzones,10))
    # 0 -- rho
    # 1 -- press
    # 2 -- eps
    # 3 -- mgrav
    # 4 -- mbary

    polyK = tovinfo.polyK
    polyG = tovinfo.polyG
     
    # central values
    if tovinfo.eostype == 1:
        tovdata[0,0] = polyK * rho_c**polyG
    elif tovinfo.eostype == 3:
        (tovdata[0,0],foo) = tabeos_press_eps(rho_c,tovinfo)
    else:
        print "eostype %d not implemented" % tovinfo.eostype
        sys.exit()

    isurf = 0
    rhos[0] = rho_c
    tovout[0,0] = rho_c
    tovout[0,1] = tovdata[0,0]
    tovout[0,2] = rhos[0]
    tovout[0,6] = 0.0
    i = 0
    while (isurf == 0 and i < nzones):

        tovdata[i+1,:] = tov_RK2(tovdata[i,:],rad[i],dr,tovinfo,rhos[i])

        # check if press below 0
        if(tovdata[i+1,0] <= tovinfo.minpress*10.0):
            tovdata[i+1,0] = tovinfo.minpress
            if(isurf == 0):
                isurf = i+1

        if (i+1 >= isurf and isurf > 0):
            tovout[i+1,3] = tovdata[i+1,1]
            tovout[i+1,4] = tovdata[i+1,2]
        else:
            tovout[i+1,3] = tovdata[i+1,1]
            tovout[i+1,4] = tovdata[i+1,2]

        # press and mass
        tovout[i+1,1] = tovdata[i+1,0]

        # compute density
        (tovout[i+1,0],tovout[i+1,2]) = \
               get_rho_eps(tovdata[i+1,0],rhos[i],tovinfo)
        # compute eps
        rhos[i+1] = tovout[i+1,0]
#        print i,rad[i+1],tovout[i,3],rhos[i+1],tovdata[i+1,0],tovinfo.minpress
        
        i+=1


#    print rad[isurf], rhos[isurf]
    
    return (tovout,isurf,rad,dr)



def tov_sequence(rho1,rho2,nsteps,tovinfo):
    nrho = nsteps
    outdata = zeros((nrho,4))
    # outdata format
    # 0 -- density
    # 1 -- grav mass
    # 2 -- bary mass
    # 3 -- areal radius

    # set up rhos
    rhos = zeros(nrho)
    dlrho = (log10(rho2)-log10(rho1)) / (nrho-1)
    for i in range(len(rhos)):
        rhos[i] = log10(rho1) + i*dlrho

    # compute TOV solutions        
    for i in range(nsteps):
        rho_c = 10.0**rhos[i]*rho_gf
        (tovout,isurf,rad,dr) = tov_integrate(rho_c,tovinfo)
        outdata[i,0] = rho_c
        outdata[i,1] = tovout[isurf,3]
        outdata[i,2] = tovout[isurf,4]
        outdata[i,3] = rad[isurf]
        print "%4d %15.6E %15.6E %15.6E %15.6E" % \
            (i,outdata[i,0],outdata[i,1],outdata[i,2],outdata[i,3])

    return outdata

sfile = open("summary.dat","aw") 

myeos = eosDriver('LS220_234r_136t_50y_analmu_20091212_SVNr26.h5')


tovinfo = tovinfoclass()
tovinfo.polyK = 100.0
tovinfo.polyG = 2.0
tovinfo.nzones = 10000
tovinfo.rmax = 50.0
tovinfo.eostype = 3

temps = [0.5,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0]
yes = [0.1,0.15,0.2,0.25,0.3]

for ii in range(len(temps)):
    for jj in range(len(yes)):
        mytype = "fixed_ye_temp"
        par1 = 0.5
        par2 = 0.15
        print "T = %5.2f, Y_e = %5.2f" % (par1,par2)
        print "Preparing EOS table: ",mytype
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

        outdata = tov_sequence(2.0e14,7.0e15,50,tovinfo)
    
        filename = mytype+"_T%06.3f_Ye%06.3f.dat" % (par1, par2)
        outfile=open(filename,"w")
        for i in range(len(outdata[:,0])):
            line = "%15.6E %15.6E\n" % (outdata[i,3]*inv_length_gf,outdata[i,1])
            outfile.write(line)
        outfile.close()

        imax = outdata[:,2].argmax() 
        print "T = %5.2f, Y_e = %5.2f" % (par1,par2)
        print "Maximum mass: M_grav = %15.6E   M_bary = %15.6E " % \
            (outdata[imax,1],outdata[imax,2])
        print "Radius at maximum mass: r_max = %15.6E cm" % (outdata[imax,3]*inv_length_gf)
        print "Maximum gravitational mass at: rho_c = %15.6E" % outdata[imax,0]
        print "                         rho_c CGS   = %15.6E" % (outdata[imax,0] * inv_rho_gf)

        outstring = "%15.6E %15.6E %15.6E %15.6E %15.6E %15.6E\n"  % \
            (par1,par2,outdata[imax,1],outdata[imax,2],outdata[imax,0]*inv_rho_gf,\
                 outdata[imax,3]*inv_length_gf)

        sfile.write(outstring)


sfile.close()




### rho_c = 7.715722E-4
### (tovout,isurf,rad,dr) = tov_integrate(rho_c,tovinfo)
### print tovout[isurf,3], tovout[isurf,4], rad[isurf]




