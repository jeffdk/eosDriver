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
        self.stopflag = False

# global constants
ggrav = 1.0
clite = 1.0
msun = 1.0
precision = 1.0e-12

def tabeos_press_eps(rho,tovinfo):
    # get pressure and eps for given density
    lrho = log10(rho)
    irho = 1 + int( (lrho - tovinfo.eoslrhomin) * \
                           tovinfo.eosdlrhoi) 

    lp = (tovinfo.eostable[irho,0] - tovinfo.eostable[irho-1,0]) \
        * tovinfo.eosdlrhoi * (lrho - tovinfo.logrhos[irho-1]) + \
        tovinfo.eostable[irho-1,0]

    leps = (tovinfo.eostable[irho,1] - tovinfo.eostable[irho-1,1]) \
        * tovinfo.eosdlrhoi * (lrho - tovinfo.logrhos[irho-1]) + \
        tovinfo.eostable[irho-1,1]
    
    press = 10.0**lp
    eps = 10.0**leps - tovinfo.eosepsshift


    return(press,eps)                  

def get_rho_eps(press,rho_old,tovinfo):
    polyK = tovinfo.polyK
    polyG = tovinfo.polyG
    tovinfo.stopflag = False
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
        lprec = precision
        while(cont):
            #print counter, rho_guess, rho_old, 10.0**tovinfo.eoslrhomin
            counter += 1
            rho_guess2 = rho_guess * 1.0001
            (xprs2,xeps2) = tabeos_press_eps(rho_guess2,tovinfo)
            (xprs,xeps) = tabeos_press_eps(rho_guess,tovinfo)
            mydpdrho = (xprs2-xprs)/(rho_guess2-rho_guess)
            if (abs(1.0-xprs/press) < lprec) :
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
#                print "count: %d rel error: %18.9E" % (counter, \
#                                                           abs(1.0-xprs/press))
            if (counter > 1000):
                fac = 0.001

            if (rho_guess <= 10.0**tovinfo.eoslrhomin):
                rho_guess = 10.0**tovinfo.eoslrhomin
                eps = 10.0**tovinfo.eostable[0,1]-tovinfo.eosepsshift
                rho = rho_guess*0.9
                tovinfo.stopflag = True
                return(rho,eps)


            if (counter > 50000):
                print "press: %18.9E" % press
                print "rho_guess: %18.9E" % (rho_guess*inv_rho_gf)
                print "rho_guess2: %18.9E" % (rho_guess2*inv_rho_gf)
                print "rel error: %18.9E" % (abs(1.0-xprs/press))
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
        if(press > tovinfo.minpress and not tovinfo.stopflag):
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
    tovinfo.stopflag = False

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
    while (isurf == 0 and i < nzones-1 and not tovinfo.stopflag):

        # in smooth part, use RK2
        if(rhos[i] > 1.0e-6): 
            tovdata[i+1,:] = tov_RK2(tovdata[i,:],rad[i],dr,tovinfo,rhos[i])
        # near the edge, use RK4 for better accuracy
        else:
            tovdata[i+1,:] = tov_RK4(tovdata[i,:],rad[i],dr,tovinfo,rhos[i])

        # check if press below 0
        if(tovdata[i+1,0] <= tovinfo.minpress or tovinfo.stopflag):
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

        if(tovinfo.stopflag):
            isurf = i+1

        # compute eps
        rhos[i+1] = tovout[i+1,0]
#        print "%d %15.6E %15.6E %15.6E %15.6E" % (i,rad[i+1],\
#                tovout[i,3],rhos[i+1]*inv_rho_gf,tovdata[i+1,0])
        
        i+=1

    if (isurf == 0):
        print "Could not solve for (entire?) TOV!"
        print "Setting isurf to nzones-1"
        isurf = nzones-1


#    print isurf, rad[isurf], rhos[isurf]
 
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
        print "working on rho_c = %18.9E" % rho_c
        (tovout,isurf,rad,dr) = tov_integrate(rho_c,tovinfo)
        outdata[i,0] = rho_c
        outdata[i,1] = tovout[isurf,3]
        outdata[i,2] = tovout[isurf,4]
        outdata[i,3] = rad[isurf]
        print "%4d %15.6E %15.6E %15.6E %15.6E" % \
            (i,outdata[i,0],outdata[i,1],outdata[i,2],outdata[i,3])

    return outdata





### rho_c = 7.715722E-4
### (tovout,isurf,rad,dr) = tov_integrate(rho_c,tovinfo)
### print tovout[isurf,3], tovout[isurf,4], rad[isurf]




