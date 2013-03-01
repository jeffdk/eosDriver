#!/opt/local/bin/python

import sys
from units import *
from scipy import *
from pylab import *


# global constants
ggrav = 1.0
clite = 1.0
msun = 1.0
min_press = 1.0e-30


def set_grid(rmax,nzones):
    rad = linspace(0.0,rmax,nzones)
    dr = rad[1]-rad[0]
    return (rad,dr)

def tov_RHS(r,data,polyG,polyK):

    rhs = zeros(9)

    mass = data[1]
    press = max(data[0],min_press)

    rho = (press / polyK)**(1.0/polyG)
    eps = press / (polyG - 1.0) / rho

    u = rho * (1.0 + eps/clite**2)

    if(r < 1.0e-10 and mass < 1.0e-10*msun):
        rhs[0] = 0.0
        rhs[1] = 4.0*pi*r**2 * u 
        rhs[2] = 4.0*pi*r**2 * rho 
        rhs[3] = 0.0
        rhs[4] = 0.0
        rhs[5] = 0.0
        rhs[6] = 0.0
    else:
        # RHS inside the star -- solver for
        # hydrostatic equilibrium
        if(press > min_press):
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

def tov_RK2(old_data,r,dr,polyG,polyK):

    ktemp = zeros((9,2))
    ktemp[:,0] = dr*tov_RHS(r,old_data,polyG,polyK)
    ktemp[:,1] = dr*\
        tov_RHS(r+0.5*dr,old_data+0.5*ktemp[:,0],polyG,polyK)

    return old_data + ktemp[:,1]
    
def tov_RK3(old_data,r,dr,polyG,polyK):
    ktemp = zeros((9,3))

    # first step
    ktemp[:,0] = dr * tov_RHS(r,old_data,polyG,polyK)

    # second_step
    ktemp[:,1] = dr * tov_RHS(r + 0.5*dr, \
                                  old_data + 0.5*ktemp[:,0],polyG,polyK)

    # third step
    ktemp[:,2] = dr *\
        tov_RHS(r + dr, old_data - ktemp[:,0] + 2.0*ktemp[:,1],\
                    polyG,polyK)

    return old_data + \
        1.0/6.0 * (ktemp[:,0] + 4.0*ktemp[:,1] + ktemp[:,2])

def tov_RK4(old_data,r,dr):
    ktemp = zeros((9,4))

    # first step
    ktemp[:,0] = dr * tov_RHS(r,old_data,polyG,polyK)

    # second_step
    ktemp[:,1] = dr * \
        tov_RHS(r + 0.5*dr, old_data + 0.5*ktemp[:,0],polyG,polyK)

    # third_step
    ktemp[:,2] = dr * \
        tov_RHS(r + 0.5*dr, old_data + 0.5*ktemp[:,1],polyG,polyK)

    # fourth step
    ktemp[:,3] = dr * \
        tov_RHS(r + dr, old_data + ktemp[:,2],polyG,polyK)

    return old_data + \
        1.0/6.0 * (ktemp[:,0] + 2*ktemp[:,1] + 2*ktemp[:,2] +\
                   ktemp[:,3])



def tov_integrate(rho_c,rmax,nzones,polyG,polyK):

    # set up grid
    (rad,dr) = set_grid(rmax,nzones)

    # initialize some variables
    tovdata = zeros((nzones,9))
    # 0 -- press
    # 1 -- mgrav
    # 2 -- mbary 

    tovout = zeros((nzones,10))
    # 0 -- rho
    # 1 -- press
    # 2 -- eps
    # 3 -- mgrav
    # 4 -- mbary

     
    # central values
    tovdata[0,0] = polyK * rho_c**polyG

    isurf = 0
    tovout[0,0] = rho_c
    tovout[0,1] = polyK * rho_c**polyG
    tovout[0,2] = tovout[0,1] / (rho_c * (polyG-1.0))
    tovout[0,6] = 0.0
    i = 0
    while (isurf == 0 and i < nzones):

        tovdata[i+1,:] = tov_RK3(tovdata[i,:],rad[i],dr,polyG,polyK)

#        print "%15.6E %15.6E %15.6E" % (rad[i],tovdata[i+1,0],tovdata[i+1,1])


        # check if press below 0
        if(tovdata[i+1,0] < 0.0):
            tovdata[i+1,0] = min_press
            if(isurf == 0):
                isurf = i+1

        if (i+1 >= isurf and isurf > 0):
            tovout[i+1,3] = tovdata[i+1,1]
            tovout[i+1,4] = tovdata[i+1,2]
            tovout[i+1,5] = tovdata[i+1,3]
        else:
            tovout[i+1,3] = tovdata[i+1,1]
            tovout[i+1,4] = tovdata[i+1,2]
            tovout[i+1,5] = tovdata[i+1,3]

        # press and mass
        tovout[i+1,1] = tovdata[i+1,0]

        # compute density
        tovout[i+1,0] = (tovdata[i+1,0]/polyK)**(1.0/polyG)
        # compute eps
        tovout[i+1,2] = tovout[i+1,1] / (polyG - 1) / tovout[i+1,0]

        i+=1
    
    return (tovout,isurf,rad,dr)



def tov_sequence(rho1,rho2,nsteps,rmax,nzones,polyG,polyK):
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
        (tovout,isurf,rad,dr) = tov_integrate(rho_c,rmax,nzones,polyG,polyK)
        outdata[i,0] = rho_c
        outdata[i,1] = tovout[isurf,3]
        outdata[i,2] = tovout[isurf,4]
        outdata[i,3] = rad[isurf]
        print "%4d %15.6E %15.6E %15.6E %15.6E" % \
            (i,outdata[i,0],outdata[i,1],outdata[i,2],outdata[i,3])

    return outdata


# grid
rmax = 100.0
nzones = 10000

# EOS
polyG = 2.00
polyK = 123.0

rho_c = 8.0e14*rho_gf

#(tovout,isurf,rad,dr) = tov_integrate(rho_c,rmax,nzones,polyG,polyK)
#print tovout[isurf,3], tovout[isurf,4]

outdata = tov_sequence(1.0e14,1.0e16,100,rmax,nzones,polyG,polyK)
    
for i in range(len(outdata[:,0])):
    print "%4d %15.6E %15.6E %15.6E %15.6E" % \
        (i,outdata[i,0],outdata[i,1],outdata[i,2],outdata[i,3])


imax = outdata[:,2].argmax() 
print "Maximum mass: M_grav = %15.6E   M_bary = %15.6E " % \
    (outdata[imax,1],outdata[imax,2])
print "Radius at maximum mass: r_max = %15.6E" % outdata[imax,3]
print "Maximum gravitational mass at: rho_c = %15.6E" % outdata[imax,0]




#plot(rad,exp(tov_star[:,5])/sqrt(tov_star[:,9]),label="lapse")
#legend()
#show()

#plot(rad,tov_star[:,8],label="scalar field")
#legend()
#show()
