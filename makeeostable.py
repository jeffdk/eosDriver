import sys
from units import *
from numpy import linspace, zeros, log10
from eosDriver import eosDriver
import numpy

def makeeostable(nrhos,rhomin,rhomax,myeos,myeosname,mytype,par1,par2):
    assert isinstance(myeos, eosDriver)

    logrhos = linspace(log10(rhomin),log10(rhomax),nrhos)
    eostable = zeros((nrhos,2))
    energy_shift = 0.0

    print mytype
    
    if(mytype == 'fixed_ye_temp'):

        ye = par2
        temp = par1
        
        eostablename = myeosname+"_eostable_Ye=%06.3f_T=%06.3f.dat" % \
            (par2,par1)
        energy_shift = myeos.h5file['energy_shift'][0]
        
        for i in range(nrhos):
            myeos.setState({'rho': 10.0**logrhos[i],
                            'ye': ye,
                            'temp': temp})
            
            press = myeos.query('logpress')

            myeos.setState({'rho': 10.0**logrhos[i],
                            'ye': ye,
                            'temp': temp})

            eps = myeos.query('logenergy')
            
            # convert units
            eostable[i,0] = log10(10.0**press * press_gf)
            eostable[i,1] = log10(10.0**eps * eps_gf)

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

    elif(mytype == 'fixed_temp_betaeq'):

        eostablename = myeosname+"_eostable_BetaEq_T=%06.3f.dat" % (par2,par1)

        energy_shift = myeos.h5file['energy_shift'][0]
        temp = par1
        for i in range(nrhos):
            ye = myeos.setBetaEqState({'rho': 10.0**logrhos[i],
                                       'temp': temp})
            
            (press,eps) = myeos.query(['logpress','logenergy'])

            # convert units
            eostable[i,0] = log10(10.0**press * press_gf)
            eostable[i,1] = log10(10.0**eps * eps_gf)

            print "Making EOS: %15.6E %15.6E %15.6E" % (10.0**logrhos[i],temp,ye)
        energy_shift = energy_shift*eps_gf

    elif(mytype == 'fixed_ye_entropy'):

        ye = par2
        entropy = par1
        
        energy_shift = myeos.h5file['energy_shift'][0]

        eostablename = "eostable_ye=%06.3f_s=%06.3f.dat" % (par2,par1)
        print eostablename
        try: 
            eosin = (numpy.loadtxt(eostablename))
            eostable[:,0] = eosin[:,1]
            eostable[:,1] = eosin[:,2]
            
        except IOError:
            for i in range(nrhos):
                temp = 0.0
                temp = myeos.getTemperatureFromQuantityTYe({'rho': 10.0**logrhos[i],
                                                         'ye': ye,
                                                         'temp': temp},'entropy', entropy)
                myeos.setState({'rho': 10.0**logrhos[i],
                                'ye' : ye,
                                'temp': temp})
                
                (press,eps) = myeos.query(['logpress','logenergy'])
                # convert units
                eostable[i,0] = log10(10.0**press * press_gf)
                eostable[i,1] = log10(10.0**eps * eps_gf)

                print "Making EOS: %15.6E %15.6E %15.6E" % (10.0**logrhos[i],temp,ye)


        energy_shift = energy_shift*eps_gf

        # write EOS table
        eosfile=open(eostablename,"w")
        for i in range(len(eostable[:,0])):
            sline = "%15.6E %15.6E %15.6E\n" % \
                (logrhos[i],eostable[i,0],eostable[i,1])
            eosfile.write(sline)
        eosfile.close()


    elif(mytype == 'fixed_entropy_betaeq'):

        energy_shift = myeos.h5file['energy_shift'][0]
        entropy = par1

        eostablename = "eostable_betaeq_s=%06.3f.dat" % (par1)
        print eostablename
        try: 
            eosin = (numpy.loadtxt(eostablename))
            eostable[:,0] = eosin[:,1]
            eostable[:,1] = eosin[:,2]

        except IOError:
            for i in range(nrhos):
                (ye,temp) = \
                    myeos.setConstQuantityAndBetaEqState({'rho': 10.0**logrhos[i]},\
                                                             'entropy',par1)
            
                (press,eps) = myeos.query(['logpress','logenergy'])
            
                # convert units
                eostable[i,0] = log10(10.0**press * press_gf)
                eostable[i,1] = log10(10.0**eps * eps_gf)

                print "Making EOS: %15.6E %15.6E %15.6E" % (10.0**logrhos[i],temp,ye)

        energy_shift = energy_shift*eps_gf

        # write EOS table
        eosfile=open(eostablename,"w")
        for i in range(len(eostable[:,0])):
            sline = "%15.6E %15.6E %15.6E\n" % \
                (logrhos[i],eostable[i,0],eostable[i,1])
            eosfile.write(sline)
        eosfile.close()


    elif(mytype == 'poly_G2_K100'):
        eostable[:,0] = log10(100.0*(10.0**logrhos[:]*rho_gf)**2.0)
        eostable[:,1] = log10(10.0**eostable[:,0]/(10.0**logrhos[:]*rho_gf) )

    else:
        print "This kind of table can't be done yet: ",mytype
        sys.exit()
        

    logrhos = log10(10.0**logrhos * rho_gf)
    dlrho = logrhos[1]-logrhos[0]
    return (eostable,energy_shift,dlrho,logrhos)

