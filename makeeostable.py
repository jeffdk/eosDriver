import sys
from units import *
from numpy import linspace, zeros, log10
from eosDriver import eosDriver, getTRollFunc, kentaDataTofLogRhoFit1, kentaDataTofLogRhoFit2
import numpy

def makeeostable(nrhos,rhomin,rhomax,myeos,myeosname,mytype,par1,par2):
    assert isinstance(myeos, eosDriver)

    logrhos = linspace(log10(rhomin),log10(rhomax),nrhos)
    eostable = zeros((nrhos,2))
    eostable2 = zeros((nrhos,4))
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

    elif(mytype == 'c30p0_fixed_Ye'):

        ye = par2

        tempfunc = getTRollFunc(30.0,0.01,14.055,0.375)

        eostablename = myeosname+"_eostable_Ye=%06.3f_c30p0.dat" % \
            (par2)
        energy_shift = myeos.h5file['energy_shift'][0]
        
        for i in range(nrhos):
            temp = tempfunc(logrhos[i])
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

    elif(mytype == 'c20p0_fixed_Ye'):

        ye = par2

        tempfunc = getTRollFunc(20.0,0.01,14.0-0.07,0.125)

        eostablename = myeosname+"_eostable_Ye=%06.3f_c20p0.dat" % \
            (par2)
        energy_shift = myeos.h5file['energy_shift'][0]
        
        for i in range(nrhos):
            temp = tempfunc(logrhos[i])
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


    elif(mytype == 'c40p0_fixed_Ye'):

        ye = par2

        tempfunc = getTRollFunc(40.0,0.01,14.25-0.07,0.5)

        eostablename = myeosname+"_eostable_Ye=%06.3f_c40p0.dat" % \
            (par2)
        energy_shift = myeos.h5file['energy_shift'][0]
        
        for i in range(nrhos):
            temp = tempfunc(logrhos[i])
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

    elif(mytype == 'c30p10_fixed_Ye'):

        ye = par2

        tempfunc = kentaDataTofLogRhoFit1()

        eostablename = myeosname+"_eostable_Ye=%06.3f_c30p10.dat" % \
            (par2)
        energy_shift = myeos.h5file['energy_shift'][0]
        
        for i in range(nrhos):
            temp = tempfunc(logrhos[i])
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

    elif(mytype == 'c30p5_fixed_Ye'):

        ye = par2

        tempfunc = kentaDataTofLogRhoFit2()

        eostablename = myeosname+"_eostable_Ye=%06.3f_c30p5.dat" % \
            (par2)
        energy_shift = myeos.h5file['energy_shift'][0]
        
        for i in range(nrhos):
            temp = tempfunc(logrhos[i])
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
        eostablename = myeosname+"_eostable_BetaEq_T=%06.3f.dat" % (par1)

        energy_shift = myeos.h5file['energy_shift'][0]
        temp = par1
        lastPress = -1e-300
        lastEps = None
        for i in range(nrhos):
            ye = myeos.setBetaEqState({'rho': 10.0**logrhos[i],
                                       'temp': temp})
            
            (press,eps) = myeos.query(['logpress','logenergy'])
            if press < lastPress:
                print "WOW, PRESSURE DECREASE WHEN INCREASING DENSITY, THAT'S FUCKED"
                print "MANUALLY SETTING VALUES TO LAST VALUES"
                press = lastPress
                eps = lastEps
            lastPress = press
            lastEps = eps
            # convert units
            eostable[i,0] = log10(10.0**press * press_gf)
            eostable[i,1] = log10(10.0**eps * eps_gf)
            eostable2[i,0] = eostable[i,0]
            eostable2[i,1] = eostable[i,1]
            eostable2[i,2] = ye
            eostable2[i,3] = temp


            print "Making EOS: %15.6E %15.6E %15.6E" % (10.0**logrhos[i],temp,ye)
        energy_shift = energy_shift*eps_gf

        # write EOS table
        eosfile=open(eostablename,"w")
        esstring = "%18.9E\n" % (energy_shift/eps_gf)
        eosfile.write(esstring)
        for i in range(len(eostable2[:,0])):
            sline = "%22.14E %22.14E %22.14E %22.14E %22.14E\n" % \
                (logrhos[i],eostable2[i,0],eostable2[i,1],eostable2[i,2],eostable2[i,3])
            eosfile.write(sline)
        eosfile.close()



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

