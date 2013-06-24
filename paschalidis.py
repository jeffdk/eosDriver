
# This script computes the EOS that is equivalent to what is
# used by Paschalidis et al. in PRD 86, 064032 (2012), but
# with the contribution of neutrinos corrected, which have
# spin statistical weight = 1 (electrons/positrons have 2)


#!/opt/local/bin/python
import sys
from units import *
from consts import CGS_H, CGS_EV, CGS_C
from numpy import linspace, zeros, log10, pi, sqrt, exp
import pylab as pl

rhomin = 1.0e6
rhomax = 5.0e15
nrhos = 400
gamma = 2
kappa = 393.9e0
massn_cgs = 1.674927211e-24
rho_trap = 10.0**(12.5)
#h_times_clite = CGS_H / (CGS_EV * 1.e6) * CGS_C
#mev_to_erg = CGS_EV * 1.e6
def press(r,Tmev):
    ppoly = kappa * (r*rho_gf)**gamma * inv_press_gf
    ppairbase = pi**5 / h_times_clite**3 * mev_to_erg * 7.0/45.0 * Tmev**4
    ppair = ppairbase * (1 + 7.0/8.0*(2 + 3*exp(-rho_trap/r)))
    ppgas = r/massn_cgs * Tmev * mev_to_erg
    return ppoly,ppair,ppgas

# print press(1.0e14,40.0)

# logrhos = linspace(log10(rhomin),log10(rhomax),nrhos)


# logpress = zeros(len(logrhos))
# logpress1 = zeros(len(logrhos))
# logpress2 = zeros(len(logrhos))
# logpress3 = zeros(len(logrhos))

# temp = 40.0
# for i in range(len(logrhos)):
#     rho = 10.0**logrhos[i]
#     res = press(rho,temp)
#     pt = res[0]+res[1]+res[2]
#     logpress[i] = log10(pt)
#     logpress1[i] = log10(res[0])
#     logpress2[i] = log10(res[1])
#     logpress3[i] = log10(res[2])

# pl.plot(logrhos,logpress)
# pl.plot(logrhos,logpress1)
# pl.plot(logrhos,logpress2)
# pl.plot(logrhos,logpress3)

# ax = pl.gca()
# ax.set_xlim([10,16])
# ax.set_ylim([30,38])

# pl.show()

