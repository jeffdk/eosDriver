from copy import deepcopy
import os
import matplotlib.pyplot as plt
import operator
import numpy
import plot_defaults
from consts import CGS_MSUN
import eosDriver
from utils import lookupIndexBisect, linInterp

fileLocation = "/home/jeff/work/tovDataThermalSupport"

files = os.listdir(fileLocation)

print files

#script = 'c30p5'
edToPlot = 2e15

modelParamsTemplate = {'M': [], 'R': [], 'ed': [], 'rhobar': [], 'rhomax': [], 'rhoMfract': []}

styles = {'3e+14': ':',
          '1e+15': '--',
          '2e+15': '-'}

colors = {'c30p0': 'g',
          'c20p0': 'b',
          'c40p0': 'r',
          'T=00010': 'm',
          'c30p5': 'c',
          'c30p10': 'k'}

symbols = {'c30p0': 's',
           'c20p0': 'v',
           'c40p0': '^',
           'T=00010': '*',
           'c30p5': 'p',
           'c30p10': 'H'}

scripts = {'c30p0': deepcopy(modelParamsTemplate.copy()),
           'c20p0': deepcopy(modelParamsTemplate.copy()),
           'c40p0': deepcopy(modelParamsTemplate.copy()),
           'T=00010': deepcopy(modelParamsTemplate.copy()),
           'c30p5': deepcopy(modelParamsTemplate.copy()),
           'c30p10': deepcopy(modelParamsTemplate.copy())}

shen = eosDriver.eosDriver('/home/jeff/work/HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5')

tfuncs = {'c30p0':  eosDriver.getTRollFunc(30.0, 0.01, 14.055, .375),
           'c20p0': eosDriver.getTRollFunc(20.0, 0.01, 13.93, .25),
           'c40p0': eosDriver.getTRollFunc(40.0, 0.01, 14.18, .5),
           'T=00010': lambda x: 0.01,
           'c30p5': eosDriver.kentaDataTofLogRhoFit2(),
           'c30p10': eosDriver.kentaDataTofLogRhoFit1()}

ye = 0.1
massFraction = .5

def thermalPressureSupport(rho, tfunc):
    shen.setState({'rho': rho, 'temp': 0.01, 'ye': ye})
    cold = shen.query('logpress', deLog10Result=True)

    shen.setState({'rho': rho, 'temp': tfunc(rho), 'ye': ye })
    hot = shen.query('logpress', deLog10Result=True)

    return hot/cold - 1.0


def findPointOfHalfM(modelData, fract=0.5):
    M = modelData['m'][-1]
    Mfract = M * fract
    #i = lookupIndexBisect(M / 2.0, modelData['m'])
    rho = linInterp(Mfract, modelData['m'], modelData['rho'])
    r = linInterp(Mfract, modelData['m'], modelData['r'])
    p = linInterp(Mfract, modelData['m'], modelData['p'])

    return {'rho': rho, 'r': r, 'p': p}

plt.figure()
for file in files:
    parts = file.split('_')
    #print parts
    ye = float(parts[3].split('=')[1])
    thisScript = ''.join((parts[4].split('.')[:-1]))
    ed = float(parts[5])
    #print ye, thisScript, ed

    #silly way to say only plot this guy
    #if not thisScript == script:
    #if thisScript == 'T=00010' or thisScript == 'c40p0' or thisScript == 'c20p0':
    #    continue
        #print "'%s' is not '%s'" % (script, thisScript)
    #    continue

    filehandle = open(fileLocation + '/' + file, 'r')
    data = {'r': [], 'rho': [], 'p': [], 'm': []}
    for line in filehandle:
        entry = line.split()
        data['r'].append(float(entry[0]) / 1.e5)
        data['rho'].append(float(entry[1]))
        data['p'].append(float(entry[2]))
        data['m'].append(float(entry[3]) / CGS_MSUN)
    M = data['m'][-1]
    R = data['r'][-1]
    rhobar = 3.* (M * CGS_MSUN) / (R * 1.e5) ** 3 / (4.0 * numpy.pi)
    rhomax = data['rho'][0]
    print thisScript, rhobar
    scripts[thisScript]['M'].append(M)
    scripts[thisScript]['R'].append(R)
    scripts[thisScript]['rhobar'].append(rhobar)
    scripts[thisScript]['ed'].append(ed)
    scripts[thisScript]['rhomax'].append(rhomax)
    scripts[thisScript]['rhoMfract'].append(findPointOfHalfM(data, massFraction)['rho'])
    label = None
    if styles[str(ed)] == '-':
        label = thisScript
    # plt.subplot(2, 2, 1 )
    # xaxis = 'r'
    # yaxis = 'rho'
    # #plt.xlabel('Radius (km)')
    # plt.ylabel(r'Density $\rho_b$ (cgs)')
    # plt.semilogy(numpy.array(data[xaxis]), data[yaxis],
    #              label=label, ls=styles[str(ed)], color=colors[thisScript])
    #
    # plt.subplot(2, 2, 2)
    # xaxis = 'r'
    # yaxis = 'rho'
    # plt.xlabel('Radius (km)')
    # plt.ylabel(r'Density $\rho_b$ (cgs)')
    # plt.plot(numpy.array(data[xaxis]), data[yaxis],
    #              label=label, ls=styles[str(ed)], color=colors[thisScript])
    #
    # # plt.subplot(2,2,3)
    # # xaxis = 'r'
    # # yaxis = 'DURRR'
    # # plt.xlabel(r'Radius (km)')
    # # plt.ylabel('Thermal Pressure Support')
    # # tfunc = numpy.frompyfunc(lambda rho: thermalPressureSupport(rho, tfuncs[thisScript]), 1, 1)
    # # plt.semilogy(numpy.array(data[xaxis]), tfunc(numpy.array(data['rho'])),
    # #              label=label, ls=styles[str(ed)], color=colors[thisScript])
    #
    # plt.subplot(2, 2, 3)
    # xaxis = 'r'
    # yaxis = 'rho'
    # plt.xlabel('Radius (km)')
    # plt.ylabel(r'$\rho_br^2$')
    # plt.plot(numpy.array(data[xaxis]), numpy.array(data[yaxis]) * numpy.array(data[xaxis]) * numpy.array(data[xaxis]),
    #              label=label, ls=styles[str(ed)], color=colors[thisScript])
    #

    #plt.subplot(2, 2, 4)
    xaxis = 'r'
    yaxis = 'rho'
    plt.xlabel('Radius (km)')
    plt.ylabel(r'$\rho_br^2$')
    plt.plot(numpy.array(data[xaxis]), numpy.array(data[yaxis]) * numpy.array(data[xaxis]) * numpy.array(data[xaxis]),
                 label=label, ls=styles[str(ed)], color=colors[thisScript])


plt.tight_layout(pad=1.05, h_pad=0.1, w_pad=0.1)
plt.subplots_adjust(left=0.06, bottom=0.08, right=0.98, top=0.98, wspace=0.15, hspace=0.17)
fig = plt.gcf()
fig.set_size_inches(24, 13.5)
plt.legend()
plt.show()

plt.clf()



#numpyfy data
for script, value in scripts.items():
    for key in value.keys():
        value[key] = numpy.array(value[key])

plt.subplot(2, 2, 1)
xaxis = 'R'
yaxis = 'rhoMfract'
plt.xlabel('Radius (km)')
plt.ylabel(r'$\rho_b$ at which $m/M = %s$' % massFraction)
for key, value in scripts.items():
    table = [col for col in value.values()]
    table = zip(*table)
    table = sorted(table, key=operator.itemgetter(0))
    table = zip(*table)
    table = numpy.array(table)
    label = key
    plt.scatter(value[xaxis], value[yaxis],
                c=colors[key], label=label, marker=symbols[key], s=100)
    plt.plot(table[2], table[5], color=colors[key])
plt.legend(loc=4)
#plt.show()

plt.subplot(2, 2, 2)
xaxis = 'R'
yaxis = 'M'
plt.xlabel('Radius (km)')
plt.ylabel(r'Mass ($M_\odot$)')
for key, value in scripts.items():
    table = [col for col in value.values()]
    table = zip(*table)
    table = sorted(table, key=operator.itemgetter(0))
    table = zip(*table)
    table = numpy.array(table)
    label = key
    plt.scatter(value[xaxis], value[yaxis],
                c=colors[key], label=label, marker=symbols[key], s=100)
    plt.plot(table[2], table[1], color=colors[key])
plt.legend()
#plt.show()

plt.subplot(2, 2, 3)
xaxis = 'rhobar'
yaxis = 'M'
plt.xlabel(r'Average density $\bar{\rho_b}$ (cgs)')
plt.ylabel(r'Mass ($M_\odot$)')
for key, value in scripts.items():
    table = [col for col in value.values()]
    table = zip(*table)
    table = sorted(table, key=operator.itemgetter(0))
    table = zip(*table)
    table = numpy.array(table)
    label = key
    plt.scatter(value[xaxis], value[yaxis],
                c=colors[key], label=label, marker=symbols[key], s=100)
    plt.plot(table[3], table[1], color=colors[key])
plt.legend(loc=2)
#plt.show()

plt.subplot(2, 2, 4)
xaxis = 'rhobar'
yaxis = 'M'
plt.xlabel(r'$\bar{\rho_b}$/$\rho_{\mathrm{max}}$')
plt.ylabel(r'Mass ($M_\odot$)')
for key, value in scripts.items():
    print value
    table = [col for col in value.values()]
    table = zip(*table)
    table = sorted(table, key=operator.itemgetter(0))
    table = zip(*table)
    table = numpy.array(table)
    label = key
    print table
    plt.scatter(numpy.array(value[xaxis]) / numpy.array(value['rhomax']), value[yaxis],
                c=colors[key], label=label, marker=symbols[key], s=100)
    plt.plot(table[3] / table[4], table[1], color=colors[key])
plt.legend(loc=2)


plt.tight_layout(pad=1.05, h_pad=0.1, w_pad=0.1)
plt.subplots_adjust(left=0.06, bottom=0.08, right=0.98, top=0.98, wspace=0.15, hspace=0.2)
fig = plt.gcf()
fig.set_size_inches(24, 13.5)
plt.show()
