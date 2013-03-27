import os
import matplotlib.pyplot as plt

fileLocation = "/home/jeff/work/tovDataThermalSupport"

files = os.listdir(fileLocation)

print files

#script = 'c30p5'
edToPlot = 2e15


for file in files:
    parts = file.split('_')
    #print parts
    ye = float(parts[3].split('=')[1])
    thisScript = ''.join((parts[4].split('.')[:-1]))
    ed = float(parts[5])
    #print ye, thisScript, ed

    #silly way to say only plot this guy
    #if not thisScript == script:
    if not ed == edToPlot:
        #print "'%s' is not '%s'" % (script, thisScript)
        continue

    filehandle = open(fileLocation + '/' + file, 'r')
    data = {'r': [], 'rho': [], 'p': [], 'm': []}
    for line in filehandle:
        entry = line.split()
        data['r'].append(float(entry[0]))
        data['rho'].append(float(entry[1]))
        data['p'].append(float(entry[2]))
        data['m'].append(float(entry[3]))
    print thisScript, data['m'][-1]/data['rho'][-1]**3
    plt.plot(data['m'], data['rho'], label=thisScript)
plt.legend()
plt.show()