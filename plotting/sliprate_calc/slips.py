# Python script to determine histograms of model slip rate populations
# Written by T.J. Craig, University of Leeds, August 2018
#
# To run in foe-linux-02, start with command 
# app setup canopy python-libs
# To run type
# python slips.py

# Print starting time:
import datetime
print (datetime.datetime.now())
startruntime = (datetime.datetime.now())


# Requires pre-installed packages
import numpy as np
import math
import copy
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import cm
import pylab

# Fix issues with font exporting to postscript and pdf files
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Close any existing plots that are open
plt.close('all')

# Set number of events, including 1 cm trench
nevs = 9

# Set number of models
nmods = 10000

# Set input filename
infile = open('slips_CAP.inp','r')
inputfile=infile.read().splitlines()

# Initialise models arrays
slip_dates = np.zeros((nmods,nevs))

# Read in model slip incremements
slipparams=map(float, inputfile[0].split())
slip_magnitudes = [0.0] * nevs
for i in range(0,nevs):
    slip_magnitudes[i] = slipparams[i]

# Read in model data

for i in range(0,nmods):
    if (i%100 == 0):
        print "Reading in model",i
    slips=map(float, inputfile[(i+1)].split())
    for j in range(0,nevs):
        slip_dates[i][j] = slips[j]
        
        
# Define slip rate and time grid
t1=0
t2=10000
tinc=2
nt=int(round((t2-t1)/tinc))

s1=0
s2=20.0
sinc=0.05
ns=int(round((s2-s1)/sinc))

counts = np.zeros((nt,ns))
avslip = [0.0] * nt
times = [0.0] * nt
discarded_up = [0] * nt
discarded_down= [0] * nt
discarded_fast= [0] * nt
negslip = [0] * nt

for i in range(0,nt):
    if (i%100 == 0):
        print "Doing slip rate calculation",i
    for k in range(0,nmods):
        runner = 0
        for l in range(0,nevs):
            if (((i*tinc)+t1) < slip_dates[k][l]):
                runner += 1
            else:
                pass

        if (runner == 0):
            discarded_up[i] +=1
        elif (runner == nevs):
            discarded_down[i] +=1 
        else:
            slip_rate = (10*slip_magnitudes[runner] / ((slip_dates[k][(runner - 1)] - slip_dates[k][runner])))

            # Discretize slip rate into increments of sinc
	    if (slip_rate < 0):
                negslip[i] +=1
            elif (slip_rate < (s2-sinc)):
                slip_rate = int(round(slip_rate / sinc))
                counts[i][slip_rate] += 1        
            else:
                discarded_fast[i] += 1 

# Open output file for amplitudes
ofile = open('output_slip_averages.dat', 'w')
for i in range(0,nt):
    times[i] = i * tinc
    counter = 0
    slip_sum = 0
    for j in range(0,ns):
        counter = counter + counts[i][j]
        slip_sum = slip_sum + counts[i][j] * ((j * sinc) + s1)
    avslip[i] = slip_sum / counter
    otext = str(times[i]) + ' ' + str(avslip[i]) + ' ' + str(slip_sum) + ' ' + str(counter) + ' ' + str(discarded_fast[i]) + ' ' + str(discarded_up[i]) + ' ' + str(discarded_down[i]) + ' ' + str(negslip[i]) + '\n'
    ofile.write(otext)
ofile.close()

plt.figure(1, figsize=(9,10.5))
plt.title('Slip rate density',fontsize=10)
counts = np.fliplr(counts)
plt.imshow(counts.T,cmap='gist_heat_r',extent=[t1,t2,s1,s2],aspect='auto',vmin=0,vmax=(0.075*counts.max()))
plt.ticklabel_format(style='plain')
plt.ylabel('Slip Rate (mm/yr)',fontsize=10)  
plt.xlabel('Time B.P. (yrs)',fontsize=10)  
plt.tick_params(axis='y',labelleft='on',labelright='off',left='on',right='off',labelsize=8)
plt.tick_params(axis='x',labelbottom='on',labeltop='off',top='on',labelsize=8)        
pylab.savefig('Slip_Rates.eps', format='ps',rasterized=True,dpi=900, papertype='a1')
pylab.savefig('Slip_rates.png', bbox_inches='tight', dpi=300)
plt.show()
        
plt.figure(2, figsize=(9,10.5))
plt.title('Slip rate density',fontsize=10)
counts = np.fliplr(counts)
plt.imshow(counts.T,cmap='gist_heat_r',extent=[t1,t2,s1,s2],aspect='auto',vmin=0,vmax=(0.075*counts.max()))
plt.plot(times,avslip,color='black',linewidth=1)
plt.ticklabel_format(style='plain')
plt.ylabel('Slip Rate (mm/yr)',fontsize=10)  
plt.xlabel('Time B.P. (yrs)',fontsize=10)  
plt.tick_params(axis='y',labelleft='on',labelright='off',left='on',right='off',labelsize=8)
plt.tick_params(axis='x',labelbottom='on',labeltop='off',top='on',labelsize=8)        
pylab.savefig('Slip_Rates_Av.eps', format='ps',rasterized=True,dpi=900, papertype='a1')
pylab.savefig('Slip_rates_Av.png', bbox_inches='tight', dpi=300)
plt.show()
        
        
# Print ending time
print(datetime.datetime.now())
print (datetime.datetime.now() - startruntime)
