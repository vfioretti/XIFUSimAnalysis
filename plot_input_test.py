"""
 plot_input_test.py  -  description
 ---------------------------------------------------------------------------------
 Plotting the input flux of the THELSim simulation
 ---------------------------------------------------------------------------------
 copyright            : (C) 2015 Valentina Fioretti
 email                : fioretti@iasfbo.inaf.it
 ----------------------------------------------
 Usage:
 python plot_input_test.py filedir N_file N_in flux "title" N_bin part_flag
 ---------------------------------------------------------------------------------
 Parameters:
 - filedir = input path (string)
 - N_file = number of simulated files
 - flux = flux in part. cm-2 s-1 sr-1
 - N_in = number of simulated particles
 - title = the title, in brackets
 - N_bin = number of bins
 - part_flag = type of input particle [0 = proton]
 --------------------------------------------------------------------------------
 Caveats:
 None
 ---------------------------------------------------------------------------------
 Modification history:
 - 2015/02/05: creation date
"""

 

import pyfits
import numpy as np
import math
import sys
import matplotlib.pyplot as plt

# Import the input parameters
arg_list = sys.argv
filedir = arg_list[1]
N_fits = int(arg_list[2])
N_in = float(arg_list[3])
flux_in = float(arg_list[4])
title_label = arg_list[5]
N_bin = float(arg_list[6])
part_flag = int(arg_list[7])

vecEnergy = []
vecEventID = []

for jfits in xrange(N_fits):

    hdulist = pyfits.open(filedir+'xyz.'+str(jfits)+'.fits.gz')

    tbdata = hdulist[1].data
    cols = hdulist[1].columns

    evt_id = tbdata.field('EVT_ID')
    e_in = tbdata.field('E_KIN_ENT')
    vol_name = tbdata.field('VOLUME_NAME')
    
    where_sphere = np.where(vol_name == 'sphere1')
    evt_id = evt_id[where_sphere]
    e_in = e_in[where_sphere]
    
    for jev in xrange(len(evt_id)):
        if jev == 0:
            vecEnergy.append(e_in[jev]/1000.) # from keV to MeV
            vecEventID.append(evt_id[jev])
        else:
            if (evt_id[jev-1] != evt_id[jev]):
                vecEnergy.append(e_in[jev]/1000.) # from keV to MeV
                vecEventID.append(evt_id[jev])



e_hit_min = np.min(vecEnergy)
e_hit_max = np.max(vecEnergy)

n_bins = np.logspace(np.log10(e_hit_min), np.log10(e_hit_max), N_bin)


# Computing the exposure

# Sphere radius
R = 61.1   # cm
# Halfe angle of aperture
q = math.pi/2.   # rad

# Normalization
A_sphere = 4.*math.pi*(R**2.)  # cm2
print 'A_sphere: ', A_sphere
omega = math.pi*((math.sin(q))**2.)   # sr
print 'omega: ', omega

r_min = math.tan(q)*R
A_min = 4.*math.pi*(r_min**2.)
F_tot = flux_in*A_sphere*omega  # part./s
print 'F_tot: ', F_tot

Time_new = N_in/F_tot   # s
print 'Time_new: ', Time_new

# histogram
N_array, bin_array = np.histogram(vecEnergy, bins = n_bins)
err_N_array = np.zeros(len(N_array))
F_array = np.zeros(len(N_array))
err_F_array = np.zeros(len(N_array))
err_ene = np.zeros(len(N_array))
ene_array = np.zeros(len(N_array))

for jn in xrange(len(N_array)):
    err_ene[jn] = (bin_array[jn+1] - bin_array[jn])/2.
    ene_array[jn] = bin_array[jn] + err_ene[jn]
    energy_bin = bin_array[jn+1] - bin_array[jn]
    err_N_array[jn] = math.sqrt(float(N_array[jn]))
    F_array[jn] = (float(N_array[jn]))/(energy_bin*Time_new*A_sphere*math.pi)
    err_F_array[jn] = (float(err_N_array[jn]))/(energy_bin*Time_new*A_sphere*math.pi)


# Open file
#f = open(gcr_file, 'r')

# Read and ignore header lines
#header1 = f.readline()

#vecEnergy_real = []
#vecFlux_real = []
# Loop over lines and extract variables of interest
#for line in f:
#    line = line.strip()
#    columns = line.split()
#    vecEnergy_real.append(float(columns[0]))
#    vecFlux_real.append(float(columns[1])/10000.)

#f.close()

# Create the input model

if (part_flag == 0):
    vecEnergy_real= [1.0000000000E+01, 1.1220184543E+01, 1.2589254118E+01, 1.4125375446E+01,1.5848931925E+01 , 1.7782794100E+01,1.9952623150E+01, 2.2387211386E+01, 2.5118864315E+01,2.8183829313E+01, 3.1622776602E+01,3.5481338923E+01,3.9810717055E+01,4.4668359215E+01,5.0118723363E+01 , 5.6234132519E+01, 6.3095734448E+01, 7.0794578438E+01,7.9432823472E+01,8.9125093813E+01,1.0000000000E+02,1.1220184543E+02,1.2589254118E+02,1.4125375446E+02,1.5848931925E+02,1.7782794100E+02, 1.9952623150E+02,2.2387211386E+02,2.5118864315E+02,2.8183829313E+02, 3.1622776602E+02, 3.5481338923E+02,3.9810717055E+02,4.4668359215E+02,5.0118723363E+02 ,  5.6234132519E+02, 6.3095734448E+02, 7.0794578438E+02, 7.9432823472E+02,8.9125093813E+02 ,  1.0000000000E+03,1.1220184543E+03,  1.2589254118E+03 ,1.4125375446E+03, 1.5848931925E+03, 1.7782794100E+03,1.9952623150E+03,2.2387211386E+03,2.5118864315E+03,2.8183829313E+03 ,3.1622776602E+03, 3.5481338923E+03,3.9810717055E+03,4.4668359215E+03,5.0118723363E+03,5.6234132519E+03,6.3095734448E+03,7.0794578438E+03, 7.9432823472E+03,8.9125093813E+03,1.0000000000E+04 ,1.1220184543E+04, 1.2589254118E+04,1.4125375446E+04, 1.5848931925E+04,1.7782794100E+04,1.9952623150E+04,2.2387211386E+04,2.5118864315E+04,2.8183829313E+04,3.1622776602E+04, 3.5481338923E+04,3.9810717055E+04,4.4668359215E+04,5.0118723363E+04,5.6234132519E+04,6.3095734448E+04,7.0794578438E+04, 7.9432823472E+04 ,  8.9125093813E+04, 1.0000000000E+05]
    
    
    #Proton fux in prot/m2/s/sr/MeV
        
    vecFlux_real = np.zeros(len(vecEnergy_real))
            
    for jen in xrange(len(vecEnergy_real)):
        vecFlux_real[jen] = ((26.63*((vecEnergy_real[jen]/1000.)**(-2.72)))/(1. + ((1.41*((vecEnergy_real[jen]/1000.)**(-3.86)) + 21.08*((vecEnergy_real[jen]/1000.)**(-1.76)))/(1. + 2.9*(10.**(-7))*((vecEnergy_real[jen]/1000.)**(-3.7))))))/10000.



# Plot the results
fig = plt.figure(1, figsize=(12, 8))
ax = fig.add_subplot(111)

# simulation
ax.errorbar(ene_array, F_array, xerr=err_ene, yerr=err_F_array, fmt='.k', lw = 2, ecolor='gray', label='Simulation')

# real
ax.plot(vecEnergy_real, vecFlux_real, '-r', lw = 2, label='Input model')

ax.set_title(title_label)
ax.set_xlabel("Energy [MeV]")
ax.set_ylabel("part. cm$^{-2}$ s$^{-1}$MeV$^{-1}$sr$^{-1}$")
ax.set_xscale('log')
ax.set_yscale('log')


ax.legend(numpoints=1)

plt.show()


hdulist.close()
