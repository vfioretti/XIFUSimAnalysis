"""
 XIFU_TESdep_raw.py  -  description
 ---------------------------------------------------------------------------------
 Plotting the output of the THELSim simulation on TES and create the ASCII file
 ---------------------------------------------------------------------------------
 copyright            : (C) 2015 Valentina Fioretti
 email                : fioretti@iasfbo.inaf.it
 ----------------------------------------------
 Usage:
 python XIFU_TESdep_raw.py filedir N_file N_in flux "title" N_bin part_flag name_out
 ---------------------------------------------------------------------------------
 Parameters:
 - filedir = input path (string)
 - N_file = number of simulated files
 - flux = flux in part. cm-2 s-1 sr-1
 - N_in = number of simulated particles
 - title = the title, in brackets
 - N_bin = number of bins
 - part_flag = type of input particle [0 = proton]
 - name_out = filename of the ASCII output
 --------------------------------------------------------------------------------
 Caveats:
 None
 ---------------------------------------------------------------------------------
 Modification history:
 - 2015/02/12: creation date
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
name_out = arg_list[8]

# XIFU settings
N_pix_X = 104
pix_side = 0.03 #cm
TES_side = pix_side*N_pix_X
e_min = 0.3 # keV
e_max = 12. # keV


# set-up
vecEnergy = []
vecEventID = []
vecTESID = []
vecVolID = []
vecMothID = []

for jfits in xrange(N_fits):
    
    hdulist = pyfits.open(filedir+'xyz.'+str(jfits)+'.fits.gz')
    
    tbdata = hdulist[1].data
    cols = hdulist[1].columns
    
    evt_id = tbdata.field('EVT_ID')
    e_dep = tbdata.field('E_DEP')
    vol_name = tbdata.field('VOLUME_NAME')
    vol_id = tbdata.field('VOLUME_ID')
    moth_id = tbdata.field('MOTHER_ID')
    
    where_tes = np.where(vol_name == 'pxlX')
    evt_id = evt_id[where_tes]
    e_dep = e_dep[where_tes]
    vol_id = vol_id[where_tes]
    moth_id = moth_id[where_tes]
 
    where_dep = np.where(e_dep > 0.)
    evt_id = evt_id[where_dep]
    e_dep = e_dep[where_dep]
    vol_id = vol_id[where_dep]
    moth_id = moth_id[where_dep]
    
    tes_id = np.zeros(len(vol_id))
    for jev in xrange(len(vol_id)):
        tes_id[jev] = moth_id[jev]*N_pix_X + vol_id[jev]

    tes_id_tot = []
    vol_id_tot = []
    moth_id_tot = []
    e_dep_tot = []
    evt_id_tot = []

    index = 0
    while 1:
        where_sameevent = np.where(evt_id == evt_id[index])
        where_sameevent = where_sameevent[0]
        
        temp_tes_id = tes_id[where_sameevent]
        temp_vol_id = vol_id[where_sameevent]
        temp_moth_id = moth_id[where_sameevent]
        temp_e_dep = e_dep[where_sameevent]

        print evt_id[index]

        vol_index = 0
        while 1:
            where_samevol = np.where(temp_tes_id == temp_tes_id[vol_index])
            where_samevol = where_samevol[0]
            print 'temp_tes_id', temp_tes_id
            print len(temp_tes_id)
            print where_samevol
            
            tes_id_tot.append(temp_tes_id[vol_index])
            vol_id_tot.append(temp_vol_id[vol_index])
            moth_id_tot.append(temp_moth_id[vol_index])
            e_dep_tot.append(np.sum(temp_e_dep[where_samevol]))
            evt_id_tot.append(evt_id[index])

            if (len(where_samevol) < len(temp_tes_id)):
                temp_tes_id = np.delete(temp_tes_id, where_samevol)
                temp_vol_id = np.delete(temp_vol_id, where_samevol)
                temp_moth_id = np.delete(temp_moth_id, where_samevol)
                temp_e_dep = np.delete(temp_e_dep, where_samevol)
                print 'temp_tes_id after', temp_tes_id
            else:
                break


        N_event_eq = len(where_sameevent)
        print 'N_event_eq', N_event_eq
        print 'where_sameevent', where_sameevent
        print 'where_sameevent[N_event_eq-1]', where_sameevent[N_event_eq-1]

        if (where_sameevent[N_event_eq-1] < (len(evt_id)-1)):
            index = where_sameevent[N_event_eq-1] + 1
        else:
            break

    for jev in xrange(len(evt_id_tot)):
            vecEnergy.append(e_dep_tot[jev])
            vecEventID.append(evt_id_tot[jev])
            vecTESID.append(tes_id_tot[jev])
            vecVolID.append(vol_id_tot[jev])
            vecMothID.append(moth_id_tot[jev])



# Select the energy range
vecEnergy = np.array(vecEnergy)
vecEventID = np.array(vecEventID)
vecTESID = np.array(vecTESID)
vecVolID = np.array(vecVolID)
vecMothID = np.array(vecMothID)

where_range = np.where((vecEnergy >= e_min) & (vecEnergy <= e_max))
where_range = where_range[0]

vecEnergy = vecEnergy[where_range]
vecEventID = vecEventID[where_range]
vecTESID = vecTESID[where_range]
vecVolID = vecVolID[where_range]
vecMothID = vecMothID[where_range]

n_bins = np.logspace(np.log10(e_min), np.log10(e_max), N_bin)


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

# detection area
A_det = TES_side**2.  # cm2


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
    F_array[jn] = (float(N_array[jn]))/(energy_bin*Time_new*A_det)
    err_F_array[jn] = (float(err_N_array[jn]))/(energy_bin*Time_new*A_det)


# Open and write file
f = open(name_out, 'wb')

f.write("# Event_ID PixelX_ID PixelY_ID GlobalVol_ID Energy\n")
for jline in xrange(vecEventID.size):
    f.write(str(int(vecEventID[jline]))+" ")
    f.write(str(int(vecVolID[jline]))+" ")
    f.write(str(int(vecMothID[jline]))+" ")
    f.write(str(int(vecTESID[jline]))+" ")
    f.write(str(vecEnergy[jline])+" ")
    f.write("\n")


f.close()


# Plot the results
fig = plt.figure(1, figsize=(12, 8))
ax = fig.add_subplot(111)

# simulation
ax.errorbar(ene_array, F_array, xerr=err_ene, yerr=err_F_array, fmt='.k', lw = 2, ecolor='gray', label='Simulation')

ax.set_title(title_label)
ax.set_xlim([e_min, e_max])
ax.set_xlabel("Energy [keV]")
ax.set_ylabel("cts cm$^{-2}$ s$^{-1}$keV$^{-1}$")
ax.set_xscale('log')
ax.set_yscale('log')


ax.legend(numpoints=1)

fig2 = plt.figure(2, figsize=(12, 8))
ax2 = fig2.add_subplot(111)

# simulation
ax2.plot(ene_array, N_array, '-r', lw = 2, label='Simulation')

ax2.set_title(title_label)
ax2.set_xlim([e_min, e_max])
ax2.set_xlabel("Energy [keV]")
ax2.set_ylabel("counts")
ax2.set_xscale('log')
ax2.set_yscale('log')



ax2.legend(numpoints=1)
plt.show()


hdulist.close()
