'''
Horizontally averaged T as a function of Height
T indexed against x, take average then plot against z
'''
# System Imports
import os
import sys
import h5py
import numpy as np
import pathlib
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm

# Dedalus Specific Imports
from dedalus import public as de
import logging
logger = logging.getLogger(__name__)

# Imports from our own scripts
from run_param_file import Lz, Pr, Ra 

Rayleigh = Ra
Prandtl = Pr
z_domain = Lz
this_dir = os.getcwd()
results_direc = "Results3h" #Change the file to the one we are investigating
save_direc = os.path.join(results_direc, "figs1_")

analysis_file = "analysis.h5"
snapshots_file = "snapshots.h5"
direc = "raw_data2/"

with h5py.File(results_direc + "/" + direc + analysis_file, mode='r') as file:
	# Why do we index the x about 0th position - NB this actually gets all the data. We have averaged over x in Dedalus already
    #L_cond_all = np.array(file['tasks']['L_cond'])[:,0,:]
    #L_conv_all = np.array(file['tasks']['L_conv'])[:,0,:]
    KE = np.array(file['tasks']['KE'])[:,0,0]
    ana_t = np.array(file['scales']['sim_time'])

with h5py.File(results_direc + "/" + direc + snapshots_file, mode='r') as file:
    T_all = np.array(file['tasks']['T'])
    snap_iter = np.array(file['scales']['iteration'])


### Calculating the average kinetic energy across the whole time period! ### 
mean_kinetic_energy = KE.mean()
print(f"Average kinetic energy is : {mean_kinetic_energy}")


# Change these dependent on what the time values were for the flux.
avg_t_start = 2.5
avg_t_stop  = 6.0

if avg_t_start <= ana_t[0] or avg_t_stop <= ana_t[0]:
    sys.exit("Average time period out of simulation range: {} -> {}".format(ana_t[0], ana_t[-1]))
if avg_t_start >= ana_t[-1] or avg_t_stop >= ana_t[-1]:
    sys.exit("Average time period out of simulation range: {} -> {}".format(ana_t[0], ana_t[-1]))

ASI = (np.abs(ana_t  - avg_t_start)).argmin()  # analysis start index
if np.isnan(avg_t_stop): # End of array if NaN value given
    AEI = -1,
else:
    AEI = (np.abs(ana_t  - avg_t_stop)).argmin()   # analysis end index
avg_t_range = ana_t[AEI] - ana_t[ASI]



z_basis = de.Chebyshev('z', 64, interval=(0,1), dealias=3/2)
z = np.array(z_basis.grid(1))

# Getting the horizontally averaged temperature profile
horizontal_average_T = np.mean(T_all, axis=1)
#print(horizontal_average_T)
time_horizontal_average_T = np.mean(horizontal_average_T, axis=0)
#print(time_horizontal_average_T)

#mean_L_cond = np.mean(np.array(L_cond_all[ASI:AEI,:]), axis=0)
#mean_L_conv = np.mean(np.array(L_conv_all[ASI:AEI,:]), axis=0)

c_m = matplotlib.cm.OrRd

s_m = matplotlib.cm.ScalarMappable(cmap=c_m)
s_m.set_array([])

plt.plot(z, time_horizontal_average_T, 'black', linestyle='-', label="Temperature Profile")
plt.colorbar(s_m)
plt.ylabel("Temperature (averaged in x and time)")
plt.xlabel("z")
plt.title(f"Temperature profile against z : Ra = {Rayleigh}")
plt.legend()
plt.savefig(save_direc + 'temperature_profile_Ra')
plt.clf()
plt.close()

plt.plot(time_horizontal_average_T, z, 'black', linestyle='-', label="Temperature Profile")
############# If you dont want the colorbar remove the following - single - line!!! #############
#plt.colorbar(s_m)
plt.xlabel("Temperature (averaged in x and time)")
plt.ylabel("z")
plt.title(f"Temperature profile against z : Ra = {Rayleigh}")
plt.legend()
plt.savefig(save_direc + 'temperature_profile_Ra_flip')
plt.clf()
plt.close()
