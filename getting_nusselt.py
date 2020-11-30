''' Waiting for Christian to reply to this, with code amendment or confirmation
	that one of these methods is the right one for obtaining the Nusselt Number. '''

# System Imports
import os
import sys
import h5py
import numpy as np
import pathlib

# Dedalus Specific Imports
from dedalus import public as de
import logging
logger = logging.getLogger(__name__)

# Imports from our own scripts
from run_param_file import Lz, Pr

Prandtl = Pr
z_domain = Lz
this_dir = os.getcwd()
results_direc = "Results3j" #Change the file to the one we are investigating
save_direc = os.path.join(results_direc, "figs1")

# Change these dependent on what the time values were for the flux.
avg_t_start = 2.5
avg_t_stop  = 6.0

analysis_file = "analysis.h5"
snapshots_file = "snapshots.h5"
direc = "raw_data2/"

with h5py.File(results_direc + "/" + direc + analysis_file, mode='r') as file:
	# Why do we index the x about 0th position - NB this actually gets all the data. We have averaged over x in Dedalus already
    L_cond_all = np.array(file['tasks']['L_cond'])[:,0,:]
    L_conv_all = np.array(file['tasks']['L_conv'])[:,0,:]
    ana_t = np.array(file['scales']['sim_time'])
    #print(ana_t)
    # Add arrays, then take the mean. And compare to take the mean, then add and divide. 
print(len(L_cond_all))
print(len(L_conv_all))

#### 1. Adding the arrays, then taking the mean for a Nusselt number #### 
nusselt1 = (L_cond_all + L_conv_all) / L_cond_all
print(f"Nusselt by adding arrays then averaging is : {nusselt1.mean()}")
print(f"Nusselt by adding arrays then taking absolute of average is : {np.absolute(nusselt1.mean())}")
# Result is : -3.9870457754552433 - take absolute and this becomes positive?




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

#### 2. Averaging the arrays over the selct time period then calculating the Nusselt Number #### 
mean_L_cond = np.mean(np.array(L_cond_all[ASI:AEI,:]), axis=0)
mean_L_conv = np.mean(np.array(L_conv_all[ASI:AEI,:]), axis=0)#The axis=0 takes the mean about the 0th axis, or time in this case. 

nusselt_2 = (mean_L_cond + mean_L_conv) / mean_L_cond
nusselt_3 = np.absolute((mean_L_cond + mean_L_conv) / mean_L_cond)
print(f"Nusselt by averaging then adding arrays is : {nusselt_2.mean()}")
print(f"Nusselt by averaging then taking absolute of added arrays is : {nusselt_3.mean()}")
# Result is : -4.981868383013649 but taking the aboslute makes this 25.855340643870626.


''' This plot results in figs1fig2.png. Which shows how Nusselt changes over z. Nusselt on x-axis.'''
import matplotlib.pyplot as plt

z_basis = de.Chebyshev('z', 64, interval=(0,1), dealias=3/2)
z = np.array(z_basis.grid(1))

mean_L_cond = np.mean(np.array(L_cond_all[ASI:AEI,:]), axis=0)
mean_L_conv = np.mean(np.array(L_conv_all[ASI:AEI,:]), axis=0)
plt.plot(((mean_L_cond+mean_L_conv)/mean_L_cond),z, 'b', linestyle='-', label="Nusselt against z")
plt.savefig(save_direc + 'fig4')
plt.clf()
plt.close()




'''
#### 3. mean_flux = (T*w*Pr).mean()	####
with h5py.File(results_direc + "/" + direc + snapshots_file, mode='r') as file:
    u_all = np.array(file['tasks']['u'])
    w_all = np.array(file['tasks']['w'])
    T_all = np.array(file['tasks']['T'])
mean_flux = (T_all*w_all*Prandtl).mean()
print(f"Mean flux is {mean_flux}")
flux_1 = np.array(T_all*w_all*Prandtl)[:,0,:]
flux_2 = np.array(T_all*w_all*Prandtl)[660:930,0,:]
flux_3 = np.array(T_all*w_all*Prandtl)[750:930,0,:]

#NB Iteration 660:930 are the iterations where the KE is constant. Hence why chosen.
print(f"Flux 1 - in array is {flux_1.mean()}")#2.5930823819515454
print(f"Flux 2 - in array is {flux_2.mean()}")#4.913245893304888
print(f"Flux 3 - in array is {flux_3.mean()}")#5.077343264141911

'''
