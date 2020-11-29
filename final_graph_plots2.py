import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm

this_dir = os.getcwd()
results_direc = "Results3i" #Change the file to the one we are investigating
save_direc = os.path.join(results_direc, "Rayleigh")

Ra = [
2e2,
5e2,
8e2,
1e3,
1.6e3,
2.2e3,
2.8e3,
3.8e3,
7e3,
1e4,
]


nusselt = [
1,
1.000000000005889,
1.059833654478345,
1.2434958055325742,
2.015390760770154,
2.754803215041899,
3.7706679543466395,
13.282422035135092,
14.723061513660493,
17.084273092673573,
]


average_KE = [
6.746461722281338e-16,
5.403527770505897e-10,
2.3693330251791784,
6.4391982472641835,
18.668431821773368,
30.043497315976545,
40.23368791576478,
54.554916013541735,
113.00751832502667,
178.34889508292403,
]

average_velocity = []
for value in average_KE:
	velocity = np.sqrt(value)
	average_velocity.append(velocity)

import math

log_nusselt = []
for item in nusselt:
	new_value = log_nusselt.append(math.log(item))

log_Ra = []
for item in Ra:
	new_value = log_Ra.append(math.log(item))

print(log_nusselt)
print(average_velocity)

nusselt_error = [
1.0028,
1.0028,
1.0190,
1.1777,
0.8816,
1.1729,
1.3656,
1.5691,
1.9501,
2.1106,
]

ln_nusselt_error = []
for item in nusselt_error:
	ln_nusselt_error.append(math.log(item))

plt.plot(Ra, nusselt, 'black', linestyle='-')
plt.scatter(Ra, nusselt)
plt.ylabel("Average Nu",fontsize=14)
plt.xlabel("Ra", fontsize=14)
plt.title("Rayleigh vs Nusselt Number",fontsize=18)
#plt.legend()
plt.savefig(this_dir + '_Nusselt')
plt.clf()
plt.close()

''' For line of best fit from the point there is convection '''
m, c = np.polyfit(log_Ra[2:],log_nusselt[2:], 1)
best_fit = []
for item in log_Ra:
	item1 = m*item + c
	best_fit.append(item1)

#plt.plot(log_Ra[2:], log_nusselt[2:], 'black', linestyle='-', label="Nusselt Number")
m, c = np.polyfit(log_Ra[2:],log_nusselt[2:], 1)
plt.plot(log_Ra[2:], best_fit[2:], 'black')
plt.plot(log_Ra[:3], [0,0,0], 'black')
plt.scatter(log_Ra, log_nusselt)
#plt.ylim(-0.1,3.5)
plt.errorbar(log_Ra, log_nusselt, yerr=ln_nusselt_error, ls="none")
plt.ylabel("ln(Nu)",fontsize=14)
plt.xlabel("ln(Ra)",fontsize=14)
plt.title(f"ln(Rayleigh) vs ln(Nusselt) Number",fontsize=18)
#plt.legend()
plt.savefig(this_dir + '_Nusselt_log')
plt.clf()
plt.close()


m1, c1 = np.polyfit(log_Ra[1:], average_velocity[1:], 1)
best_fit_velocity = []
for item in log_Ra:
	item1 = m1*item + c1
	best_fit_velocity.append(item1)


plt.plot(Ra, average_velocity, 'black', linestyle='-')
y=x**1/3
plt.plot(y)
plt.scatter(Ra, average_velocity)
plt.ylabel("Average Velocity",fontsize=14)
plt.xlabel("Ra",fontsize=14)
plt.title(f"Rayleigh vs Average Velocity",fontsize=18)
#plt.legend()
plt.savefig(this_dir + '_average_velocity')
plt.clf()
plt.close()


plt.plot(log_Ra[1:], best_fit_velocity[1:], 'black', linestyle="-")
plt.plot(log_Ra[:2], [0,0], 'black')
plt.scatter(log_Ra, average_velocity)
plt.ylabel("Average Velocity",fontsize=14)
plt.xlabel("ln(Ra)",fontsize=14)
plt.ylim(-0.5,15)
plt.title(f"ln(Rayleigh) vs Average Velocity",fontsize=18)
#plt.legend()
plt.savefig(this_dir + '_average_velocity_log')
plt.clf()
plt.close()
