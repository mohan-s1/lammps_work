#/* -----------------------------------------------------------------------------------------------
#   Egorov Group
#   Mohan Shankar, mjs7eek@virginia.edu

#   thermo_data.py
#   "This file reads in LAMMPS logfiles (and anything else) to grab relevant thermodynamic data"
#----------------------------------------------------------------------------------------------- */
# DEPENDENCIES 
import lammps_logfile
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
#----------------------------------------------------------------------------------------------- */
# IMPORT DATA 
log = lammps_logfile.File("PATH/log.lammps") # read in LAMMPS log file

#----------------------------------------------------------------------------------------------- */
# SAVE VARIABLES FROM LAMMPS LOGFILE
t = log.get("Step") # save the step data as `t`
pres = log.get("Press") # save the pressure data as `pres`

#----------------------------------------------------------------------------------------------- */
# SAVE VARIABLES FROM .TXT FILE

mu = np.loadtxt('PATH', usecols = 0) # save first col from txt file as mu 

avg_pres = np.mean(pres) # average the pressure array

avg_mu = np.mean(mu) # average the chemical potential array

density = 3.5 # we specify density 

temp = 0.5 # we specify the reduced temperature 

free_energy = (avg_mu - (avg_pres/density)) / (1/temp) # equation for Helmholtz free energy

print("Free Energy", free_energy)

#----------------------------------------------------------------------------------------------- */
# PLOT 
fig, ax = plt.subplots(2, 1)

ax[0].plot(t, mu, color = 'royalblue')
ax[0].set_ylabel(r'$\mu / \varepsilon$')

ax[1].plot(t, pres, color = 'k')
ax[1].set_ylabel(r'$P \sigma^3 / \varepsilon$')
ax[1].set_ylim(71, 74)

for j in range(2):
    ax[j].tick_params(bottom=True, top=True, left=True, right=True)
    ax[j].tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)

    ax[j].tick_params(axis="x",direction = "in",length=5)
    ax[j].tick_params(axis="y",direction="in",length=5)

    for axis in ['top','bottom','left','right']:
        ax[j].spines[axis].set_linewidth(1.4)

    ax[j].tick_params(width=1.4)
    ax[j].xaxis.set_minor_locator(AutoMinorLocator())
    ax[j].tick_params(top=True, which='minor', length=3,direction='in',width = 1.4)
    ax[j].yaxis.set_minor_locator(AutoMinorLocator())
    ax[j].tick_params(left=True,right = True, which='minor', length=3,direction='in',width = 1.4)

fig.supxlabel(r'$Step$')
plt.show();
#----------------------------------------------------------------------------------------------- */
print("Average Pressure: {}".format(avg_pres))
print("Average Chemical Potential: {}".format(avg_mu))
#----------------------------------------------------------------------------------------------- */