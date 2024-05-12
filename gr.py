#/* --------------------------------------------------------------------------------
#   Egorov Group
#   Mohan Shankar, mjs7eek@virginia.edu

#   gr.py
#   "This file reads in g(r) data directly from lammpstrj files and plots them"
#-------------------------------------------------------------------------------- */
# DEPENDENCIES 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import pandas as pd

#-------------------------------------------------------------------------------- */
# NOT SURE WHY THESE ARE NEEDED, BUT THEY ARE
# might have to run `pip install PyQt6` or `pip3 install PyQt6`
import sys
from PyQt6.QtWidgets import QApplication
app = QApplication(sys.argv)
#-------------------------------------------------------------------------------- */
# OVITO IMPORTS 

from ovito.io import import_file
from ovito.modifiers import CoordinationAnalysisModifier, TimeAveragingModifier
#-------------------------------------------------------------------------------- */
# READ IN DATA 

# Load a particle dataset, apply the modifier, and evaluate pipeline.
pipeline = import_file("/Users/mohan/Desktop/Research/likos/p7.5/dump.lammpstrj") # PATH to file


#-------------------------------------------------------------------------------- */
# PERFORM g(r) CALCULATIONS 
pipeline.modifiers.append(CoordinationAnalysisModifier(cutoff = 4.5, number_of_bins = 200, partial = True))

# Insert the time-averaging modifier into the pipeline, which accumulates
# the instantaneous DataTable produced by the previous modifier and computes a mean histogram.
pipeline.modifiers.append(TimeAveragingModifier(operate_on='table:coordination-rdf'))

# Data export method 1: Convert to NumPy array and write data to a text file:
total_rdf = pipeline.compute().tables['coordination-rdf[average]'].xy()

# Access the output DataTable:
rdf_table = pipeline.compute().tables['coordination-rdf']

# The y-property of the data points of the DataTable is now a vectorial property.
# Each vector component represents one partial RDF.
rdf_names = rdf_table.y.component_names

#-------------------------------------------------------------------------------- */
# SAVE DATA IN RELEVANT VARIABLES 
r = total_rdf[:, 0]
one_one = total_rdf[:, 1]
one_two = total_rdf[:, 2]
two_two = total_rdf[:, 3]

# np.savetxt('gr.txt', np.transpose((r, one_one, one_two, two_two)), delimiter = ',', header="r, one-one, one-two, two-two", fmt='%1.5e') # uncomment to save data

#-------------------------------------------------------------------------------- */
# UNCOMMENT THIS SECTION IF YOU'RE READING IN g(r) DATA

# df = pd.read_csv("PATH/gr.txt", sep = ',') # use this line of code to read in saved data if needed

# r_s = df.iloc[:, 0] # the first column is r
# one_one_s = df.iloc[:, 1] # one-one interaction
# one_two_s = df.iloc[:, 2] # one-two interaction
# two_two_s = df.iloc[:, 3] # two-two interaction
#-------------------------------------------------------------------------------- */
# PLOT g(r)

fig, ax = plt.subplots(3, figsize=(12, 7))

ax[0].plot(r, one_one, color = 'RoyalBlue', label = 'OVITO') # plot 1-1 interaction g(r)
ax[1].plot(r, one_two, color = 'r', label = 'OVITO') # plot 1-2 interaction g(r)
ax[2].plot(r, two_two, color = 'k', label = 'OVITO') # plot 2-2 interaction g(r)


for j in range(3):
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
    ax[j].legend()

ax[0].set_ylabel(r'$\bf{g_{11}(r)}$')
ax[0].axhline(y = 1.0, color = 'k')

ax[1].set_ylabel(r'$\bf{g_{12}(r)}$')
ax[1].axhline(y = 1.0, color = 'k')

ax[2].set_ylabel(r'$\bf{g_{22}(r)}$')
ax[2].axhline(y = 1.0, color = 'k')

fig.supxlabel(r'$\bf{r / \sigma}$')
plt.show();
#-------------------------------------------------------------------------------- */