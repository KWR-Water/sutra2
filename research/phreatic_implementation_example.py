#%% ----------------------------------------------------------------------------
# A. Hockin, March 2021
# KWR BO 402045-247
# ZZS verwijdering bodempassage
# AquaPriori - Transport Model
# With Martin Korevaar, Martin vd Schans, Steven Ros
#
# ------------------------------------------------------------------------------

#### Notes ####
# Example for use in the read the docs 
####

#%% ----------------------------------------------------------------------------
# INITIALISATION OF PYTHON e.g. packages, etc.
# ------------------------------------------------------------------------------

# %reset -f #reset all variables for each run, -f 'forces' reset, !! 
# only seems to work in Python command window...

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
# from pandas import read_excel
from pandas import read_csv
from pandas import read_excel
# import pyarrow.parquet as pq
import math
from scipy.special import kn as besselk

from pathlib import Path
try:
    from project_path import module_path #the dot says look in the current folder, this project_path.py file must be in the folder here
except ModuleNotFoundError:
    from project_path import module_path

from greta.Analytical_Well import *
from greta.Substance_Transport import *
from testing.test_transatomic import *
# get directory of this file
path = Path(__file__).parent #os.getcwd() #path of working directory

# path = os.getcwd()  # path of working directory

# %%
phreatic_schematisation = HydroChemicalSchematisation(schematisation_type='phreatic',
                                      well_discharge=7500, #m3/day
                                      ground_surface = 22.,
                                      redox_vadose_zone='anoxic',
                                      redox_shallow_aquifer='anoxic',
                                      redox_target_aquifer='deeply_anoxic',
                                      pH_target_aquifer=7.,
                                      temperature=11.,
                                      diffuse_input_concentration = 100, #ug/L
                                      )

print(phreatic_schematisation.schematisation_type)
print(phreatic_schematisation.well_discharge)
print(phreatic_schematisation.porosity_shallow_aquifer)

#%%

phreatic_well = AnalyticalWell(phreatic_schematisation)
phreatic_well.phreatic() 

radial_plot = phreatic_well.plot_travel_time_versus_radial_distance(xlim=[0, 2000] )#, ylim=[1e3, 1e6])

cumulative_plot = phreatic_well.plot_travel_time_versus_cumulative_abstracted_water(xlim=[0, 1])#, ylim=[1e3, 1e6])

# radial_plot
# cumulative_plot
#%%

# test_substance = Substance(substance_name='benzene')

phreatic_concentration = SubstanceTransport(phreatic_well, substance = 'OMP-X')
phreatic_concentration.compute_omp_removal()
phreatic_concentration.plot_concentration(ylim=[0,100 ])

#%%

phreatic_well = AnalyticalWell(phreatic_schematisation)
phreatic_well.phreatic() 
phreatic_concentration = SubstanceTransport(phreatic_well, substance = 'benzene')
phreatic_concentration.compute_omp_removal()
phreatic_concentration.plot_concentration(ylim=[0,10 ])

phreatic_well = AnalyticalWell(phreatic_schematisation)
phreatic_well.phreatic() 
phreatic_concentration = SubstanceTransport(phreatic_well, substance = 'AMPA')
phreatic_concentration.compute_omp_removal()
phreatic_concentration.plot_concentration( ylim=[0,1 ])
# %%
phreatic_well.plot_travel_time_versus_radial_distance(xlim=[0, 2000], ylim=[1e3, 1e6])

phreatic_well.plot_travel_time_versus_cumulative_abstracted_water(xlim=[0, 1], ylim=[1e3, 1e6])

# %%


time_array = np.arange(0, 505, 1)*365.24

phreatic_concentration.df_flowline['well_conc'] = phreatic_concentration.df_flowline['breakthrough_concentration'] * phreatic_concentration.df_flowline['discharge']/ phreatic_concentration.df_flowline['well_discharge']
df = phreatic_concentration.df_flowline
well_conc = []
for i in range(len(time_array)):
    t = time_array[i]
    well_conc.append(sum(df['well_conc'].loc[df['total_breakthrough_travel_time'] <= t]))

plt.plot(time_array/365.24, well_conc)



# %%
