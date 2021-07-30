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
from tqdm import tqdm  # tqdm gives a progress bar for the simultation
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
                                      vertical_anistropy_shallow_aquifer = 0.0006,
                                      porosity_vadose_zone=0.38,
                                      porosity_shallow_aquifer=0.35,
                                      porosity_target_aquifer=0.35,
                                      recharge_rate=0.00082, #m/day
                                      moisture_content_vadose_zone=0.15,
                                      ground_surface = 22.,
                                      thickness_vadose_zone_at_boundary=5.,
                                      thickness_shallow_aquifer=10.,
                                      thickness_target_aquifer=40.,
                                      hor_permeability_target_aquifer=35.,
                                      thickness_full_capillary_fringe=0.4,
                                      redox_vadose_zone='anoxic',
                                      redox_shallow_aquifer='anoxic',
                                      redox_target_aquifer='deeply_anoxic',
                                      pH_vadose_zone=5.,
                                      pH_shallow_aquifer=6.,
                                      pH_target_aquifer=7.,
                                      dissolved_organic_carbon_vadose_zone=10., 
                                      dissolved_organic_carbon_shallow_aquifer=4., 
                                      dissolved_organic_carbon_target_aquifer=2.,
                                      fraction_organic_carbon_vadose_zone=0.001,
                                      fraction_organic_carbon_shallow_aquifer=0.0005,
                                      fraction_organic_carbon_target_aquifer=0.0005, 
                                      temperature=11.,
                                      solid_density_vadose_zone= 2.650, 
                                      solid_density_shallow_aquifer= 2.650, 
                                      solid_density_target_aquifer= 2.650, 
                                      diameter_borehole = 0.75,
                                      diffuse_input_concentration = 100, #ug/L
                                    #   halflife_anoxic= 650,
                                      )

phreatic_well = AnalyticalWell(phreatic_schematisation)
phreatic_well.phreatic() 

# test_substance = Substance(substance_name='benzene')

phreatic_concentration = SubstanceTransport(phreatic_well, substance = 'OMP-X')
phreatic_concentration.compute_omp_removal()
phreatic_concentration.plot_concentration(ylim=[0,100 ])

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
