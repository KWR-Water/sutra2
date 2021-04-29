#%% ----------------------------------------------------------------------------
# A. Hockin, March 2021
# KWR BO 402045-247
# ZZS verwijdering bodempassage
# AquaPriori - Transport Model
# With Martin Korevaar, Martin vd Schans, Steven Ros
#
# ------------------------------------------------------------------------------

#### CHANGE LOG ####

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
from draft_transport_function import *


from pathlib import Path
try:
    from project_path import module_path #the dot says looik in the current folder, this project_path.py file must be in the folder here
except ModuleNotFoundError:
    from project_path import module_path

from testing.test_transatomic import *
# get directory of this file
path = Path(__file__).parent #os.getcwd() #path of working directory


# path = os.getcwd()  # path of working directory

#%%
# Test
# test_travel_time_distribution_phreatic()

# %%
scheme1 = HydroChemicalSchematisation(schematisation='phreatic',
                                      well_discharge_m3hour=319.4,
                                      vertical_resistance_aquitard=500,
                                      porosity_vadose_zone=0.38,
                                      porosity_shallow_aquifer=0.35,
                                      porosity_target_aquifer=0.35,
                                      recharge_rate=0.3,
                                      soil_moisture_content_vadose_zone=0.15,
                                      thickness_vadose_zone=5,
                                      thickness_shallow_aquifer=10,
                                      thickness_target_aquifer=40,
                                      KD=1400,
                                      thickness_full_capillary_fringe=0.4,
                                      redox_vadose_zone=0.,
                                      redox_shallow_aquifer=0.,
                                      redox_target_aquifer=0.,
                                      pH_vadose_zone=5,
                                      pH_shallow_aquifer=6,
                                      pH_target_aquifer=7,
                                      dissolved_organic_carbon_vadose_zone=10, 
                                      dissolved_organic_carbon_shallow_aquifer=4, 
                                      dissolved_organic_carbon_target_aquifer=2,
                                      fraction_organic_carbon_vadose_zone=0.001,
                                      fraction_organic_carbon_shallow_aquifer=0.0005,
                                      fraction_organic_carbon_target_aquifer=0.0005, 
                                      input_concentration = 100,)


# phreatic_dict = scheme1.make_dictionary()  
well1 = AnalyticalWell(scheme1)
well1.phreatic()     
df_flowline, df_particle = well1.export_to_df(what_to_export = 'omp_parameters')
df_particle

conc1 = Concentration(well1) #, df_particle, df_flowline)
conc1.compute_omp_removal()
conc1.df_particle.input_concentration #.omp_half_life_temperature_corrected
conc1.df_particle.steady_state_concentration_vadose_zone

# %%
well1.plot_travel_time_versus_radial_distance(xlim=[0, 4000])

well1.plot_travel_time_versus_cumulative_abstracted_water(xlim=[0, 100])
#%%
