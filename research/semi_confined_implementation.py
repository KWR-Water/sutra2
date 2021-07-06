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
import ast

from pathlib import Path
try:
    from project_path import module_path #the dot says look in the current folder, this project_path.py file must be in the folder here
except ModuleNotFoundError:
    from project_path import module_path

from greta.Analytical_Well import *
from greta.Substance_Transport import *
# if change classes to seperate files, then import them seperately here AH

from testing.test_transatomic import *

# get directory of this file
path = Path(__file__).parent #os.getcwd() #path of working directory

path = os.getcwd() #path of working directory


#%%
# Test
test_travel_time_distribution_semiconfined()
test_retardation_temp_koc_correction(substance='benzene', schematisation_type='semiconfined')
test_retardation_temp_koc_correction(substance='benzo(a)pyrene', schematisation_type='semiconfined')
test_retardation_temp_koc_correction(substance='AMPA', schematisation_type='semiconfined')
test_steady_concentration_temp_koc_correction_semiconfined(substance='benzene')
test_steady_concentration_temp_koc_correction_semiconfined(substance='benzo(a)pyrene')
test_steady_concentration_temp_koc_correction_semiconfined(substance='AMPA')

#%%

semiconfined_scheme = HydroChemicalSchematisation(schematisation_type='semiconfined',
                                      what_to_export='omp_parameters',
                                      well_discharge=319.4*24,
                                      vertical_anistropy_shallow_aquifer = (10/(0.02*500)),
                                      porosity_vadose_zone=0.38,
                                      porosity_shallow_aquifer=0.35,
                                      porosity_target_aquifer=0.35,
                                      recharge_rate=0.3/365.25,
                                      moisture_content_vadose_zone=0.15,
                                      ground_surface = 22,
                                      thickness_vadose_zone_at_boundary=5,
                                      thickness_shallow_aquifer=10,
                                      thickness_target_aquifer=40,
                                      hor_permeability_target_aquifer=35,
                                      hor_permeability_shallow_aquifer = 0.02,
                                      thickness_full_capillary_fringe=0.4,
                                      redox_vadose_zone='anoxic', #'suboxic',
                                      redox_shallow_aquifer='anoxic',
                                      redox_target_aquifer='deeply_anoxic',
                                      pH_vadose_zone=5,
                                      pH_shallow_aquifer=6,
                                      pH_target_aquifer=7,
                                      dissolved_organic_carbon_vadose_zone=10, 
                                      dissolved_organic_carbon_shallow_aquifer=4, 
                                      dissolved_organic_carbon_target_aquifer=2,
                                      fraction_organic_carbon_vadose_zone=0.001,
                                      fraction_organic_carbon_shallow_aquifer=0.0005,
                                      fraction_organic_carbon_target_aquifer=0.0005, 
                                      # diffuse_input_concentration = 100, #ug/L
                                      concentration_point_contamination = 100,
                                      distance_point_contamination_from_well = 25, #5.45045, #
                                      depth_point_contamination =21, #m ASL
                                      discharge_point_contamination=1000,
                                      temperature=11,
                                      solid_density_vadose_zone= 2.650, 
                                      solid_density_shallow_aquifer= 2.650, 
                                      solid_density_target_aquifer= 2.650, 
                                      diameter_borehole = 0.75,
                                      # substance = 'benzo(a)pyrene',
                                      # substance = 'benzene',
                                      # halflife_oxic=600,
                                      # partition_coefficient_water_organic_carbon = 3.3,
                                    )


# semiconfined_well_dict = semiconfined_scheme.make_dictionary()  
semiconfined_well = AnalyticalWell(semiconfined_scheme) #.semiconfined()
semiconfined_well.semiconfined()   
semiconfined_conc = SubstanceTransport(semiconfined_well, substance = 'OMP-X')

# # semiconfined_conc = SubstanceTransport(semiconfined_well, substance = 'benzo(a)pyrene')
# semiconfined_conc = SubstanceTransport(semiconfined_well, substance = 'benzene')

semiconfined_conc.compute_omp_removal()
semiconfined_conc.df_flowline
semiconfined_conc.df_particle
semiconfined_conc.plot_concentration(xlim=[0, 100], ylim=[0,1 ])

# semiconfined_conc.compute_omp_removal()
# # semiconfined_conc.df_particle #.steady_state_concentration
# semiconfined_conc.df_particle
# semiconfined_conc.substance_dict

# %%
semiconfined_well.plot_travel_time_versus_radial_distance(xlim=[0, 4000], ylim=[1e3, 1e6])
semiconfined_well.plot_travel_time_versus_cumulative_abstracted_water(xlim=[0, 1], ylim=[1e3, 1e6])
#%%
# Export dicts for steven

all_dicts = { 'simulation_paramters' : semiconfined_scheme.simulation_paramters,
        'geo_parameters' : semiconfined_scheme.geo_parameters,
        'ibound_parameters' : semiconfined_scheme.ibound_parameters,
        'recharge_parameters' : semiconfined_scheme.recharge_parameters,
        'well_parameters' : semiconfined_scheme.well_parameters,
        'point_parameters' : semiconfined_scheme.point_parameters,
        'substance_parameters' : semiconfined_scheme.substance_parameters,
        'bas_parameters' : semiconfined_scheme.bas_parameters,
}


f = open("semiconfined_dict.txt","w")
f.write( str(all_dicts))
f.close()
#%%
# import the dictionary
file = open("semiconfined_dict.txt", "r")
contents = file.read()
dictionary = ast.literal_eval(contents)
file.close()

dictionary
# recharge_parameters = dictionary['recharge_parameters']
#%%