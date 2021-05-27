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

from pathlib import Path
try:
    from project_path import module_path #the dot says looik in the current folder, this project_path.py file must be in the folder here
except ModuleNotFoundError:
    from project_path import module_path

from greta.draft_transport_function import *

from testing.test_transatomic import *
# get directory of this file
path = Path(__file__).parent #os.getcwd() #path of working directory

#%%
# Test
test_travel_time_distribution_semiconfined()
test_retardation_temp_koc_correction(substance='benzene', schematisation_type='semiconfined')
# test_retardation_temp_koc_correction(substance='benzo(a)pyrene', schematisation_type='semiconfined')
# test_retardation_temp_koc_correction(substance='AMPA', schematisation_type='semiconfined')
test_steady_concentration_temp_koc_correction_semiconfined(substance='benzene')
# test_steady_concentration_temp_koc_correction_semiconfined(substance='benzo(a)pyrene')
# test_steady_concentration_temp_koc_correction_semiconfined(substance='AMPA')

#%%

semiconfined_scheme = HydroChemicalSchematisation(schematisation_type='semiconfined',
                                    what_to_export='omp_parameters',
                                    # biodegradation_sorbed_phase = False,
                                      well_discharge=319.4*24,
                                      # vertical_resistance_aquitard=500,
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
                                      vertical_anistropy_shallow_aquifer = (10/(0.02*500)),
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
                                      input_concentration = 100,
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


# phreatic_dict = scheme1.make_dictionary()  
semiconfined_well = AnalyticalWell(semiconfined_scheme) #.semiconfined()
semiconfined_well.semiconfined()   
# semiconfined_conc = Concentration(semiconfined_well, substance = 'benzo(a)pyrene')
semiconfined_conc = Concentration(semiconfined_well, substance = 'benzene')


semiconfined_conc.compute_omp_removal()
# semiconfined_conc.df_particle #.steady_state_concentration
semiconfined_conc.df_particle
semiconfined_conc.substance_dict

# %%
semiconfined_well.plot_travel_time_versus_radial_distance(xlim=[0, 4000], ylim=[1e3, 1e6])
semiconfined_well.plot_travel_time_versus_cumulative_abstracted_water(xlim=[0, 100], ylim=[1e3, 1e6])
#%%
substance_parameters = {
    'benzene': {
        'log_Koc': 7,
        'pKa': None,
        'omp_half_life': {
            'suboxic': 450,
            'anoxic': 620,
            'deeply_anoxic': None
            },
        }
    }
# Substance dict here as placeholder for the actual database
substances_dict = { 
    'benzene': {
        'log_Koc': 1.92,
        'molar_mass': 78.1, 
        'pKa': 99,
        'omp_half_life': {
            'suboxic': 10.5,
            'anoxic': 420,
            'deeply_anoxic': 1e99,
            },
        },
    'AMPA': {
        'log_Koc': -0.36,
        'molar_mass': 111.04 , 
        'pKa': 0.4,
        'omp_half_life': {
            'suboxic': 46,
            'anoxic': 46,
            'deeply_anoxic': 1e99,
            },
        },
    'benzo(a)pyrene': {
        'log_Koc': 6.43,
        'molar_mass': 252.3, 
        'pKa': 99,
        'omp_half_life': {
            'suboxic': 530,
            'anoxic': 2120,
            'deeply_anoxic': 2120,
            },
        },
    }

#%%
substance = 'benzene'
old = substance_parameters.copy()
new = substances_dict[substance].copy() #['benzene'].copy()

# old.update( (k,v) for k,v in new.items() if v is None)

old

# %%
for key1, value1 in old.items():
  if old[key1] == new[key1]:
    for key, value in old.items():
      if type(value) is dict:
        for tkey, cvalue in value.items():
          if cvalue is None:
              old[key][tkey]= new[key][tkey]
      else:
          if value is None:
            old[key] = new[key]

old
# %%
