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
from numpy.lib.histograms import _search_sorted_inclusive
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
import copy


from pathlib import Path
try:
    from project_path import module_path #the dot says looik in the current folder, this project_path.py file must be in the folder here
except ModuleNotFoundError:
    from project_path import module_path

from greta.Analytical_Well import *
from greta.Substance_Transport import *
from testing.test_transatomic import *
# get directory of this file
path = Path(__file__).parent #os.getcwd() #path of working directory

#FIX
# FIX THE DIAMETERS
# FIX THE DICIOTNARY FOR the UBSTANCE, ADD SUBSTANCE NAME,
# FIX R_START not r = zero (default) for steven
 
# IS GRAVEL RMIN GRVEL = RMAX WELL SCREEN
# CHECK THE DIAMETERS CAREFULLY FOR THESE INTERLINKING FUNCTION OF HTE WELLS

# #%%
# # Test
# test_travel_time_distribution_semiconfined()
# test_retardation_temp_koc_correction(substance='benzene', schematisation_type='semiconfined')
# test_retardation_temp_koc_correction(substance='benzo(a)pyrene', schematisation_type='semiconfined')
# test_retardation_temp_koc_correction(substance='AMPA', schematisation_type='semiconfined')
# test_steady_concentration_temp_koc_correction_semiconfined(substance='benzene')
# test_steady_concentration_temp_koc_correction_semiconfined(substance='benzo(a)pyrene')
# test_steady_concentration_temp_koc_correction_semiconfined(substance='AMPA')

#%%

semiconfined_scheme = HydroChemicalSchematisation(schematisation_type='semiconfined',
                                    computation_method = 'modpath',
                                    what_to_export='omp_parameters',
                                    # biodegradation_sorbed_phase = False,
                                      well_discharge=319.4*24,
                                      # vertical_resistance_aquitard=500,
                                      porosity_vadose_zone=0.38,
                                      porosity_shallow_aquifer=0.35,
                                      porosity_target_aquifer=0.35,
                                      recharge_rate=0.3/365.25,
                                      moisture_content_vadose_zone=0.15,
                                      ground_surface = 22.0,
                                      thickness_vadose_zone_at_boundary=5.0,
                                      thickness_shallow_aquifer=10.0,
                                      thickness_target_aquifer=40.0,
                                      hor_permeability_target_aquifer=35.0,
                                      hor_permeability_shallow_aquifer = 0.02,
                                      vertical_anistropy_shallow_aquifer = (10/(0.02*500)),
                                      thickness_full_capillary_fringe=0.4,
                                      redox_vadose_zone='anoxic', #'suboxic',
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
                                      diffuse_input_concentration = 100, #ug/L
                                      temperature=11.,
                                      solid_density_vadose_zone= 2.650, 
                                      solid_density_shallow_aquifer= 2.650, 
                                      solid_density_target_aquifer= 2.650, 
                                      diameter_borehole = 0.75,
                                      substance = 'benzo(a)pyrene',
                                      halflife_oxic= 530,
                                      halflife_anoxic= 2120,
                                      halflife_deeply_anoxic= 2120,
                                      partition_coefficient_water_organic_carbon= 6.43,
                                      dissociation_constant= 99,
                                      # diameter_filterscreen = 0.2,
                                      concentration_point_contamination = 100.,
                                      discharge_point_contamination = 100.,#made up value
                                      top_clayseal = 17,
                                      compute_contamination_for_date='2020-01-01',
                                      # substance = 'benzene',
                                      # halflife_oxic=600,
                                      # partition_coefficient_water_organic_carbon = 3.3,
                                    )

semiconfined_well_dict = semiconfined_scheme.make_dictionary()  
# Export dicts for steven

semi_dict_1  = { 'simulation_parameters' : semiconfined_scheme.simulation_parameters,
        'geo_parameters' : semiconfined_scheme.geo_parameters,
        'ibound_parameters' : semiconfined_scheme.ibound_parameters,
        'recharge_parameters' : semiconfined_scheme.recharge_parameters,
        'well_parameters' : semiconfined_scheme.well_parameters,
        'point_parameters' : semiconfined_scheme.point_parameters,
        'substance_parameters' : semiconfined_scheme.substance_parameters,
        'bas_parameters' : semiconfined_scheme.bas_parameters,
}

semi_dict_1 
# f = open("semiconfined_dict_nogravel.txt","w")
# f.write( str(semi_dict_1))
# f.close()

#%%
semiconfined_scheme = HydroChemicalSchematisation(schematisation_type='semiconfined',
                                    computation_method = 'modpath',
                                    what_to_export='omp_parameters',
                                    # biodegradation_sorbed_phase = False,
                                      well_discharge=319.4*24,
                                      # vertical_resistance_aquitard=500,
                                      porosity_vadose_zone=0.38,
                                      porosity_shallow_aquifer=0.35,
                                      porosity_target_aquifer=0.35,
                                      recharge_rate=0.3/365.25,
                                      moisture_content_vadose_zone=0.15,
                                      ground_surface = 22.0,
                                      thickness_vadose_zone_at_boundary=5.0,
                                      thickness_shallow_aquifer=10.0,
                                      thickness_target_aquifer=40.0,
                                      hor_permeability_target_aquifer=35.0,
                                      hor_permeability_shallow_aquifer = 0.02,
                                      vertical_anistropy_shallow_aquifer = (10/(0.02*500)),
                                      thickness_full_capillary_fringe=0.4,
                                      redox_vadose_zone='anoxic', #'suboxic',
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
                                      input_concentration = 100.,
                                      temperature=11.,
                                      solid_density_vadose_zone= 2.650, 
                                      solid_density_shallow_aquifer= 2.650, 
                                      solid_density_target_aquifer= 2.650, 
                                      diameter_borehole = 0.75,
                                      substance = 'benzo(a)pyrene',
                                      halflife_oxic= 530,
                                      halflife_anoxic= 2120,
                                      halflife_deeply_anoxic= 2120,
                                      partition_coefficient_water_organic_carbon= 6.43,
                                      dissociation_constant= 99,
                                      diameter_filterscreen = 0.2,
                                      concentration_point_contamination = 100.,
                                      discharge_point_contamination = 100.,#made up value
                                      top_clayseal = 17,
                                      compute_contamination_for_date='2020-01-01',
                                      # substance = 'benzene',
                                      # halflife_oxic=600,
                                      # partition_coefficient_water_organic_carbon = 3.3,
                                    )


semiconfined_well_dict = semiconfined_scheme.make_dictionary()  
# Export dicts for steven

semi_dict_2 = { 'simulation_parameters' : semiconfined_scheme.simulation_parameters,
        'geo_parameters' : semiconfined_scheme.geo_parameters,
        'ibound_parameters' : semiconfined_scheme.ibound_parameters,
        'recharge_parameters' : semiconfined_scheme.recharge_parameters,
        'well_parameters' : semiconfined_scheme.well_parameters,
        'point_parameters' : semiconfined_scheme.point_parameters,
        'substance_parameters' : semiconfined_scheme.substance_parameters,
        'bas_parameters' : semiconfined_scheme.bas_parameters,
}

f = open("semiconfined_dict_gravelpack.txt","w")
f.write( str(semi_dict_2))
f.close()
#%%
# # import the dictionary
# file = open("semiconfined_dict.txt", "r")
# contents = file.read()
# dictionary = ast.literal_eval(contents)
# file.close()

# dictionary
# recharge_parameters = dictionary['recharge_parameters']
#%%

# Test
test_travel_time_distribution_phreatic()
test_retardation_temp_koc_correction(substance='benzene', schematisation_type='phreatic')
test_steady_concentration_temp_koc_correction_phreatic(substance='benzene')
test_steady_concentration_temp_koc_correction_phreatic(substance='benzo(a)pyrene')
test_steady_concentration_temp_koc_correction_phreatic(substance='AMPA')
# %%
phreatic_scheme= HydroChemicalSchematisation(schematisation_type='phreatic',
                                    computation_method = 'modpath',
                                    what_to_export='omp_parameters',
                                    # biodegradation_sorbed_phase = False,
                                      well_discharge=319.4*24,
                                      # vertical_resistance_aquitard=500,
                                      porosity_vadose_zone=0.38,
                                      porosity_shallow_aquifer=0.35,
                                      porosity_target_aquifer=0.35,
                                      recharge_rate=0.3/365.25,
                                      moisture_content_vadose_zone=0.15,
                                      ground_surface = 22.0,
                                      thickness_vadose_zone_at_boundary=5.0,
                                      thickness_shallow_aquifer=10.0,
                                      thickness_target_aquifer=40.0,
                                      hor_permeability_target_aquifer=35.0,
                                      hor_permeability_shallow_aquifer = 0.02,
                                      vertical_anistropy_shallow_aquifer = (10/(0.02*500)),
                                      thickness_full_capillary_fringe=0.4,
                                      redox_vadose_zone='anoxic', #'suboxic',
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
                                      substance = 'benzo(a)pyrene',
                                      halflife_oxic= 530,
                                      halflife_anoxic= 2120,
                                      halflife_deeply_anoxic= 2120,
                                      partition_coefficient_water_organic_carbon= 6.43,
                                      dissociation_constant= 99,
                                      # diameter_filterscreen = 0.2,
                                      concentration_point_contamination = 100.,
                                      discharge_point_contamination = 100.,#made up value
                                      top_clayseal = 17,
                                      compute_contamination_for_date='2020-01-01',

                                      # substance = 'benzene',
                                      # halflife_oxic=600,
                                      # partition_coefficient_water_organic_carbon = 3.3,
                                    )

phreatic_scheme.make_dictionary()  

phreatic_dict_1 = { 'simulation_parameters' : phreatic_scheme.simulation_parameters,
        'geo_parameters' : phreatic_scheme.geo_parameters,
        'ibound_parameters' : phreatic_scheme.ibound_parameters,
        'recharge_parameters' : phreatic_scheme.recharge_parameters,
        'well_parameters' : phreatic_scheme.well_parameters,
        'point_parameters' : phreatic_scheme.point_parameters,
        'substance_parameters' : phreatic_scheme.substance_parameters,
        'bas_parameters' : phreatic_scheme.bas_parameters,
}
phreatic_dict_1
# f = open("phreatic_dict_nogravel.txt","w")
# f.write( str(phreatic_dict_1))
# f.close()

#%%
phreatic_scheme= HydroChemicalSchematisation(schematisation_type='phreatic',
                                    computation_method = 'modpath',
                                    what_to_export='omp_parameters',
                                    # biodegradation_sorbed_phase = False,
                                      well_discharge=319.4*24,
                                      # vertical_resistance_aquitard=500,
                                      porosity_vadose_zone=0.38,
                                      porosity_shallow_aquifer=0.35,
                                      porosity_target_aquifer=0.35,
                                      recharge_rate=0.3/365.25,
                                      moisture_content_vadose_zone=0.15,
                                      ground_surface = 22.0,
                                      thickness_vadose_zone_at_boundary=5.0,
                                      thickness_shallow_aquifer=10.0,
                                      thickness_target_aquifer=40.0,
                                      hor_permeability_target_aquifer=35.0,
                                      hor_permeability_shallow_aquifer = 0.02,
                                      vertical_anistropy_shallow_aquifer = (10/(0.02*500)),
                                      thickness_full_capillary_fringe=0.4,
                                      redox_vadose_zone='anoxic', #'suboxic',
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
                                      input_concentration = 100.,
                                      temperature=11.,
                                      solid_density_vadose_zone= 2.650, 
                                      solid_density_shallow_aquifer= 2.650, 
                                      solid_density_target_aquifer= 2.650, 
                                      diameter_borehole = 0.75,
                                      substance = 'benzo(a)pyrene',
                                      halflife_oxic= 530,
                                      halflife_anoxic= 2120,
                                      halflife_deeply_anoxic= 2120,
                                      partition_coefficient_water_organic_carbon= 6.43,
                                      dissociation_constant= 99,
                                      diameter_filterscreen = 0.2,
                                      concentration_point_contamination = 100.,
                                      discharge_point_contamination = 100.,#made up value
                                      top_clayseal = 17,
                                      compute_contamination_for_date='2020-01-01',
                                      # substance = 'benzene',
                                      # halflife_oxic=600,
                                      # partition_coefficient_water_organic_carbon = 3.3,
                                    )

phreatic_scheme.make_dictionary()  

phreatic_dict_2 = { 'simulation_parameters' : phreatic_scheme.simulation_parameters,
        'geo_parameters' : phreatic_scheme.geo_parameters,
        'ibound_parameters' : phreatic_scheme.ibound_parameters,
        'recharge_parameters' : phreatic_scheme.recharge_parameters,
        'well_parameters' : phreatic_scheme.well_parameters,
        'point_parameters' : phreatic_scheme.point_parameters,
        'substance_parameters' : phreatic_scheme.substance_parameters,
        'bas_parameters' : phreatic_scheme.bas_parameters,
}

f = open("phreatic_dict_gravelpack.txt","w")
f.write( str(phreatic_dict_2))
f.close()

#%%