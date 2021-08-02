#%% ----------------------------------------------------------------------------
# A. Hockin, March 2021
# KWR BO 402045-247
# ZZS verwijdering bodempassage
# AquaPriori - Transport Model
# With Martin Korevaar, Martin vd Schans, Steven Ros
#
# ------------------------------------------------------------------------------

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

#%%
# Test
# test_travel_time_distribution_phreatic()
# test_retardation_temp_koc_correction(substance='benzene', schematisation_type='phreatic')
# test_steady_concentration_temp_koc_correction_phreatic(substance='benzene')
# test_steady_concentration_temp_koc_correction_phreatic(substance='benzo(a)pyrene')
# test_steady_concentration_temp_koc_correction_phreatic(substance='AMPA')
# test_start_end_dates_contamination()
# test_compute_for_date_start_dates_contamination()
# test_compute_for_date_start_date_well()
# %%
phreatic_scheme = HydroChemicalSchematisation(schematisation_type='phreatic',
                                      computation_method= 'analytical', 
                                      what_to_export='omp',
                                      well_discharge=319.4*24, #m3/day
                                      vertical_anistropy_shallow_aquifer = (10/(35*500)),
                                      porosity_vadose_zone=0.38,
                                      porosity_shallow_aquifer=0.35,
                                      porosity_target_aquifer=0.35,
                                      recharge_rate=0.3/365.25, #m/day
                                      moisture_content_vadose_zone=0.15,
                                      ground_surface=22,
                                      thickness_vadose_zone_at_boundary=5,
                                      thickness_shallow_aquifer=10,
                                      thickness_target_aquifer=40,
                                      hor_permeability_target_aquifer=35,
                                      thickness_full_capillary_fringe=0.4,
                                      redox_vadose_zone='suboxic',
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
                                      temperature=11,
                                      solid_density_vadose_zone=2.650, 
                                      solid_density_shallow_aquifer=2.650, 
                                      solid_density_target_aquifer=2.650, 
                                      diameter_borehole=0.75,
                                      #diffuse parameters
                                      diffuse_input_concentration=100, #ug/L
                                      #point paramters
                                      concentration_point_contamination=100,
                                      distance_point_contamination_from_well=25, #5.45045, #
                                      depth_point_contamination=21, #m ASL
                                      discharge_point_contamination=1000,
                                      #dates
                                      start_date_well='1968-01-01', 
                                      start_date_contamination= '1966-01-01',
                                      compute_contamination_for_date='2050-01-01',
                                      end_date_contamination='1990-01-01',
                                      )

# phreatic_scheme.make_dictionary()  
phreatic_well = AnalyticalWell(phreatic_scheme)
    
phreatic_well.phreatic() 
#^ @MartinK in the functions run this? how to avoid errors running 'semiconfined' 
# here when the schematisation_type was defined as 'phreatic'?

phreatic_conc = SubstanceTransport(phreatic_well, substance = 'OMP-X')
# # phreatic_conc = SubstanceTransport(phreatic_well, substance = 'benzo(a)pyrene')
# phreatic_conc = SubstanceTransport(phreatic_well, substance = 'AMPA')

phreatic_conc.compute_omp_removal()

df_well_concentration = phreatic_conc.compute_concentration_in_well_at_date()

phreatic_conc.plot_concentration()

# phreatic_conc.plot_concentration(x_axis='Time') 
# %%


#%%

phreatic_scheme = HydroChemicalSchematisation(schematisation_type='phreatic',
                                    computation_method= 'analytical', 
                                    what_to_export='omp',
                                    well_discharge=319.4*24, #m3/day
                                    vertical_anistropy_shallow_aquifer = (10/(35*500)),
                                    porosity_vadose_zone=0.38,
                                    porosity_shallow_aquifer=0.35,
                                    porosity_target_aquifer=0.35,
                                    recharge_rate=0.3/365.25, #m/day
                                    moisture_content_vadose_zone=0.15,
                                    ground_surface=22,
                                    thickness_vadose_zone_at_boundary=5,
                                    thickness_shallow_aquifer=10,
                                    thickness_target_aquifer=40,
                                    hor_permeability_target_aquifer=35,
                                    thickness_full_capillary_fringe=0.4,
                                    redox_vadose_zone='suboxic',
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
                                    temperature=11,
                                    solid_density_vadose_zone=2.650, 
                                    solid_density_shallow_aquifer=2.650, 
                                    solid_density_target_aquifer=2.650, 
                                    diameter_borehole=0.75,
                                    #diffuse parameters
                                    diffuse_input_concentration=100, #ug/L
                                    #point paramters
                                    # concentration_point_contamination=100,
                                    # distance_point_contamination_from_well=25,
                                    # depth_point_contamination=21, #m ASL
                                    # discharge_point_contamination=1000,
                                    #dates
                                    start_date_well='1968-01-01', 
                                    start_date_contamination= '1966-01-01',
                                    compute_contamination_for_date='2050-01-01',
                                    end_date_contamination='1990-01-01',
                                    )
phreatic_well = AnalyticalWell(phreatic_scheme)
phreatic_well.phreatic() 
phreatic_conc = SubstanceTransport(phreatic_well, substance = 'OMP-X')
phreatic_conc.compute_omp_removal()
df_well_concentration = phreatic_conc.compute_concentration_in_well_at_date()

df_well_concentration.to_csv('phreatic_diffuse_only_test.csv')
df_well_concentration_test = read_csv(path / 'phreatic_diffuse_only_test.csv', index_col=0)
df_well_concentration_test = df_well_concentration_test.astype({'time': 'int32', 'date': 'datetime64[ns]', 'total_concentration_in_well': 'float64'}) 

assert_frame_equal(df_well_concentration, df_well_concentration_test, check_dtype=False)

#%%


phreatic_scheme = HydroChemicalSchematisation(schematisation_type='phreatic',
                                    computation_method= 'analytical', 
                                    what_to_export='omp',
                                    well_discharge=319.4*24, #m3/day
                                    vertical_anistropy_shallow_aquifer = (10/(35*500)),
                                    porosity_vadose_zone=0.38,
                                    porosity_shallow_aquifer=0.35,
                                    porosity_target_aquifer=0.35,
                                    recharge_rate=0.3/365.25, #m/day
                                    moisture_content_vadose_zone=0.15,
                                    ground_surface=22,
                                    thickness_vadose_zone_at_boundary=5,
                                    thickness_shallow_aquifer=10,
                                    thickness_target_aquifer=40,
                                    hor_permeability_target_aquifer=35,
                                    thickness_full_capillary_fringe=0.4,
                                    redox_vadose_zone='suboxic',
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
                                    temperature=11,
                                    solid_density_vadose_zone=2.650, 
                                    solid_density_shallow_aquifer=2.650, 
                                    solid_density_target_aquifer=2.650, 
                                    diameter_borehole=0.75,
                                    #diffuse parameters
                                    diffuse_input_concentration=0, #ug/L
                                    #point paramters
                                    concentration_point_contamination=100,
                                    distance_point_contamination_from_well=25,
                                    depth_point_contamination=21, #m ASL
                                    discharge_point_contamination=1000,
                                    #dates
                                    start_date_well='1968-01-01', 
                                    start_date_contamination= '1966-01-01',
                                    compute_contamination_for_date='2050-01-01',
                                    end_date_contamination='1990-01-01',
                                    )
phreatic_well = AnalyticalWell(phreatic_scheme)
phreatic_well.phreatic() 
phreatic_conc = SubstanceTransport(phreatic_well, substance = 'OMP-X')
phreatic_conc.compute_omp_removal()
df_well_concentration = phreatic_conc.compute_concentration_in_well_at_date()

df_well_concentration.to_csv('phreatic_point_only_test.csv')
df_well_concentration_test = read_csv(path / 'phreatic_point_only_test.csv', index_col=0)
df_well_concentration_test = df_well_concentration_test.astype({'time': 'int32', 'date': 'datetime64[ns]', 'total_concentration_in_well': 'float64'}) 

assert_frame_equal(df_well_concentration, df_well_concentration_test, check_dtype=False)
