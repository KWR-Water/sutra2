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

from sutra2.Analytical_Well import *
from sutra2.Transport_Removal import *
from testing.test_transatomic import *
# get directory of this file
path = Path(__file__).parent #os.getcwd() #path of working directory


#%%

#DICTIONARIES START HERE
#%%
#%%
# PHREATIC TEST 

phreatic_scheme = HydroChemicalSchematisation(schematisation_type='phreatic',
                                      computation_method= 'analytical',
                                      what_to_export='omp',
                                      well_discharge=-319.4*24, #m3/day
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
                                      point_input_concentration=100,
                                      distance_point_contamination_from_well=25,
                                      depth_point_contamination=21, #m ASL
                                      discharge_point_contamination=-1000,
                                      #dates
                                      start_date_well=dt.datetime.strptime('1968-01-01',"%Y-%m-%d"),
                                      start_date_contamination= dt.datetime.strptime('1966-01-01',"%Y-%m-%d"),
                                      compute_contamination_for_date=dt.datetime.strptime('2050-01-01',"%Y-%m-%d"),
                                      end_date_contamination=dt.datetime.strptime('1990-01-01',"%Y-%m-%d"),

                                      )
phreatic_well = AnalyticalWell(phreatic_scheme)
phreatic_well.phreatic()
phreatic_conc = Transport(phreatic_well, substance = 'OMP-X')
phreatic_conc.compute_omp_removal()
phreatic_conc.compute_concentration_in_well_at_date()

# df_particle =phreatic_test_well.df_particle
df_flowline = phreatic_conc.df_flowline
df_flowline.head(10)
# df_particle.to_excel('phreatic_df_particle_Modflow.xlsx')
# df_flowline.to_excel('phreatic_df_flowline_Modflow.xlsx')
#%%
# SEMICONFINED TEST CASE SAME AS FOR THE UNIT TEST FOR SEMICONFINED TTD

semiconfined_test_scheme = HydroChemicalSchematisation(schematisation_type='semiconfined',
                                        computation_method= 'analytical',
                                        what_to_export='omp',
                                        well_discharge=-319.4*24,
                                        # vertical_resistance_shallow_aquifer=500,
                                        hor_permeability_shallow_aquifer = 0.02,
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
                                        # KD=1400,
                                        thickness_full_capillary_fringe=0.4,
                                        temperature=11,
                                        solid_density_vadose_zone= 2.650,
                                        solid_density_shallow_aquifer= 2.650,
                                        solid_density_target_aquifer= 2.650,
                                        diameter_borehole = 0.75,)

semiconfined_test_scheme.make_dictionary()  

semiconfined_test_scheme_dict= { 'simulation_parameters' : phreatic_scheme.simulation_parameters,
        'endpoint_id': semiconfined_test_scheme.endpoint_id,
        'mesh_refinement': semiconfined_test_scheme.mesh_refinement,
        'geo_parameters' : phreatic_scheme.geo_parameters,
        'ibound_parameters' : phreatic_scheme.ibound_parameters,
        'recharge_parameters' : phreatic_scheme.recharge_parameters,
        'well_parameters' : phreatic_scheme.well_parameters,
        'point_parameters' : phreatic_scheme.point_parameters,
        'substance_parameters' : phreatic_scheme.substance_parameters,
        'bas_parameters' : phreatic_scheme.bas_parameters,
}

f = open("semiconfined_test_scheme_dict.txt","w")
f.write( str(semiconfined_test_scheme_dict))
f.close()
semiconfined_test_scheme_dict

semiconfined_test_well= AnalyticalWell(semiconfined_test_scheme)

semiconfined_test_well.semiconfined()

df_particle =semiconfined_test_well.df_particle
df_flowline = semiconfined_test_well.df_flowline

df_particle.to_excel('semiconfined_df_particle_Modflow.xlsx')
df_flowline.to_excel('semiconfined_df_flowline_Modflow.xlsx')

#%%