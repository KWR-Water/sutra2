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

semiconfined_scheme = HydroChemicalSchematisation(schematisation_type='semiconfined',
                                    computation_method = 'modpath',
                                    what_to_export='omp',
                                    # biodegradation_sorbed_phase = False,
                                      well_discharge=-319.4*24,
                                      # vertical_resistance_shallow_aquifer=500,
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
                                      # diameter_filterscreen = 0.2,
                                      point_input_concentration = 100.,
                                      discharge_point_contamination = 100.,#made up value
                                      top_clayseal = 17,
                                      compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),

                                      # substance = 'benzene',
                                      # halflife_suboxic=600,
                                      # partition_coefficient_water_organic_carbon = 3.3,
                                    )


semiconfined_well_dict = semiconfined_scheme.make_modflow_schematisation()
# Export dicts for steven

semi_dict_1  = { 'simulation_parameters' : semiconfined_scheme.simulation_parameters,
        'endpoint_id': semiconfined_scheme.endpoint_id,
        'mesh_refinement': semiconfined_scheme.mesh_refinement,
        'geo_parameters' : semiconfined_scheme.geo_parameters,
        'ibound_parameters' : semiconfined_scheme.ibound_parameters,
        'recharge_parameters' : semiconfined_scheme.recharge_parameters,
        'well_parameters' : semiconfined_scheme.well_parameters,
        'point_parameters' : semiconfined_scheme.point_parameters,
        'bas_parameters' : semiconfined_scheme.bas_parameters,
}


f = open("semiconfined_dict_nogravel.txt","w")
f.write( str(semi_dict_1))
f.close()
semi_dict_1

#%%
semiconfined_scheme = HydroChemicalSchematisation(schematisation_type='semiconfined',
                                    computation_method = 'modpath',
                                    what_to_export='omp',
                                    # biodegradation_sorbed_phase = False,
                                      well_discharge=-319.4*24,
                                      # vertical_resistance_shallow_aquifer=500,
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
                                      diameter_filterscreen = 0.2,
                                      point_input_concentration = 100.,
                                      discharge_point_contamination = 100.,#made up value
                                      top_clayseal = 17,
                                      compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),
                                      # substance = 'benzene',
                                      # halflife_suboxic=600,
                                      # partition_coefficient_water_organic_carbon = 3.3,
                                    )


semiconfined_well_dict = semiconfined_scheme.make_modflow_schematisation()
# Export dicts for steven

semi_dict_2 = { 'simulation_parameters' : semiconfined_scheme.simulation_parameters,
        'endpoint_id': semiconfined_scheme.endpoint_id,
        'mesh_refinement': semiconfined_scheme.mesh_refinement,
        'geo_parameters' : semiconfined_scheme.geo_parameters,
        'ibound_parameters' : semiconfined_scheme.ibound_parameters,
        'recharge_parameters' : semiconfined_scheme.recharge_parameters,
        'well_parameters' : semiconfined_scheme.well_parameters,
        'point_parameters' : semiconfined_scheme.point_parameters,
        'bas_parameters' : semiconfined_scheme.bas_parameters,
}

f = open("semiconfined_dict_gravelpack.txt","w")
f.write( str(semi_dict_2))
f.close()
semi_dict_2
#%%
# # import the dictionary
# file = open("semiconfined_dict.txt", "r")
# contents = file.read()
# dictionary = ast.literal_eval(contents)
# file.close()

# dictionary
# recharge_parameters = dictionary['recharge_parameters']
#%%

# %%
phreatic_scheme= HydroChemicalSchematisation(schematisation_type='phreatic',
                                    computation_method = 'modpath',
                                    what_to_export='omp',
                                    # biodegradation_sorbed_phase = False,
                                      well_discharge=-319.4*24,
                                      # vertical_resistance_shallow_aquifer=500,
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
                                      # diameter_filterscreen = 0.2,
                                      point_input_concentration = 100.,
                                      discharge_point_contamination = 100.,#made up value
                                      top_clayseal = 17,
                                      compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),

                                      # substance = 'benzene',
                                      # halflife_suboxic=600,
                                      # partition_coefficient_water_organic_carbon = 3.3,
                                    )

phreatic_scheme.make_modflow_schematisation()

phreatic_dict_1 = { 'simulation_parameters' : phreatic_scheme.simulation_parameters,
        'endpoint_id': phreatic_scheme.endpoint_id,
        'mesh_refinement': phreatic_scheme.mesh_refinement,
        'geo_parameters' : phreatic_scheme.geo_parameters,
        'ibound_parameters' : phreatic_scheme.ibound_parameters,
        'recharge_parameters' : phreatic_scheme.recharge_parameters,
        'well_parameters' : phreatic_scheme.well_parameters,
        'point_parameters' : phreatic_scheme.point_parameters,
        'bas_parameters' : phreatic_scheme.bas_parameters,
}
f = open("phreatic_dict_nogravel.txt","w")
f.write( str(phreatic_dict_1))
f.close()
phreatic_dict_1


#%%
phreatic_scheme= HydroChemicalSchematisation(schematisation_type='phreatic',
                                    computation_method = 'modpath',
                                    what_to_export='omp',
                                    # biodegradation_sorbed_phase = False,
                                      well_discharge=-319.4*24,
                                      # vertical_resistance_shallow_aquifer=500,
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
                                      diameter_filterscreen = 0.2,
                                      point_input_concentration = 100.,
                                      discharge_point_contamination = 100.,#made up value
                                      top_clayseal = 17,
                                      compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),
                                      # substance = 'benzene',
                                      # halflife_suboxic=600,
                                      # partition_coefficient_water_organic_carbon = 3.3,
                                    )

phreatic_scheme.make_modflow_schematisation()

phreatic_dict_2 = { 'simulation_parameters' : phreatic_scheme.simulation_parameters,
        'endpoint_id': phreatic_scheme.endpoint_id,
        'mesh_refinement': phreatic_scheme.mesh_refinement,
        'geo_parameters' : phreatic_scheme.geo_parameters,
        'ibound_parameters' : phreatic_scheme.ibound_parameters,
        'recharge_parameters' : phreatic_scheme.recharge_parameters,
        'well_parameters' : phreatic_scheme.well_parameters,
        'point_parameters' : phreatic_scheme.point_parameters,
        'bas_parameters' : phreatic_scheme.bas_parameters,
}

f = open("phreatic_dict_gravelpack.txt","w")
f.write( str(phreatic_dict_2))
f.close()
phreatic_dict_2

#%%
# PHREATIC TEST CASE SAME AS FOR THE UNIT TEST FOR PHREATIC TTD

phreatic_test_scheme = HydroChemicalSchematisation(schematisation_type='phreatic',
                                      computation_method= 'analytical',
                                      what_to_export='omp', # @alex: what_to_export sounds very cryptic and ad-hoc. maybe we can think of something better
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
                                      diameter_borehole = 0.75,
                                      )

phreatic_test_scheme.make_modflow_schematisation()

phreatic_test_scheme_dict= { 'simulation_parameters' : phreatic_scheme.simulation_parameters,
        'endpoint_id': semiconfined_scheme.endpoint_id,
        'mesh_refinement': semiconfined_scheme.mesh_refinement,
        'geo_parameters' : phreatic_scheme.geo_parameters,
        'ibound_parameters' : phreatic_scheme.ibound_parameters,
        'recharge_parameters' : phreatic_scheme.recharge_parameters,
        'well_parameters' : phreatic_scheme.well_parameters,
        'point_parameters' : phreatic_scheme.point_parameters,
        'bas_parameters' : phreatic_scheme.bas_parameters,
}

f = open("phreatic_test_scheme_dict.txt","w")
f.write( str(phreatic_test_scheme_dict))
f.close()
phreatic_test_scheme_dict

phreatic_test_well= AnalyticalWell(phreatic_test_scheme)
phreatic_test_well.phreatic()

df_particle =phreatic_test_well.df_particle
df_flowline = phreatic_test_well.df_flowline

df_particle.to_excel('phreatic_df_particle_Modflow.xlsx')
df_flowline.to_excel('phreatic_df_flowline_Modflow.xlsx')
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

semiconfined_test_scheme.make_modflow_schematisation()

semiconfined_test_scheme_dict= { 'simulation_parameters' : phreatic_scheme.simulation_parameters,
        'endpoint_id': semiconfined_scheme.endpoint_id,
        'mesh_refinement': semiconfined_scheme.mesh_refinement,
        'geo_parameters' : phreatic_scheme.geo_parameters,
        'ibound_parameters' : phreatic_scheme.ibound_parameters,
        'recharge_parameters' : phreatic_scheme.recharge_parameters,
        'well_parameters' : phreatic_scheme.well_parameters,
        'point_parameters' : phreatic_scheme.point_parameters,
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