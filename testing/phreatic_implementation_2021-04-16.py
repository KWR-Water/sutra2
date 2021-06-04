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

# path = os.getcwd()  # path of working directory

#%%
# Test
test_travel_time_distribution_phreatic()
test_retardation_temp_koc_correction(substance='benzene', schematisation_type='phreatic')
test_steady_concentration_temp_koc_correction_phreatic(substance='benzene')
test_steady_concentration_temp_koc_correction_phreatic(substance='benzo(a)pyrene')
test_steady_concentration_temp_koc_correction_phreatic(substance='AMPA')
# %%
phreatic_scheme = HydroChemicalSchematisation(schematisation_type='phreatic',
                                    what_to_export='omp_parameters',
                                    # biodegradation_sorbed_phase = False,
                                      well_discharge=319.4*24,
                                    #   vertical_resistance_aquitard=500,
                                    vertical_anistropy_shallow_aquifer = (10/(35*500)),
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
                                      )

phreatic_scheme.make_dictionary()  
# phreatic_well = AnalyticalWell(phreatic_scheme)
# phreatic_well.phreatic()     

# phreatic_conc = Concentration(phreatic_well, substance = 'benzene')
# conc1 = Concentration(well1, substance = 'benzo(a)pyrene')
# conc1 = Concentration(well1, substance = 'AMPA')
# phreatic_conc.compute_omp_removal()
# phreatic_conc.df_particle #.steady_state_concentration

# plt.plot(conc1.df_flowline.flowline_id, conc1.df_flowline.total_breakthrough_travel_time/365.25)
# plt.xlabel('Flowline ID')
# plt.ylabel('Breakthrough time (years)')
# %%
phreatic_well.plot_travel_time_versus_radial_distance(xlim=[0, 4000], ylim=[1e3, 1e6])

phreatic_well.plot_travel_time_versus_cumulative_abstracted_water(xlim=[0, 100], ylim=[1e3, 1e6])
#%%
# Export dicts for steven

all_dicts = { 'simulation_parameters' : phreatic_scheme.simulation_parameters,
        'geo_parameters' : phreatic_scheme.geo_parameters,
        'ibound_parameters' : phreatic_scheme.ibound_parameters,
        'recharge_parameters' : phreatic_scheme.recharge_parameters,
        'well_parameters' : phreatic_scheme.well_parameters,
        'point_parameters' : phreatic_scheme.point_parameters,
        'substance_parameters' : phreatic_scheme.substance_parameters,
        'bas_parameters' : phreatic_scheme.bas_parameters,
}


f = open("phreatic_dict.txt","w")
f.write( str(all_dicts))
f.close()
#%%
# import the dictionary
file = open("phreatic_dict.txt", "r")
contents = file.read()
ph_dictionary = ast.literal_eval(contents)
file.close()

ph_dictionary

#%%