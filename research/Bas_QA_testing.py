#%% ----------------------------------------------------------------------------
# A. Hockin, March 2021
# KWR BO 402045-247
# ZZS verwijdering bodempassage
# AquaPriori - Transport Model
# With Martin Korevaar, Martin vd Schans, Steven Ros
#
# ------------------------------------------------------------------------------

#### Notes ####
# Example for use in the read the docs and for Bas for 
# Quality assurance of the analytical code
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
from pandas import read_csv
from pandas import read_excel
import math
from scipy.special import kn as besselk

from pathlib import Path
try:
    from project_path import module_path #the dot says look in the current folder, this project_path.py file must be in the folder here
except ModuleNotFoundError:
    from project_path import module_path

from sutra2.Analytical_Well import *
from sutra2.Transport_Removal import *
from testing.test_transatomic import *
# get directory of this file
path = Path(__file__).parent #os.getcwd() #path of working directory

# path = os.getcwd()  # path of working directory

# %%
###########
# Example 1
###########
# Phreatic example using many of the default values
# Shows: 
# * how to make a basic schematisation for a phreatic aquifer
# * how to create an ANalyticalWell object
# * change the default OMP parameters
# * compute and plot the OMP concentration in the well 
#%% Step 1: Define the HydroChemicalSchematisation
# ==============================================
# The first step is to define the hydrogeochemistry of the system using the HydroChemicalSchematisation class.
# In this class you specify the:
    # * Computational method ('analytical' or 'modpath').
    # * The schematisation type ('phreatic', 'semiconfined').
    # .. ('riverbankfiltration', 'basinfiltration' coming soon).
    # * The removal function ('omp' )
    # .. *  coming soon, 'pathogen).
    # * Input the relevant parameters for the porous media, the hydrochemistry, hydrology and the contamination of interest


#Lets start with a simple example defining a HydroChemicalSchematisation object for a phreatic aquifer:
#%% SCHEME 1
phreatic_schematisation = HydroChemicalSchematisation(schematisation_type='phreatic',
                                      well_discharge=- 319.4*24, #7500, #m3/day
                                      recharge_rate=0.0008, #m/day
                                      thickness_vadose_zone_at_boundary=5,
                                      thickness_shallow_aquifer=10,
                                      thickness_target_aquifer=40,
                                      hor_permeability_target_aquifer=35,
                                      redox_vadose_zone='anoxic',
                                      redox_shallow_aquifer='anoxic',
                                      redox_target_aquifer='deeply_anoxic',
                                      pH_target_aquifer=7.,
                                      temperature=11.,
                                      substance='benzene',
                                      diffuse_input_concentration = 100, #ug/L
                                      )

phreatic_well = AnalyticalWell(phreatic_schematisation)

phreatic_well.phreatic() 

# radial_plot = phreatic_well.plot_travel_time_versus_radial_distance(xlim=[0, 2000], ylim=[1e3, 1e6])
# cumulative_plot = phreatic_well.plot_travel_time_versus_cumulative_abstracted_water(xlim=[0, 1], ylim=[1e3, 1e6])

# # Save the plots
# radial_plot.savefig('travel_time_versus_radial_distance_phreatic.png', dpi=300, bbox_inches='tight')
# cumulative_plot.savefig('travel_time_versus_cumulative_abs_water_phreatic.png', dpi=300, bbox_inches='tight')

phreatic_well.df_particle.head(10)
# phreatic_well.df_flowline.head(10)

#%% SCHEME 2

phreatic_schematisation = HydroChemicalSchematisation(schematisation_type='phreatic',
                                      well_discharge=-7500, #m3/day
                                      recharge_rate=0.0008, #m/day
                                      thickness_vadose_zone_at_boundary=1,
                                      thickness_shallow_aquifer=10,
                                      thickness_target_aquifer=20,
                                      hor_permeability_target_aquifer=35,
                                      redox_vadose_zone='anoxic',
                                      redox_shallow_aquifer='anoxic',
                                      redox_target_aquifer='deeply_anoxic',
                                      pH_target_aquifer=7.,
                                      temperature=11.,
                                      substance='benzene',
                                      diffuse_input_concentration = 100, #ug/L
                                      )

phreatic_well = AnalyticalWell(phreatic_schematisation)

phreatic_well.phreatic() 

# radial_plot = phreatic_well.plot_travel_time_versus_radial_distance(xlim=[0, 2000], ylim=[1e3, 1e6])
# cumulative_plot = phreatic_well.plot_travel_time_versus_cumulative_abstracted_water(xlim=[0, 1], ylim=[1e3, 1e6])

# # Save the plots
# radial_plot.savefig('travel_time_versus_radial_distance_phreatic.png', dpi=300, bbox_inches='tight')
# cumulative_plot.savefig('travel_time_versus_cumulative_abs_water_phreatic.png', dpi=300, bbox_inches='tight')

phreatic_well.df_particle.head(10)
phreatic_well.df_flowline.head(10)

#%% SCHEME 3

phreatic_schematisation = HydroChemicalSchematisation(schematisation_type='phreatic',
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
                                        thickness_vadose_zone_at_boundary=1,
                                        thickness_shallow_aquifer=1,
                                        thickness_target_aquifer=20,
                                        hor_permeability_target_aquifer=35,
                                        thickness_full_capillary_fringe=0.4,
                                        temperature=11,
                                        solid_density_vadose_zone= 2.650,
                                        solid_density_shallow_aquifer= 2.650,
                                        solid_density_target_aquifer= 2.650,
                                        diameter_borehole = 0.75,

                                      )
#%%

phreatic_well = AnalyticalWell(phreatic_schematisation)

phreatic_well.phreatic() 

radial_plot = phreatic_well.plot_travel_time_versus_radial_distance(xlim=[0, 2000], ylim=[1e3, 1e6])
cumulative_plot = phreatic_well.plot_travel_time_versus_cumulative_abstracted_water(xlim=[0, 1], ylim=[1e3, 1e6])

# # Save the plots
# radial_plot.savefig('travel_time_versus_radial_distance_phreatic.png', dpi=300, bbox_inches='tight')
# cumulative_plot.savefig('travel_time_versus_cumulative_abs_water_phreatic.png', dpi=300, bbox_inches='tight')

phreatic_well.df_particle.head(10)

df_particle =phreatic_well.df_particle

# df_particle.to_excel('bas_testing_QA_df_particle.xlsx')

##%% Plotting the situation modelled

crosssection_plot = phreatic_well.plot_depth_aquifers()
               
#%%