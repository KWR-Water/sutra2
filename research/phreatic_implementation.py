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
    from project_path import module_path #the dot says look in the current folder, this project_path.py file must be in the folder here
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
                                      well_discharge=319.4*24, #m3/day
                                      vertical_anistropy_shallow_aquifer = (10/(35*500)),
                                      porosity_vadose_zone=0.38,
                                      porosity_shallow_aquifer=0.35,
                                      porosity_target_aquifer=0.35,
                                      recharge_rate=0.3/365.25, #m/day
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
                                      #  diffuse_input_concentration = 100, #ug/L
                                      temperature=11,
                                      solid_density_vadose_zone= 2.650, 
                                      solid_density_shallow_aquifer= 2.650, 
                                      solid_density_target_aquifer= 2.650, 
                                      diameter_borehole = 0.75,

                                      concentration_point_contamination = 600,
                                      distance_point_contamination_from_well = 5.45045, #300,
                                      depth_point_contamination = 22, #m ASL
                                      discharge_point_contamination=10,
                                      )

# phreatic_scheme.make_dictionary()  
phreatic_well = AnalyticalWell(phreatic_scheme)
    
phreatic_well.phreatic() #@MartinK in the functions run this? how to avoid errors running 'semiconfined' here when the schematisation_type was defined as 'phreatic'?
# phreatic_well.df_particle
phreatic_conc = Concentration(phreatic_well, substance = 'benzene')
# # phreatic_conc = Concentration(phreatic_well, substance = 'benzo(a)pyrene')
# phreatic_conc = Concentration(phreatic_well, substance = 'AMPA')


phreatic_conc.compute_omp_removal()
phreatic_conc.df_flowline
# phreatic_conc.df_particle['steady_state_concentration']
# phreatic_conc.df_particle #.steady_state_concentration

# plt.plot(phreatic_conc.df_flowline.flowline_id, phreatic_conc.df_flowline.total_breakthrough_travel_time/365.25)
# plt.xlabel('Flowline ID')
# plt.ylabel('Breakthrough time (years)')

# phreatic_conc.df_particle
# %%
phreatic_well.plot_travel_time_versus_radial_distance(xlim=[0, 4000], ylim=[1e3, 1e6])

# phreatic_well.plot_travel_time_versus_cumulative_abstracted_water(xlim=[0, 100], ylim=[1e3, 1e6])

#%%

# %%
def calculate_unsaturated_zone_travel_time (phreatic_scheme, distance):
  head = (phreatic_scheme.groundwater_level - phreatic_scheme.well_discharge
          / (2 * math.pi * phreatic_scheme.KD) 
          * np.log(phreatic_scheme.radial_distance_recharge / distance )) #phreatic_scheme.radial_distance))

  thickness_vadose_zone_drawdown = (phreatic_scheme.groundwater_level 
                                    + phreatic_scheme.thickness_vadose_zone_at_boundary) - head

  travel_time_unsaturated = ((((phreatic_scheme.groundwater_level + thickness_vadose_zone_drawdown)
                              - phreatic_scheme.groundwater_level
                              - phreatic_scheme.thickness_full_capillary_fringe)
                              * phreatic_scheme.moisture_content_vadose_zone
                              + phreatic_scheme.thickness_full_capillary_fringe
                              * phreatic_scheme.porosity_vadose_zone)
                          / phreatic_scheme.recharge_rate)
  return travel_time_unsaturated, head

# Calculate travel times in aquifer for WATER
travel_time_unsaturated, head = calculate_unsaturated_zone_travel_time (phreatic_scheme = phreatic_scheme, 
                        distance = phreatic_scheme.distance_point_contamination_from_well)

travel_time_shallow_aquifer = ((phreatic_scheme.thickness_shallow_aquifer - (phreatic_scheme.groundwater_level - head))
                    * phreatic_scheme.porosity_shallow_aquifer / phreatic_scheme.recharge_rate)

# travel_time_target_aquifer = (phreatic_scheme.porosity_target_aquifer * phreatic_scheme.thickness_target_aquifer / phreatic_scheme.recharge_rate
#                     * np.log(1 / (1 -phreatic_scheme.fraction_flux)))

horizontal_travel_time_target_aquifer =((phreatic_scheme.porosity_target_aquifer * phreatic_scheme.thickness_target_aquifer) / phreatic_scheme.recharge_rate 
                            * math.log(phreatic_scheme.well_discharge / (phreatic_scheme.well_discharge - math.pi * phreatic_scheme.recharge_rate * phreatic_scheme.distance_point_contamination_from_well ** 2 ) )) 

# Calculate travel time for OMP (add RETARDATION)
travel_time_unsaturated_omp = travel_time_unsaturated * phreatic_conc.df_particle['retardation'].iloc[1]
travel_time_shallow_aquifer_omp = travel_time_shallow_aquifer * phreatic_conc.df_particle['retardation'].iloc[2]
travel_time_target_aquifer_omp = horizontal_travel_time_target_aquifer * phreatic_conc.df_particle['retardation'].iloc[3]

# Total travel time for OMP [days]
total_travel_time_omp = (travel_time_unsaturated_omp + travel_time_shallow_aquifer_omp
                    + travel_time_target_aquifer_omp)

#Calculate the concentrations
#after unsaturated zone... assume 100% for now
# Eq. 4.11 in report
fraction_unsat_in = 1
fraction_concentration_end_unsaturated_zone = fraction_unsat_in * 2 **(-1 * travel_time_unsaturated * phreatic_conc.df_particle['retardation'].iloc[1] / phreatic_conc.df_particle['omp_half_life'].iloc[1]  )
fraction_concentration_end_shallow = fraction_concentration_end_unsaturated_zone * 2 **(-1 * travel_time_shallow_aquifer * phreatic_conc.df_particle['retardation'].iloc[2] / phreatic_conc.df_particle['omp_half_life'].iloc[2]  )

fraction_concentration_end_target = fraction_concentration_end_shallow* 2 **(-1 * horizontal_travel_time_target_aquifer* phreatic_conc.df_particle['retardation'].iloc[3] / phreatic_conc.df_particle['omp_half_life'].iloc[3]  )

# %%
#%%

cumulative_fraction_abstracted_water =phreatic_well.cumulative_fraction_abstracted_water

# distance = phreatic_scheme.distance_point_contamination_from_well
distance = phreatic_well.radial_distance

column_names = ["total_travel_time", "travel_time_unsaturated", 
                "travel_time_shallow_aquifer", "travel_time_target_aquifer",
                "radial_distance", "head", "cumulative_fraction_abstracted_water", 
                "flowline_discharge"
                ]

if isinstance(cumulative_fraction_abstracted_water, float):
# if len(cumulative_fraction_abstracted_water) ==1 :
  phreatic_conc.flowline_discharge = cumulative_fraction_abstracted_water * phreatic_scheme.well_discharge
  # print(1)
else:
  phreatic_conc.flowline_discharge = (np.diff(np.insert(cumulative_fraction_abstracted_water,0,0., axis=0)))*phreatic_scheme.well_discharge
  print(2)
  distance = phreatic_well.radial_distance

# AH check the well_discharge calculations for the streamline... Pr*Q? 
data = [phreatic_conc.total_travel_time, phreatic_conc.travel_time_unsaturated, 
            phreatic_conc.travel_time_shallow_aquifer, phreatic_conc.travel_time_target_aquifer,
            distance, phreatic_conc.head, phreatic_conc.cumulative_fraction_abstracted_water,  
            phreatic_conc.flowline_discharge ]

# phreatic_conc.df_output = pd.DataFrame (data = np.transpose(data), columns=column_names)
phreatic_conc.df_output = pd.DataFrame ([data], columns=column_names)
phreatic_conc.df_output
# %%
