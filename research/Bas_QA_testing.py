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

from greta.Analytical_Well import *
from greta.Substance_Transport import *
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
                                      well_discharge= 319.4*24, #7500, #m3/day
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

radial_plot = phreatic_well.plot_travel_time_versus_radial_distance(xlim=[0, 2000], ylim=[1e3, 1e6])
cumulative_plot = phreatic_well.plot_travel_time_versus_cumulative_abstracted_water(xlim=[0, 1], ylim=[1e3, 1e6])

# # Save the plots
# radial_plot.savefig('travel_time_versus_radial_distance_phreatic.png', dpi=300, bbox_inches='tight')
# cumulative_plot.savefig('travel_time_versus_cumulative_abs_water_phreatic.png', dpi=300, bbox_inches='tight')

phreatic_well.df_particle.head(10)
# phreatic_well.df_flowline.head(10)

#%% SCHEME 2

phreatic_schematisation = HydroChemicalSchematisation(schematisation_type='phreatic',
                                      well_discharge=7500, #m3/day
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

radial_plot = phreatic_well.plot_travel_time_versus_radial_distance(xlim=[0, 2000], ylim=[1e3, 1e6])
cumulative_plot = phreatic_well.plot_travel_time_versus_cumulative_abstracted_water(xlim=[0, 1], ylim=[1e3, 1e6])

# # Save the plots
# radial_plot.savefig('travel_time_versus_radial_distance_phreatic.png', dpi=300, bbox_inches='tight')
# cumulative_plot.savefig('travel_time_versus_cumulative_abs_water_phreatic.png', dpi=300, bbox_inches='tight')

phreatic_well.df_particle.head(10)
phreatic_well.df_flowline.head(10)

#%% SCHEME 3

phreatic_schematisation = HydroChemicalSchematisation(schematisation_type='phreatic',
                                        computation_method= 'analytical',
                                        what_to_export='omp', # @alex: what_to_export sounds very cryptic and ad-hoc. maybe we can think of something better
                                        well_discharge=319.4*24,
                                        # vertical_resistance_aquitard=500,
                                        hor_permeability_shallow_aquifer = 0.02,
                                        vertical_anisotropy_shallow_aquifer = (10/(0.02*500)),
                                        porosity_vadose_zone=0.38,
                                        porosity_shallow_aquifer=0.35,
                                        porosity_target_aquifer=0.35,
                                        recharge_rate=0.3/365.25,
                                        moisture_content_vadose_zone=0.15,
                                        ground_surface = 22,
                                        # thickness_vadose_zone_at_boundary=1,
                                        # thickness_shallow_aquifer=1,
                                        # thickness_target_aquifer=20,
                                        # hor_permeability_target_aquifer=35,
                                        # thickness_full_capillary_fringe=0.4,
                                        # temperature=11,
                                        # solid_density_vadose_zone= 2.650,
                                        # solid_density_shallow_aquifer= 2.650,
                                        # solid_density_target_aquifer= 2.650,
                                        # diameter_borehole = 0.75,

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
df_x = pd.DataFrame({"fraction_flux": phreatic_well.schematisation.fraction_flux,
        "total_travel_time": phreatic_well.total_travel_time})

# fig = plt.figure(figsize=[10, 5])
# plt.plot(df_x.cumulative_fraction_abstracted_water, df_x.total_travel_time)
# plt.yscale('log')
# plt.xlabel('Cumulative fraction of abstracted water')
# plt.ylabel('Total travel time (days)')

# df_x.sort_values(by=['total_travel_time'])
df_x = pd.DataFrame({"fraction_flux": phreatic_well.schematisation.fraction_flux,
        "total_travel_time": phreatic_well.total_travel_time})

dfx = df_x.sort_values(by=['total_travel_time'], ascending=True)

percent_diffs = np.diff(phreatic_well.schematisation.fraction_flux)
percent_diffs = np.insert(percent_diffs,0,0., axis=0)

dfx['percent_diffs'] = percent_diffs

dfx['cumulative_fraction_abstracted_water'] = dfx['percent_diffs'].cumsum()

fig = plt.figure(figsize=[10, 5])
plt.plot(dfx.cumulative_fraction_abstracted_water, dfx.total_travel_time)
plt.yscale('log')
plt.xlabel('Cumulative fraction of abstracted water')
plt.ylabel('Total travel time (days)')

# df_x
#%%
# Step 3: View the Substance class (Optional)
# ===========================================
# You can retrieve the default substance parameters used to calculate the removal in the
# SubstanceTransport class. The data are stored in a dictionary

test_substance = Substance(substance_name='benzene')
test_substance.substance_dict

#%% Step 4: Run the SubstanceTransport class
# ========================================
# To calculate the removal and the steady-state concentration in each zone, create a concentration
# object by running the SubstanceTransport class with the phreatic_well object and specifying
# the OMP (or pathogen) of interest.

# In this example we use benzene. First we create the object and view the substance properties:
phreatic_concentration = SubstanceTransport(phreatic_well, substance = 'benzene')
phreatic_concentration.substance_dict

# Optional: You may specify a different value for the substance parameters, for example
# a different half-life for the anoxic redox zone. This can be input in the HydroChemicalSchematisation
# and this will be used in the calculation for the removal for the OMP. The AnalyticalWell and 
# phreatic() functions must be rerun:

phreatic_schematisation = HydroChemicalSchematisation(schematisation_type='phreatic',
                                      well_discharge=7500, #m3/day
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
                                      partition_coefficient_water_organic_carbon=2,
                                      dissociation_constant=1,
                                      halflife_suboxic=12, 
                                      halflife_anoxic=420, 
                                      halflife_deeply_anoxic=6000,
                                      )
phreatic_well = AnalyticalWell(phreatic_schematisation)
phreatic_well.phreatic() 
phreatic_concentration = SubstanceTransport(phreatic_well, substance = 'benzene')

# If you have specified a values for the substance (e.g. half-life, pKa, log_Koc),
# the default value is overriden and used in the calculation of the removal. You can
# view the updated substance dictionary from the concentration object:

phreatic_concentration.substance_dict

# Then we compute the removal by running the 'compute_omp_removal' function:
# phreatic_concentration.compute_omp_removal()

phreatic_concentration.compute_omp_removal()

# Once the removal has been calculated, you can view the steady-state concentration
# and breakthrough time per zone for the OMP in the df_particle:

phreatic_concentration.df_particle[['flowline_id', 'zone', 'steady_state_concentration', 'breakthrough_travel_time']].head(4)

# View the steady-state concentration of the flowline or the steady-state
# contribution of the flowline to the concentration in the well
phreatic_concentration.df_flowline[['flowline_id', 'breakthrough_concentration', 'total_breakthrough_travel_time']].head(5)

beznene_plot = phreatic_concentration.plot_concentration(ylim=[0,10 ])
beznene_plot.savefig('benzene_plot.png', dpi=300, bbox_inches='tight')

#%%
phreatic_well = AnalyticalWell(phreatic_schematisation)
phreatic_well.phreatic() 
phreatic_concentration = SubstanceTransport(phreatic_well, substance = 'OMP-X')
phreatic_concentration.compute_omp_removal()
omp_x_plot = phreatic_concentration.plot_concentration(ylim=[0,100 ])
omp_x_plot.savefig('omp_x_plot.png', dpi=300, bbox_inches='tight')


phreatic_well = AnalyticalWell(phreatic_schematisation)
phreatic_well.phreatic() 
phreatic_concentration = SubstanceTransport(phreatic_well, substance = 'benzo(a)pyrene')
phreatic_concentration.compute_omp_removal()
benzo_plot = phreatic_concentration.plot_concentration(ylim=[0,1])
benzo_plot.savefig('benzo_plot.png', dpi=300, bbox_inches='tight')


phreatic_well = AnalyticalWell(phreatic_schematisation)
phreatic_well.phreatic() 
phreatic_concentration = SubstanceTransport(phreatic_well, substance = 'AMPA')
phreatic_concentration.compute_omp_removal()
ampa_plot = phreatic_concentration.plot_concentration( ylim=[0,1])
ampa_plot.savefig('ampa_plot.png', dpi=300, bbox_inches='tight')

#%%
###########
# Example 2
###########
# Phreatic exmaple with both a diffuse and point source
# Also specifies the start, end dates for the wells and contamination

phreatic_scheme = HydroChemicalSchematisation(schematisation_type='phreatic',
                                      computation_method= 'analytical', 
                                      what_to_export='omp',
                                      well_discharge=319.4*24, #m3/day
                                      vertical_anisotropy_shallow_aquifer = (10/(35*500)),
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
                                      start_date_well=dt.datetime.strptime('1968-01-01', '%Y-%m-%d'),
                                      start_date_contamination= dt.datetime.strptime('1966-01-01',"%Y-%m-%d"),
                                      compute_contamination_for_date=dt.datetime.strptime('2050-01-01',"%Y-%m-%d"),
                                      end_date_contamination=dt.datetime.strptime('1990-01-01',"%Y-%m-%d"),
                                      )

phreatic_well = AnalyticalWell(phreatic_scheme)
phreatic_well.phreatic() 
phreatic_conc = SubstanceTransport(phreatic_well, substance = 'OMP-X')
# phreatic_conc = SubstanceTransport(phreatic_well, substance = 'benzene')
# phreatic_conc = SubstanceTransport(phreatic_well, substance = 'benzo(a)pyrene')
# phreatic_conc = SubstanceTransport(phreatic_well, substance = 'AMPA')
phreatic_conc.compute_omp_removal()
conc_plot = phreatic_conc.plot_concentration(x_axis='Time') 
# conc_plot.savefig('phreatic_diffuse_and_point_source.png', dpi=300, bbox_inches='tight')  
#%%
###########
# Example 3
###########
# Phreatic exmaple with both only point source
# Also specifies the start, end dates for the wells and contamination
phreatic_scheme = HydroChemicalSchematisation(schematisation_type='phreatic',
                                      computation_method= 'analytical', 
                                      what_to_export='omp',
                                      well_discharge=319.4*24, #m3/day
                                      vertical_anisotropy_shallow_aquifer = (10/(35*500)),
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
                                      distance_point_contamination_from_well=25, #5.45045, #
                                      depth_point_contamination=21, #m ASL
                                      discharge_point_contamination=1000,
                                      #dates
                                      start_date_well=dt.datetime.strptime('1968-01-01', '%Y-%m-%d'),
                                      start_date_contamination= dt.datetime.strptime('1966-01-01',"%Y-%m-%d"),
                                      compute_contamination_for_date=dt.datetime.strptime('2050-01-01',"%Y-%m-%d"),
                                      end_date_contamination=dt.datetime.strptime('1990-01-01',"%Y-%m-%d"),
                                      )

phreatic_well = AnalyticalWell(phreatic_scheme)
phreatic_well.phreatic() 
phreatic_conc = SubstanceTransport(phreatic_well, substance = 'OMP-X')
# phreatic_conc = SubstanceTransport(phreatic_well, substance = 'benzene')
# phreatic_conc = SubstanceTransport(phreatic_well, substance = 'benzo(a)pyrene')
# phreatic_conc = SubstanceTransport(phreatic_well, substance = 'AMPA')
phreatic_conc.compute_omp_removal()
conc_plot = phreatic_conc.plot_concentration(x_axis='Time') 
# conc_plot.savefig('phreatic_point_source_only.png', dpi=300, bbox_inches='tight') 

#%%
###########
# Example 4
###########
# Semiconfined exmaple with both diffuse and point source
# Also specifies the start, end dates for the wells and contamination

semiconfined_scheme = HydroChemicalSchematisation(schematisation_type='semiconfined',
                                      computation_method= 'analytical', 
                                      what_to_export='omp',
                                      well_discharge=319.4*24, #m3/day
                                      vertical_anisotropy_shallow_aquifer = (10/(35*500)),
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
                                      distance_point_contamination_from_well=25, #5.45045, #
                                      depth_point_contamination=21, #m ASL
                                      discharge_point_contamination=1000,
                                      #dates
                                      start_date_well=dt.datetime.strptime('1968-01-01', '%Y-%m-%d'),
                                      start_date_contamination= dt.datetime.strptime('1966-01-01',"%Y-%m-%d"),
                                      compute_contamination_for_date=dt.datetime.strptime('2050-01-01',"%Y-%m-%d"),
                                      end_date_contamination=dt.datetime.strptime('1990-01-01',"%Y-%m-%d"),
                                      )

semiconfined_well = AnalyticalWell(semiconfined_scheme)
semiconfined_well.semiconfined()   
semiconfined_well.df_flowline
semi_plot1 = semiconfined_well.plot_travel_time_versus_radial_distance(xlim=[0, 4000], ylim=[1e3, 1e6])
semi_plot2 = semiconfined_well.plot_travel_time_versus_cumulative_abstracted_water(xlim=[0, 1], ylim=[1e3, 1e10])
semiconfined_conc = SubstanceTransport(semiconfined_well, substance = 'OMP-X')
semiconfined_conc.compute_omp_removal()

# Plot the concentration by the start dates or by time since start (year = 0)
date_plot = semiconfined_conc.plot_concentration()
time_plot = semiconfined_conc.plot_concentration(x_axis='Time') 