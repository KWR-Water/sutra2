
#%%
import pytest
from pandas import RangeIndex, read_csv
import datetime as dt
import numpy as np
import pandas as pd
import os
import ast  # abstract syntax trees
import sys
import copy
# path = os.getcwd()  # path of working directory
from pathlib import Path

from zmq import zmq_version_info

# Plotting modules
import matplotlib.pyplot as plt
import matplotlib 
import matplotlib.colors as colors

# sutra2 modules    
import sutra2.Analytical_Well as AW
import sutra2.ModPath_Well as mpw
import sutra2.Transport_Removal as TR

from pandas._testing import assert_frame_equal

# get directory of this file
path = Path(__file__).parent #os.getcwd() #path of working directory

# Create test files dir
testfiles_dir = os.path.join(path,"test_files")
if not os.path.exists(testfiles_dir):
    os.makedirs(testfiles_dir)

# modflow executable location
mf_exe = os.path.join(path, "mf2005.exe")
# modpath executable location
mp_exe = os.path.join(path, "mpath7.exe")

# <<<<<<< HEAD
#%%
# def test_modpath_run_phreatic_nogravelpack(organism_name = "MS2"):
#     ''' Phreatic scheme without gravelpack: modpath run.'''

#     #@Steven: mag weg uit testing: of final concentration toevoegen (geen echte 'assert'). Deze test (of variant ervan) naar readthedocs --> toon de dataframes
#     test_phrea = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
#                                     computation_method = 'modpath',
#                                     what_to_export='all',
#                                     removal_function = 'mbo',
#                                     # biodegradation_sorbed_phase = False,
#                                       well_discharge=-319.4*24,
#                                       # vertical_resistance_shallow_aquifer=500,
#                                       porosity_vadose_zone=0.38,
#                                       porosity_shallow_aquifer=0.35,
#                                       porosity_target_aquifer=0.35,
#                                       recharge_rate=0.3/365.25,
#                                       moisture_content_vadose_zone=0.15,
#                                       ground_surface = 22.0,
#                                       thickness_vadose_zone_at_boundary=5.0,
#                                       thickness_shallow_aquifer=10.0,
#                                       thickness_target_aquifer=40.0,
#                                       hor_permeability_target_aquifer=35.0,
#                                       hor_permeability_shallow_aquifer = 0.02,
#                                       thickness_full_capillary_fringe=0.4,
#                                       redox_vadose_zone='anoxic', #'suboxic',
#                                       redox_shallow_aquifer='anoxic',
#                                       redox_target_aquifer='deeply_anoxic',
#                                       pH_vadose_zone=5.,
#                                       pH_shallow_aquifer=6.,
#                                       pH_target_aquifer=7.,
#                                       dissolved_organic_carbon_vadose_zone=10., 
#                                       dissolved_organic_carbon_shallow_aquifer=4., 
#                                       dissolved_organic_carbon_target_aquifer=2.,
#                                       fraction_organic_carbon_vadose_zone=0.001,
#                                       fraction_organic_carbon_shallow_aquifer=0.0005,
#                                       fraction_organic_carbon_target_aquifer=0.0005, 
#                                       temp_water=11.,
#                                       solid_density_vadose_zone= 2.650, 
#                                       solid_density_shallow_aquifer= 2.650, 
#                                       solid_density_target_aquifer= 2.650, 
#                                       diameter_borehole = 0.75,
#                                       name = organism_name,
#                                       # diameter_filterscreen = 0.2,
#                                       point_input_concentration = 100.,
#                                       discharge_point_contamination = 100.,#made up value
#                                       top_clayseal = 17,
#                                       compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),
#                                       ncols_near_well = 20,
#                                       ncols_far_well = 80,
#                                     )

#     test_phrea.make_dictionary()
#     # Remove/empty point_parameters
#     test_phrea.point_parameters = {}

#     # print(test_phrea.__dict__)
#     modpath_phrea = mpw.ModPathWell(test_phrea,
#                             workspace = os.path.join(path,"test2_phrea_nogp"),
#                             modelname = "phreatic",
#                             bound_left = "xmin",
#                             bound_right = "xmax")
#     # print(modpath_phrea.__dict__)
#     # Run phreatic schematisation
#     modpath_phrea.run_model(run_mfmodel = True,
#                         run_mpmodel = True)

    # # Pollutant to define removal parameters for
    # organism = TR.MicrobialOrganism(organism_name = organism_name)

    # # Calculate advective microbial removal
    # modpath_removal = TR.Transport(modpath_phrea,
    #                                         pollutant = organism)
 
#     # Calculate advective microbial removal
#     # Final concentration per endpoint_id
#     C_final = {}
#     for endpoint_id in modpath_phrea.schematisation_dict.get("endpoint_id"):
#         df_particle, df_flowline, C_final[endpoint_id] = modpath_removal.calc_advective_microbial_removal(
#                                             modpath_phrea.df_particle, modpath_phrea.df_flowline, 
#                                             endpoint_id = endpoint_id,
#                                             trackingdirection = modpath_phrea.trackingdirection,
#                                             mu1 = 0.1151, grainsize = 0.00025, alpha0 = 0.037e-2, pH0 = 7.5,
#                                             temp_water = 11., rho_water = 999.703, organism_diam = 2.731e-6,
#                                             conc_start = 1., conc_gw = 0.)

#     # df_particle file name 
#     particle_fname = os.path.join(modpath_phrea.dstroot,modpath_phrea.schematisation_type + "_df_particle_microbial_removal.csv")
#     # Save df_particle 
#     df_particle.to_csv(particle_fname)
    
#     # df_flowline file name
#     flowline_fname = os.path.join(modpath_phrea.dstroot,modpath_phrea.schematisation_type + "_df_flowline_microbial_removal.csv")
#     # Save df_flowline
#     df_flowline.to_csv(flowline_fname)
    
#     assert modpath_phrea.success_mp

def test_modpath_run_horizontal_flow_points(organism_name = "MS2"):
    ''' Horizontal flow test in target_aquifer: modpath run.'''

    # @Steven: deze test hier laten staan (evt. nieuwe naamgeving workspace)
    # Wel mogen dataframes hier weg?

    # microbial removal properties
    # organism_name = 'MS2'
    alpha0 = {"suboxic": 1.e-3, "anoxic": 1.e-5, "deeply_anoxic": 1.e-5}
    pH0 = {"suboxic": 6.6, "anoxic": 6.8, "deeply_anoxic": 6.8}
    organism_diam =  2.33e-8
    mu1 = {"suboxic": 0.039,"anoxic": 0.023,"deeply_anoxic": 0.023}

    removal_parameters = {organism_name: 
                    {"organism_name": organism_name,
                        "alpha0": alpha0,
                        "pH0": pH0,
                        "organism_diam": organism_diam,
                        "mu1": mu1
                    }
                }
    # Removal parameters organism
    rem_parms = removal_parameters[organism_name]

    # test dataframe of final concentrations to assert the test
    summary_conc_fname_test = os.path.join(testfiles_dir,"Final_concentrations_horizontal_flow_test.csv")
    df_conc_horflow_test = pd.read_csv(summary_conc_fname_test,
                                        index_col = 0)
    df_conc_horflow_test.index.name = "Boundary_distance"

    # dataframes
    df_particle, df_flowline = {}, {}

    # Distances to boundary
    dist_boundary = list(range(10,100,5)) + list(range(100,560,10)) # [260] #  + 
    # [10,20, 50, 100, 150, 200,250,300,350,400,450,500,550]

    df_conc = pd.DataFrame(index=dist_boundary, columns = ["Final_concentration"])
    df_conc.index.name = "Boundary_distance"

    for iDist in dist_boundary:

        
        # well discharge
        well_discharge = -1000.
        # distance to boundary
        distance_boundary = float(iDist)
        # Thickness target aquifer
        thickness_target_aquifer = 20.
        # Center depth of target aquifer
        z_point = 0 - 0.1 - 10.
        x_point = distance_boundary - 0.25
        # Phreatic scheme without gravelpack: modpath run.
        test_conf_hor = AW.HydroChemicalSchematisation(schematisation_type='semiconfined',
                                    computation_method = 'modpath',
                                    what_to_export='all',
                                    removal_function = 'mbo',
                                    # biodegradation_sorbed_phase = False,
                                    well_discharge= well_discharge,
                                    # vertical_resistance_shallow_aquifer=500,
                                    porosity_vadose_zone=0.33,
                                    porosity_shallow_aquifer=0.33,
                                    porosity_target_aquifer=0.33,
                                    recharge_rate=0.,
                                    moisture_content_vadose_zone=0.15,
                                    ground_surface = 0.0,
                                    thickness_vadose_zone_at_boundary=0.,
                                    thickness_shallow_aquifer=0.1,
                                    thickness_target_aquifer= thickness_target_aquifer,
                                    hor_permeability_target_aquifer=10.0,
                                    hor_permeability_shallow_aquifer = 1.,
                                    vertical_anisotropy_shallow_aquifer = 10.,
                                    thickness_full_capillary_fringe=0.4,
                                    redox_vadose_zone='anoxic', #'suboxic',
                                    redox_shallow_aquifer='anoxic',
                                    redox_target_aquifer='anoxic',
                                    pH_vadose_zone=7.5,
                                    pH_shallow_aquifer=7.5,
                                    pH_target_aquifer=7.5,
                                    dissolved_organic_carbon_vadose_zone=1., 
                                    dissolved_organic_carbon_shallow_aquifer=1., 
                                    dissolved_organic_carbon_target_aquifer=1.,
                                    fraction_organic_carbon_vadose_zone=0.001,
                                    fraction_organic_carbon_shallow_aquifer=0.001,
                                    fraction_organic_carbon_target_aquifer=0.001, 
                                    temp_water=12.,
                                    solid_density_vadose_zone= 2.650, 
                                    solid_density_shallow_aquifer= 2.650, 
                                    solid_density_target_aquifer= 2.650, 
                                    diameter_borehole = 0.2,
                                    name = organism_name,
                                    diameter_filterscreen = 0.2,
                                    top_clayseal = 0,
                                    compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),
                                    model_radius = distance_boundary,
                                    ncols_near_well = 20,
                                    ncols_far_well = int((distance_boundary-thickness_target_aquifer)/0.5),
                                    # Point contamination
                                    point_input_concentration = 1.,
                                    discharge_point_contamination = abs(well_discharge),
                                    distance_point_contamination_from_well = x_point,
                                    depth_point_contamination = z_point,
                                    )

        test_conf_hor.make_dictionary()
        

        # Remove "vadose" layer from geo_parameters
        test_conf_hor.geo_parameters.pop("vadose_zone")

        ### Adjust ibound_parameters to add horizontal flow ###
        # Confined top boundary ; no recharge_parameters
        test_conf_hor.ibound_parameters.pop("top_boundary1")
        test_conf_hor.ibound_parameters.pop("inner_boundary_shallow_aquifer")

        # Add outer boundary for horizontal flow test
        test_conf_hor.ibound_parameters["outer_boundary_target_aquifer"] = {
                                'top': test_conf_hor.bottom_shallow_aquifer,
                                'bot': test_conf_hor.bottom_target_aquifer,
                                'xmin': test_conf_hor.model_radius - 0.5,
                                'xmax': test_conf_hor.model_radius,
                                'ibound': -1
                                }
        # Remove/empty concentration_boundary_parameters
        test_conf_hor.concentration_boundary_parameters = {}

        # ModPath well
        modpath_hor = mpw.ModPathWell(test_conf_hor, #test_phrea,
                                workspace = os.path.join(path,"test3_conf_hor_pnts"),
                                modelname = "confined_hor",
                                bound_left = "xmin",
                                bound_right = "xmax")
        # print(modpath_phrea.__dict__)
        obj_dict = modpath_hor.__dict__

        # Run phreatic schematisation
        modpath_hor.run_model(run_mfmodel = True,
                            run_mpmodel = True)

        # Pollutant to define removal parameters for
        organism = TR.MicrobialOrganism(organism_name = organism_name,
                                alpha0_suboxic = rem_parms["alpha0"]["suboxic"],
                                alpha0_anoxic = rem_parms["alpha0"]["anoxic"],
                                alpha0_deeply_anoxic =rem_parms["alpha0"]["deeply_anoxic"],
                                pH0_suboxic =rem_parms["pH0"]["suboxic"],
                                pH0_anoxic =rem_parms["pH0"]["anoxic"],
                                pH0_deeply_anoxic =rem_parms["pH0"]["deeply_anoxic"],
                                mu1_suboxic = rem_parms["mu1"]["suboxic"],
                                mu1_anoxic = rem_parms["mu1"]["anoxic"],
                                mu1_deeply_anoxic = rem_parms["mu1"]["deeply_anoxic"],
                                organism_diam = rem_parms["organism_diam"]
                                )

        # Calculate advective microbial removal
        modpath_removal = TR.Transport(modpath_hor,
                                pollutant = organism)
    
        # Calculate advective microbial removal
        # Final concentration per endpoint_id
        C_final = {}
        # Start conc
        conc_start = 1.
        for endpoint_id in modpath_hor.schematisation_dict.get("endpoint_id"):
            df_particle[iDist], df_flowline[iDist], C_final[endpoint_id] = modpath_removal.calc_advective_microbial_removal(
                                                modpath_hor.df_particle, modpath_hor.df_flowline, 
                                                endpoint_id = endpoint_id,
                                                trackingdirection = modpath_hor.trackingdirection,
                                                conc_start = conc_start, conc_gw = 0.)

            # print("Final concentration at " + str(iDist) + " m is: " + str(round(C_final[endpoint_id],4)))
            # Add final concentration to summary dataframe
            df_conc.loc[iDist,"Final_concentration"] = C_final[endpoint_id]

        # df_particle file name 
        particle_fname = os.path.join(modpath_hor.dstroot,modpath_hor.schematisation_type + "_df_particle_microbial_removal" + str(iDist) + "m.csv")
        # Save df_particle 
        df_particle[iDist].to_csv(particle_fname)
        
        # df_flowline file name
        flowline_fname = os.path.join(modpath_hor.dstroot,modpath_hor.schematisation_type + "_df_flowline_microbial_removal" + str(iDist) + "m.csv")
        # Save df_flowline
        df_flowline[iDist].to_csv(flowline_fname)

    # Save summary dataframe of final concentration (dependent on boudary type)
    summary_conc_fname = os.path.join(modpath_hor.dstroot,"Final_concentrations.csv")
    df_conc.to_csv(summary_conc_fname)

    assert modpath_hor.success_mp

    # relative tolerance
    rtol = 5.e-3
    # absolute tollerance
    atol = rtol * conc_start

    # Assert that dataframes are equal
    assert_frame_equal(df_conc.loc[df_conc.index.isin(df_conc_horflow_test.index),:],
                        df_conc_horflow_test.loc[df_conc_horflow_test.index.isin(df_conc.index),:],
                        check_dtype=False,rtol = rtol,atol = atol)

# def test_modpath_run_horizontal_flow_diffuse(organism_name = "MS2"):

#     ''' Horizontal flow test in target_aquifer: modpath run.'''

#     # @Steven: deze test in read the docs zodat je zelf kunt visualisren of code klopt
#     # Test mag hier weg (geen assert)

#     # well discharge
#     well_discharge = -1000.
#     # distance to boundary
#     distance_boundary = 50.
#     # Center depth of target aquifer
#     z_point = 0 - 0.1 - 10.
#     x_point = distance_boundary - 0.5
#     # Phreatic scheme without gravelpack: modpath run.
#     test_conf_hor = AW.HydroChemicalSchematisation(schematisation_type='semiconfined',
#                                 computation_method = 'modpath',
#                                 what_to_export='all',
#                                 removal_function = 'mbo',
#                                 # biodegradation_sorbed_phase = False,
#                                 well_discharge= well_discharge,
#                                 # vertical_resistance_shallow_aquifer=500,
#                                 porosity_vadose_zone=0.33,
#                                 porosity_shallow_aquifer=0.33,
#                                 porosity_target_aquifer=0.33,
#                                 recharge_rate=0.,
#                                 moisture_content_vadose_zone=0.15,
#                                 ground_surface = 0.0,
#                                 thickness_vadose_zone_at_boundary=0.,
#                                 thickness_shallow_aquifer=0.1,
#                                 thickness_target_aquifer=20.0,
#                                 hor_permeability_target_aquifer=10.0,
#                                 hor_permeability_shallow_aquifer = 1.,
#                                 vertical_anisotropy_shallow_aquifer = 1000.,
#                                 thickness_full_capillary_fringe=0.4,
#                                 redox_vadose_zone='anoxic', #'suboxic',
#                                 redox_shallow_aquifer='anoxic',
#                                 redox_target_aquifer='anoxic',
#                                 pH_vadose_zone=7.5,
#                                 pH_shallow_aquifer=7.5,
#                                 pH_target_aquifer=7.5,
#                                 dissolved_organic_carbon_vadose_zone=1., 
#                                 dissolved_organic_carbon_shallow_aquifer=1., 
#                                 dissolved_organic_carbon_target_aquifer=1.,
#                                 fraction_organic_carbon_vadose_zone=0.001,
#                                 fraction_organic_carbon_shallow_aquifer=0.001,
#                                 fraction_organic_carbon_target_aquifer=0.001, 
#                                 temp_water=12.,
#                                 solid_density_vadose_zone= 2.650, 
#                                 solid_density_shallow_aquifer= 2.650, 
#                                 solid_density_target_aquifer= 2.650, 
#                                 diameter_borehole = 0.2,
#                                 name = organism_name,
#                                 diameter_filterscreen = 0.2,
#                                 top_clayseal = 0,
#                                 compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),
#                                 model_radius = distance_boundary,
#                                 ncols_near_well = 10,
#                                 ncols_far_well = int(distance_boundary/5),
#                                 # Point contamination
#                                 point_input_concentration = 1.,
#                                 discharge_point_contamination = abs(well_discharge),
#                                 distance_point_contamination_from_well = x_point,
#                                 depth_point_contamination = z_point,
#                                 )

#     test_conf_hor.make_dictionary()

    
#     # Remove "vadose" layer from geo_parameters
#     test_conf_hor.geo_parameters.pop("vadose_zone")

#     ### Adjust ibound_parameters to add horizontal flow ###
    
#     # Confined top boundary ; no recharge_parameters
#     test_conf_hor.ibound_parameters.pop("top_boundary1")
#     test_conf_hor.ibound_parameters.pop("inner_boundary_shallow_aquifer")
    
#     # Width (delr) outer boundary cells
#     delr_outer = 1.

#     # Add outer boundary for horizontal flow test
#     test_conf_hor.ibound_parameters["outer_boundary_target_aquifer"] = {
#                             'top': test_conf_hor.bottom_shallow_aquifer,
#                             'bot': test_conf_hor.bottom_target_aquifer,
#                             'xmin': test_conf_hor.model_radius - delr_outer,
#                             'xmax': test_conf_hor.model_radius,
#                             'ibound': -1
#                             }
#     # concentration_boundary_parameters
#     test_conf_hor.concentration_boundary_parameters = {'source1': {
#                 'organism_name': organism_name,
#                 'xmin': test_conf_hor.model_radius - delr_outer * 0.5,
#                 'xmax': test_conf_hor.model_radius,
#                 'zmax': test_conf_hor.bottom_shallow_aquifer,
#                 'zmin': test_conf_hor.bottom_target_aquifer,
#                 'input_concentration': test_conf_hor.diffuse_input_concentration
#                 }
#             # 'source2' :{}> surface water (BAR & RBF) #@MartinvdS come back to this when we start this module
#             }
#     # Remove/empty point_parameters
#     test_conf_hor.point_parameters = {}

#     # Add point parameters                       
#     startpoint_id = ["outer_boundary_target_aquifer"]
    

#     # ModPath well
#     modpath_hor = mpw.ModPathWell(test_conf_hor, #test_phrea,
#                             workspace = "test4_conf_hor_diffuse",
#                             modelname = "confined_hor",
#                             bound_left = "xmin",
#                             bound_right = "xmax")
#     # print(modpath_phrea.__dict__)
#     # Run phreatic schematisation
#     modpath_hor.run_model(run_mfmodel = True,
#                         run_mpmodel = True)

#     # microbial removal properties
#     # organism_name = 'MS2'
#     alpha0 = {"suboxic": 1.e-3, "anoxic": 1.e-5, "deeply_anoxic": 1.e-5}
#     pH0 = {"suboxic": 6.6, "anoxic": 6.8, "deeply_anoxic": 6.8}
#     organism_diam =  2.33e-8
#     mu1 = {"suboxic": 0.149,"anoxic": 0.023,"deeply_anoxic": 0.023}

#     removal_parameters = {organism_name: 
#                     {"organism_name": organism_name,
#                         "alpha0": alpha0,
#                         "pH0": pH0,
#                         "organism_diam": organism_diam,
#                         "mu1": mu1
#                     }
#                 }
#     # Removal parameters organism
#     rem_parms = removal_parameters[organism_name]

    # Pollutant to define removal parameters for
    organism = TR.MicrobialOrganism(organism_name=organism_name)
    # Calculate advective microbial removal
    modpath_removal = TR.Transport(modpath_hor,
                                    pollutant = organism)
 
#     # Calculate advective microbial removal
#     # Final concentration per endpoint_id
#     C_final = {}
#     for endpoint_id in modpath_hor.schematisation_dict.get("endpoint_id"):
#         df_particle, df_flowline, C_final[endpoint_id] = modpath_removal.calc_advective_microbial_removal(
#                                             modpath_hor.df_particle, modpath_hor.df_flowline, 
#                                             endpoint_id = endpoint_id,
#                                             trackingdirection = modpath_hor.trackingdirection,
#                                             conc_start = 1., conc_gw = 0.)

#     # df_particle file name 
#     particle_fname = os.path.join(modpath_hor.dstroot,modpath_hor.schematisation_type + "_df_particle_microbial_removal.csv")
#     # Save df_particle 
#     df_particle.to_csv(particle_fname)
    
#     # df_flowline file name
#     flowline_fname = os.path.join(modpath_hor.dstroot,modpath_hor.schematisation_type + "_df_flowline_microbial_removal.csv")
#     # Save df_flowline
#     df_flowline.to_csv(flowline_fname)


#     # Create travel time plots
#     fpath_scatter_times_log = os.path.join(modpath_hor.dstroot,"log_travel_times_" + endpoint_id + ".png")
#     fpath_scatter_times = os.path.join(modpath_hor.dstroot,"travel_times_" + endpoint_id + ".png")
#     # # df particle
#     # df_particle = modpath_removal.df_particle
#     # time limits
#     tmin, tmax = 0.1, 10000.
#     # xcoord bounds
#     xmin, xmax = 0., min(50., x_point)

#     # Create travel time plots (lognormal)
#     modpath_removal.plot_age_distribution(df_particle=df_particle,
#             vmin = tmin,vmax = tmax,
#             fpathfig = fpath_scatter_times_log, figtext = None,x_text = 0,
#             y_text = 0, lognorm = True, xmin = xmin, xmax = xmax,
#             line_dist = 1, dpi = 192, trackingdirection = "forward",
#             cmap = 'viridis_r')

#     # Create travel time plots (linear)
#     modpath_removal.plot_age_distribution(df_particle=df_particle,
#             vmin = 0.,vmax = tmax,
#             fpathfig = fpath_scatter_times, figtext = None,x_text = 0,
#             y_text = 0, lognorm = False, xmin = xmin, xmax = xmax,
#             line_dist = 1, dpi = 192, trackingdirection = "forward",
#             cmap = 'viridis_r')

#     # Create concentration plots
#     fpath_scatter_removal_log = os.path.join(modpath_hor.dstroot,"log_removal_" + endpoint_id + ".png")

#     # relative conc limits
#     cmin, cmax = 1.e-21, 1.
#     # xcoord bounds
#     xmin, xmax = 0., min(50., x_point)

#     # Create travel time plots (lognormal)
#     modpath_removal.plot_logremoval(df_particle=df_particle,
#             df_flowline=df_flowline,
#             vmin = cmin,vmax = cmax,
#             fpathfig = fpath_scatter_removal_log, figtext = None,x_text = 0,
#             y_text = 0, lognorm = True, xmin = xmin, xmax = xmax,
#             line_dist = 1, dpi = 192, trackingdirection = "forward",
#             cmap = 'viridis_r')

#     assert modpath_hor.success_mp

#%%
# def test_modpath_run_phreatic_withgravelpack_traveltimes(organism_name = "MS2"):
#     ''' Phreatic scheme with gravelpack: modpath run.'''

#     #@ Steven: test mag hier weg; documentatie mag in readthedocs structuur
#     test_phrea_gp = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
#                                 computation_method = 'modpath',
#                                 what_to_export='all',
#                                 removal_function = 'mbo',
#                                 # biodegradation_sorbed_phase = False,
#                                 well_discharge=-319.4*24,
#                                 # vertical_resistance_shallow_aquifer=500,
#                                 porosity_vadose_zone=0.38,
#                                 porosity_shallow_aquifer=0.35,
#                                 porosity_target_aquifer=0.35,
#                                 recharge_rate=0.3/365.25,
#                                 moisture_content_vadose_zone=0.15,
#                                 ground_surface = 22.0,
#                                 thickness_vadose_zone_at_boundary=5.0,
#                                 thickness_shallow_aquifer=10.0,
#                                 thickness_target_aquifer=40.0,
#                                 hor_permeability_target_aquifer=35.0,
#                                 hor_permeability_shallow_aquifer = 0.02,
#                                 thickness_full_capillary_fringe=0.4,
#                                 redox_vadose_zone='anoxic', #'suboxic',
#                                 redox_shallow_aquifer='anoxic',
#                                 redox_target_aquifer='deeply_anoxic',
#                                 pH_vadose_zone=5.,
#                                 pH_shallow_aquifer=6.,
#                                 pH_target_aquifer=7.,
#                                 dissolved_organic_carbon_vadose_zone=10., 
#                                 dissolved_organic_carbon_shallow_aquifer=4., 
#                                 dissolved_organic_carbon_target_aquifer=2.,
#                                 fraction_organic_carbon_vadose_zone=0.001,
#                                 fraction_organic_carbon_shallow_aquifer=0.0005,
#                                 fraction_organic_carbon_target_aquifer=0.0005, 
#                                 temp_water=11.,
#                                 solid_density_vadose_zone= 2.650, 
#                                 solid_density_shallow_aquifer= 2.650, 
#                                 solid_density_target_aquifer= 2.650, 
#                                 diameter_borehole = 0.75,
#                                 name = organism_name,
#                                 diameter_filterscreen = 0.2,
#                                 point_input_concentration = 100.,
#                                 discharge_point_contamination = 100.,#made up value
#                                 top_clayseal = 17,
#                                 compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),
#                                 # Modpath grid parms
#                                 ncols_near_well = 20,
#                                 ncols_far_well = 80,
#                             )

#     test_phrea_gp.make_dictionary()

#     # Remove/empty point_parameters
#     test_phrea_gp.point_parameters = {}

#     # print(test_phrea.__dict__)
#     modpath_phrea = mpw.ModPathWell(test_phrea_gp, #test_phrea,
#                             workspace = "test5_phrea_gp",
#                             modelname = "phreatic",
#                             bound_left = "xmin",
#                             bound_right = "xmax")
#     # print(modpath_phrea.__dict__)
#     # Run phreatic schematisation
#     modpath_phrea.run_model(run_mfmodel = True,
#                         run_mpmodel = True)
    
#     # Create travel time plots
#     fpath_scatter_times_log = os.path.join(modpath_phrea.dstroot,"log_travel_times_test.png")
#     fpath_scatter_times = os.path.join(modpath_phrea.dstroot,"travel_times_test.png")
#     # df particle
#     df_particle = modpath_phrea.df_particle
#     # time limits
#     tmin, tmax = 0.1, 10000.
#     # xcoord bounds
#     xmin, xmax = 0., 50.

#     # Create travel time plots (lognormal)
#     modpath_phrea.plot_pathtimes(df_particle = df_particle, 
#             vmin = tmin,vmax = tmax,
#             fpathfig = fpath_scatter_times_log, figtext = None,x_text = 0,
#             y_text = 0, lognorm = True, xmin = xmin, xmax = xmax,
#             line_dist = 1, dpi = 192, trackingdirection = "forward",
#             cmap = 'viridis_r')

#     # Create travel time plots (linear)
#     modpath_phrea.plot_pathtimes(df_particle = df_particle, 
#             vmin = 0.,vmax = tmax,
#             fpathfig = fpath_scatter_times, figtext = None,x_text = 0,
#             y_text = 0, lognorm = False, xmin = xmin, xmax = xmax,
#             line_dist = 1, dpi = 192, trackingdirection = "forward",
#             cmap = 'viridis_r')

#     assert modpath_phrea.success_mp

#%%
def test_modpath_run_phreatic_withgravelpack_removal(organism_name = "MS2"):
    ''' Phreatic scheme with gravelpack: modpath run.'''

    # @Steven: code met plots e.d. graag elders.
    # @Steven: Voeg assert toe in relatie tot final concentration == value X....

    test_phrea_gp = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                computation_method = 'modpath',
                                what_to_export='all',
                                removal_function = 'mbo',
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
                                temp_water=11.,
                                solid_density_vadose_zone= 2.650, 
                                solid_density_shallow_aquifer= 2.650, 
                                solid_density_target_aquifer= 2.650, 
                                diameter_borehole = 0.75,
                                name = organism_name,
                                diameter_filterscreen = 0.2,
                                point_input_concentration = 100.,
                                discharge_point_contamination = 100.,#made up value
                                top_clayseal = 17,
                                compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),
                                # Modpath grid parms
                                ncols_near_well = 20,
                                ncols_far_well = 80,
                            )

    test_phrea_gp.make_dictionary()

    # Remove/empty point_parameters
    test_phrea_gp.point_parameters = {}

    # print(test_phrea.__dict__)
    modpath_phrea = mpw.ModPathWell(test_phrea_gp, #test_phrea,
                            workspace = os.path.join(path,"test6_phrea_gp_removal"),
                            modelname = "phreatic",
                            bound_left = "xmin",
                            bound_right = "xmax")
    # print(modpath_phrea.__dict__)
    # Run phreatic schematisation
    modpath_phrea.run_model(run_mfmodel = True,
                        run_mpmodel = True)

    # microbial removal properties
    # organism_name = 'MS2'
    alpha0 = {"suboxic": 1.e-3, "anoxic": 1.e-5, "deeply_anoxic": 1.e-5}
    pH0 = {"suboxic": 6.6, "anoxic": 6.8, "deeply_anoxic": 6.8}
    organism_diam =  2.33e-8
    mu1 = {"suboxic": 0.149,"anoxic": 0.023,"deeply_anoxic": 0.023}

    removal_parameters = {organism_name: 
                    {"organism_name": organism_name,
                        "alpha0": alpha0,
                        "pH0": pH0,
                        "organism_diam": organism_diam,
                        "mu1": mu1
                    }
                }
    # Removal parameters organism
    rem_parms = removal_parameters[organism_name]

    # Pollutant to define removal parameters for
    organism = TR.MicrobialOrganism(organism_name=organism_name)

    # Calculate advective microbial removal
    modpath_removal = TR.Transport(modpath_phrea,
                            pollutant = organism)
 
    # Calculate advective microbial removal
    # Final concentration per endpoint_id
    C_final = {}
    for endpoint_id in modpath_phrea.schematisation_dict.get("endpoint_id"):
        df_particle, df_flowline, C_final[endpoint_id] = modpath_removal.calc_advective_microbial_removal(
                                            modpath_phrea.df_particle, modpath_phrea.df_flowline, 
                                            endpoint_id = endpoint_id,
                                            trackingdirection = modpath_phrea.trackingdirection,
                                            conc_start = 1., conc_gw = 0.)
    
        # Create travel time plots
        fpath_scatter_times_log = os.path.join(modpath_phrea.dstroot,"log_travel_times_" + endpoint_id + ".png")
        fpath_scatter_times = os.path.join(modpath_phrea.dstroot,"travel_times_" + endpoint_id + ".png")
        # # df particle
        # df_particle = modpath_removal.df_particle
        # time limits
        tmin, tmax = 0.1, 10000.
        # xcoord bounds
        xmin, xmax = 0., 50.

        # Create travel time plots (lognormal)
        modpath_removal.plot_age_distribution(df_particle=df_particle,
                vmin = tmin,vmax = tmax,
                fpathfig = fpath_scatter_times_log, figtext = None,x_text = 0,
                y_text = 0, lognorm = True, xmin = xmin, xmax = xmax,
                line_dist = 1, dpi = 192, trackingdirection = "forward",
                cmap = 'viridis_r')

        # Create travel time plots (linear)
        modpath_removal.plot_age_distribution(df_particle=df_particle,
                vmin = 0.,vmax = tmax,
                fpathfig = fpath_scatter_times, figtext = None,x_text = 0,
                y_text = 0, lognorm = False, xmin = xmin, xmax = xmax,
                line_dist = 1, dpi = 192, trackingdirection = "forward",
                cmap = 'viridis_r')

        # Create concentration plots
        fpath_scatter_removal_log = os.path.join(modpath_phrea.dstroot,"log_removal_" + endpoint_id + ".png")

        # relative conc limits
        cmin, cmax = 1.e-21, 1.
        # xcoord bounds
        xmin, xmax = 0., 50.

        # Create travel time plots (lognormal)
        modpath_removal.plot_logremoval(df_particle=df_particle,
                df_flowline=df_flowline,
                vmin = cmin,vmax = cmax,
                fpathfig = fpath_scatter_removal_log, figtext = None,x_text = 0,
                y_text = 0, lognorm = True, xmin = xmin, xmax = xmax,
                line_dist = 1, dpi = 192, trackingdirection = "forward",
                cmap = 'viridis_r')

    assert modpath_phrea.success_mp


#%%

def test_modpath_run_semiconfined_nogravelpack_traveltimes(organism_name = "MS2"):
    ''' Phreatic scheme with gravelpack: modpath run.'''

    # @Steven: test overbodig? --> plots en andere uitvoer hier weg
    # Voeg assert toe in relatie tot final concentration == value X....

    test_semiconf = AW.HydroChemicalSchematisation(schematisation_type='semiconfined',
                                    computation_method = 'modpath',
                                    what_to_export='all',
                                    removal_function = 'mbo',
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
                                      temp_water=11.,
                                      solid_density_vadose_zone= 2.650, 
                                      solid_density_shallow_aquifer= 2.650, 
                                      solid_density_target_aquifer= 2.650, 
                                      diameter_borehole = 0.75,
                                      name = organism_name,
                                      # diameter_filterscreen = 0.2,
                                      point_input_concentration = 100.,
                                      discharge_point_contamination = 100.,#made up value
                                      top_clayseal = 17,
                                      compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),
                                      
                                      ncols_near_well = 20,
                                      ncols_far_well = 20,
                                      nlayers_target_aquifer = 20,

                                    )

    test_semiconf.make_dictionary()

    # Test concentration_boundary_parameters (check line 796 in ModPath_functions)
    concentration_boundary_parameters = test_semiconf.recharge_parameters

    # Remove/empty point_parameters
    test_semiconf.point_parameters = {}

    modpath_semiconf = mpw.ModPathWell(test_semiconf, # semiconf_dict_1, #test_phrea,
                            workspace = os.path.join(path,"test7_semiconf_nogp"),
                            modelname = "semi_conf_nogp",
                            bound_left = "xmin",
                            bound_right = "xmax")

    # Run phreatic schematisation
    modpath_semiconf.run_model(run_mfmodel = True,
                        run_mpmodel = True)

    # Pollutant to define removal parameters for
    organism = TR.MicrobialOrganism(organism_name=organism_name)

    # Calculate advective microbial removal
    modpath_removal = TR.Transport(modpath_semiconf,
                                            pollutant = organism,
                                            )
    # modpath_removal.compute_omp_removal()
    # Final concentration per endpoint_id
    C_final = {}
    for endpoint_id in modpath_semiconf.schematisation_dict.get("endpoint_id"):
        df_particle, df_flowline, C_final[endpoint_id] = modpath_removal.calc_advective_microbial_removal(
                                            modpath_semiconf.df_particle, modpath_semiconf.df_flowline, 
                                            endpoint_id = endpoint_id,
                                            trackingdirection = modpath_semiconf.trackingdirection,
                                            mu1 = 0.023, grainsize = 0.00025, alpha0 = 1.E-5, pH0 = 6.8,
                                            temp_water = 11., rho_water = 999.703, organism_diam = 2.33e-8,
                                            conc_start = 1., conc_gw = 0.)

    # df_particle file name 
    particle_fname = os.path.join(modpath_semiconf.dstroot,modpath_semiconf.schematisation_type + "_df_particle_microbial_removal.csv")
    # Save df_particle 
    df_particle.to_csv(particle_fname)
    
    # df_flowline file name
    flowline_fname = os.path.join(modpath_semiconf.dstroot,modpath_semiconf.schematisation_type + "_df_flowline_microbial_removal.csv")
    # Save df_flowline
    df_flowline.to_csv(flowline_fname)

    
    # Create travel time plots
    fpath_scatter_times_log = os.path.join(modpath_semiconf.dstroot,"log_travel_times_test.png")
    fpath_scatter_times = os.path.join(modpath_semiconf.dstroot,"travel_times_test.png")
    # df particle
    df_particle = modpath_semiconf.df_particle
    # time limits
    tmin, tmax = 0.1, 10000.
    # xcoord bounds
    xmin, xmax = 0., 100.

    # Create travel time plots (lognormal)
    modpath_semiconf.plot_pathtimes(df_particle = df_particle, 
            vmin = tmin,vmax = tmax,
            fpathfig = fpath_scatter_times_log, figtext = None,x_text = 0,
            y_text = 0, lognorm = True, xmin = xmin, xmax = xmax,
            line_dist = 1, dpi = 192, trackingdirection = "forward",
            cmap = 'viridis_r')

    # Create travel time plots (linear)
    modpath_semiconf.plot_pathtimes(df_particle = df_particle, 
            vmin = 0.,vmax = tmax,
            fpathfig = fpath_scatter_times, figtext = None,x_text = 0,
            y_text = 0, lognorm = False, xmin = xmin, xmax = xmax,
            line_dist = 1, dpi = 192, trackingdirection = "forward",
            cmap = 'viridis_r')

    assert modpath_semiconf.success_mp

# #%%

# def test_diffuse_modpath_run_semiconfined_nogravelpack_traveltimes(organism_name = "MS2"):

#     ''' Phreatic scheme without gravelpack: modpath run.'''

#     # @Steven: test mag weg --> komt overeen met test 'test_modpath_run_semiconfined_nogravelpack_traveltimes'
#     test_semiconf = AW.HydroChemicalSchematisation(schematisation_type='semiconfined',
#                                     computation_method = 'modpath',
#                                     what_to_export='all',
#                                     removal_function = 'mbo',
#                                     # biodegradation_sorbed_phase = False,
#                                       well_discharge=-319.4*24,
#                                       # vertical_resistance_shallow_aquifer=500,
#                                       porosity_vadose_zone=0.38,
#                                       porosity_shallow_aquifer=0.35,
#                                       porosity_target_aquifer=0.35,
#                                       recharge_rate=0.3/365.25,
#                                       moisture_content_vadose_zone=0.15,
#                                       ground_surface = 22.0,
#                                       thickness_vadose_zone_at_boundary=5.0,
#                                       thickness_shallow_aquifer=10.0,
#                                       thickness_target_aquifer=40.0,
#                                       hor_permeability_target_aquifer=35.0,
#                                       hor_permeability_shallow_aquifer = 0.02,
#                                       thickness_full_capillary_fringe=0.4,
#                                       redox_vadose_zone='anoxic', #'suboxic',
#                                       redox_shallow_aquifer='anoxic',
#                                       redox_target_aquifer='deeply_anoxic',
#                                       pH_vadose_zone=5.,
#                                       pH_shallow_aquifer=6.,
#                                       pH_target_aquifer=7.,
#                                       dissolved_organic_carbon_vadose_zone=10., 
#                                       dissolved_organic_carbon_shallow_aquifer=4., 
#                                       dissolved_organic_carbon_target_aquifer=2.,
#                                       fraction_organic_carbon_vadose_zone=0.001,
#                                       fraction_organic_carbon_shallow_aquifer=0.0005,
#                                       fraction_organic_carbon_target_aquifer=0.0005, 
#                                       diffuse_input_concentration = 100, #ug/L
#                                       temp_water=11.,
#                                       solid_density_vadose_zone= 2.650, 
#                                       solid_density_shallow_aquifer= 2.650, 
#                                       solid_density_target_aquifer= 2.650, 
#                                       diameter_borehole = 0.75,
#                                       name = organism_name,
#                                       # diameter_filterscreen = 0.2,
#                                       point_input_concentration = 100.,
#                                       discharge_point_contamination = 100.,#made up value
#                                       top_clayseal = 17,
#                                       compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),
#                                       ncols_near_well = 20,
#                                       ncols_far_well = 20,
#                                       nlayers_target_aquifer = 20,
#                                     )

#     test_semiconf.make_dictionary()

#     # Remove/empty point_parameters
#     test_semiconf.point_parameters = {}

#     modpath_semiconf = mpw.ModPathWell(test_semiconf, #test_phrea,
#                             workspace = "test8_semiconf",
#                             modelname = "semi_conf_no_gp",
#                             bound_left = "xmin",
#                             bound_right = "xmax")

#     # Run phreatic schematisation
#     modpath_semiconf.run_model(run_mfmodel = True,
#                         run_mpmodel = True)
    
#     # Create travel time plots
#     fpath_scatter_times_log = os.path.join(modpath_semiconf.dstroot,"log_travel_times_test.png")
#     fpath_scatter_times = os.path.join(modpath_semiconf.dstroot,"travel_times_test.png")
#     # df particle
#     df_particle = modpath_semiconf.df_particle
#     # time limits
#     tmin, tmax = 0.1, 10000.
#     # xcoord bounds
#     xmin, xmax = 0., 100.

#     # Create travel time plots (lognormal)
#     modpath_semiconf.plot_pathtimes(df_particle = df_particle, 
#             vmin = tmin,vmax = tmax,
#             fpathfig = fpath_scatter_times_log, figtext = None,x_text = 0,
#             y_text = 0, lognorm = True, xmin = xmin, xmax = xmax,
#             line_dist = 1, dpi = 192, trackingdirection = "forward",
#             cmap = 'viridis_r')

#     # Create travel time plots (linear)
#     modpath_semiconf.plot_pathtimes(df_particle = df_particle, 
#             vmin = 0.,vmax = tmax,
#             fpathfig = fpath_scatter_times, figtext = None,x_text = 0,
#             y_text = 0, lognorm = False, xmin = xmin, xmax = xmax,
#             line_dist = 1, dpi = 192, trackingdirection = "forward",
#             cmap = 'viridis_r')

#     assert modpath_semiconf.success_mp

#%%

def test_phreatic_scheme_withgravelpack_dictinput(organism_name = "MS2"):
    ''' Check writing and reading of dictionary files: phreatic scheme.'''
    
    phreatic_scheme= AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                    computation_method = 'modpath',
                                    what_to_export='all',
                                    removal_function = 'mbo',
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
                                      temp_water=11.,
                                      solid_density_vadose_zone= 2.650, 
                                      solid_density_shallow_aquifer= 2.650, 
                                      solid_density_target_aquifer= 2.650, 
                                      diameter_borehole = 0.75,
                                      name = organism_name,
                                      diameter_filterscreen = 0.2,
                                      point_input_concentration = 100.,
                                      discharge_point_contamination = 100.,#made up value
                                      top_clayseal = 17,
                                      compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),
                                    )

    phreatic_scheme.make_dictionary()  
    
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

    fpath_research = os.path.abspath(os.path.join(path,os.pardir,"research"))
    fpath_phreatic_dict_check2 = os.path.join(fpath_research,"phreatic_dict_withgravel_test.txt")
    with open (fpath_phreatic_dict_check2, "w") as dict_file:
        dict_file.write(str(phreatic_dict_2))

    with open(fpath_phreatic_dict_check2,"r") as dict_file:
        dict_raw = dict_file.read()
        phreatic_dict_check2 = ast.literal_eval(dict_raw)  # ast --> abstract syntax trees
        # pd.read_csv(fpath_phreatic_dict_check2, delimiter=" ", header = None)

    assert phreatic_dict_2 == phreatic_dict_check2


#%%

def test_travel_time_distribution_phreatic_analytical_plus_modpath(organism_name = "MS2"):
    ''' Compare AnalyticalWell.py and ModpathWell.py travel times distribution.'''

    # @ Steven: plots naar readthedocs --> niet in test toevoegen
    # assert_frame_equal analytical vs Modpath df_uitvoer behouden
    
    output_phreatic = pd.read_csv(path / 'phreatic_test_mp.csv')
    output_phreatic = output_phreatic.round(7) #round to 7 digits (or any digit), keep same as for the output for the model to compare
    output_phreatic.index = RangeIndex(start = 1,stop=len(output_phreatic.index)+1)

    test_phrea = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                        computation_method= 'analytical',
                                        what_to_export='all',
                                        removal_function = 'mbo',
                                        well_discharge=-319.4*24,
                                        # vertical_resistance_shallow_aquifer=500,
                                        hor_permeability_shallow_aquifer = 35.,
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
                                        ncols_near_well = 20,
                                        ncols_far_well = 100,
                                        
                                        thickness_full_capillary_fringe=0.4,
                                        temp_water=11,
                                        solid_density_vadose_zone= 2.650,
                                        solid_density_shallow_aquifer= 2.650,
                                        solid_density_target_aquifer= 2.650,
                                        diameter_borehole = 0.75,
                                        name=organism_name
                                        )

    # AnalyticalWell object
    well1_an = AW.AnalyticalWell(test_phrea)
    well1_an.phreatic()
    df_particle_an = well1_an.df_particle
    df_particle_an = df_particle_an.round(7)

    # ModPath well object
    well1_mp = mpw.ModPathWell(test_phrea, #test_phrea,
                            workspace = os.path.join(path,"test11_phrea"),
                            modelname = "phreatic",
                            bound_left = "xmin",
                            bound_right = "xmax")

    # Run phreatic schematisation
    well1_mp.run_model(run_mfmodel = True,
                        run_mpmodel = True)

    df_particle_mp = well1_mp.df_particle
    # Obtain "flowline_id" from current index ;Add RangeIndex
    df_particle_mp["flowline_id"] = df_particle_mp.index
    df_particle_mp.index = RangeIndex(start = 1,stop=len(df_particle_mp.index)+1)
    # df_particle_mp = df_particle_mp[["total_travel_time", "travel_time_unsaturated",
    #                  "travel_time_shallow_aquifer", "travel_time_target_aquifer",
    #                  "radial_distance", ]]
    df_particle_mp = df_particle_mp.round(7)

    # Export output of analytical and modpath calculations
    
    # df_particle file name (analytical well data)
    particle_fname_an = os.path.join(well1_mp.dstroot,well1_mp.schematisation_type + "_df_particle_analytical.csv")
    # Save df_particle 
    df_particle_an.to_csv(particle_fname_an)
   
    # df_particle file name 
    particle_fname_mp = os.path.join(well1_mp.dstroot,well1_mp.schematisation_type + "_df_particle_modpath.csv")
    # Save df_particle 
    df_particle_mp.to_csv(particle_fname_mp)


    # Columns in summary dataframe traveltimes
    summary_columns = output_phreatic.columns
    # flowline indices
    flowline_ids_an = list(df_particle_an.loc[:,"flowline_id"].unique())
    # Zones to compare
    zones = ["vadose_zone","shallow_aquifer","target_aquifer"]
    
                    #  ["total_travel_time",
                    # "travel_time_vadose_zone",
                    # "travel_time_shallow_aquifer",
                    # "travel_time_target_aquifer",
                    # "xcoord"]

    # Summary traveltimes (analytical)
    summary_traveltimes_an = pd.DataFrame(index = flowline_ids_an,
                                        columns = summary_columns)
    # Fill summary df         
    for fid in flowline_ids_an:
        # Total traveltime
        summary_traveltimes_an.loc[fid,"total_travel_time"] = \
                                df_particle_an.loc[df_particle_an["flowline_id"] == fid,:].sort_values(by="total_travel_time", 
                                        ascending = True).loc[:,"total_travel_time"].iloc[-1]
        # Starting location
        summary_traveltimes_an.loc[fid,"xcoord"] = \
                        df_particle_an.loc[df_particle_an["flowline_id"] == fid,:].sort_values(by="total_travel_time", 
                                ascending = True).loc[:,"xcoord"].iloc[0]
        for iZone in zones:
            # Travel time per zone
            summary_traveltimes_an.loc[fid,"travel_time_" + iZone] = df_particle_an.loc[(df_particle_an["flowline_id"] == fid) & (df_particle_an["zone"] == iZone),"travel_time"].values[0]

    # Summary traveltimes (Numerical / modpath)
    # flowline indices (modpath)
    flowline_ids_mp = list(df_particle_mp.loc[:,"flowline_id"].unique())
    summary_traveltimes_mp = pd.DataFrame(index = flowline_ids_mp,
                                        columns = summary_columns)

    # Fill summary df         
    for fid in flowline_ids_mp:
        # Total traveltime
        summary_traveltimes_mp.loc[fid,"total_travel_time"] = \
                                df_particle_mp.loc[df_particle_mp["flowline_id"] == fid,:].sort_values(by="total_travel_time", 
                                        ascending = True).loc[:,"total_travel_time"].iloc[-1]
        # Starting location
        summary_traveltimes_mp.loc[fid,"xcoord"] = \
                        df_particle_mp.loc[df_particle_mp["flowline_id"] == fid,:].sort_values(by="total_travel_time", 
                                ascending = True).loc[:,"xcoord"].iloc[0]
        
        ## Travel time per zone ##
        # Vadose zone
        summary_traveltimes_mp.loc[fid,"travel_time_vadose_zone"] = abs(df_particle_mp.loc[(df_particle_mp["flowline_id"] == fid) & (df_particle_mp["zone"] == "shallow_aquifer"),"total_travel_time"].values.min() - \
                                                                    df_particle_mp.loc[(df_particle_mp["flowline_id"] == fid) & (df_particle_mp["zone"] == "vadose_zone"),"total_travel_time"].values.min())
        # Shallow aquifer
        summary_traveltimes_mp.loc[fid,"travel_time_shallow_aquifer"] = abs(df_particle_mp.loc[(df_particle_mp["flowline_id"] == fid) & (df_particle_mp["zone"] == "target_aquifer"),"total_travel_time"].values.min() - \
                                                                    df_particle_mp.loc[(df_particle_mp["flowline_id"] == fid) & (df_particle_mp["zone"] == "shallow_aquifer"),"total_travel_time"].values.min())
        # Target aquifer
        summary_traveltimes_mp.loc[fid,"travel_time_target_aquifer"] = abs(df_particle_mp.loc[(df_particle_mp["flowline_id"] == fid) & (df_particle_mp["zone"] == "well1"),"total_travel_time"].values.min() - \
                                                                    df_particle_mp.loc[(df_particle_mp["flowline_id"] == fid) & (df_particle_mp["zone"] == "target_aquifer"),"total_travel_time"].values.min())

    # df_particle file name (analytical well data)
    summary_fname_an = os.path.join(well1_mp.dstroot,well1_mp.schematisation_type + "_summary_traveltimes_analytical.csv")
    # Save df_particle 
    summary_traveltimes_an.to_csv(summary_fname_an)
   
    # df_particle file name 
    summary_fname_mp = os.path.join(well1_mp.dstroot,well1_mp.schematisation_type + "_summary_traveltimes_modpath.csv")
    # Save df_particle 
    summary_traveltimes_mp.to_csv(summary_fname_mp)

    # Pollutant to define removal parameters for
    organism = TR.MicrobialOrganism(organism_name=organism_name)

    # Create traveltime distribution plot using Substance Transport class
    modpath_removal = TR.Transport(well1_mp,
                            pollutant = organism)

    modpath_removal.plot_travel_time_distribution(modpath_removal.df_particle, times_col = "total_travel_time",
                                      distance_col = "xcoord", index_col = "flowline_id",
                                      fpath_fig = None)

    ## analytical solution ##
    # Distance points
    distance_points_an = np.empty((len(flowline_ids_an)), dtype = 'float')
    # Time points
    time_points_an = np.empty((len(flowline_ids_an)), dtype = 'float')

    for idx, fid in enumerate(flowline_ids_an):

        # Fill distance_points array
        distance_points_an[idx] = summary_traveltimes_an.loc[summary_traveltimes_an.index == fid, ["xcoord","total_travel_time"]].sort_values(
                                    by = "total_travel_time", ascending = True).iloc[-1,0]

        # Fill time_points array
        time_points_an[idx] = summary_traveltimes_an.loc[summary_traveltimes_an.index == fid, ["xcoord","total_travel_time"]].sort_values(
                                    by = "total_travel_time", ascending = True).iloc[-1,1]

    ## modpath solution ##
    # Distance points
    distance_points_mp = np.empty((len(flowline_ids_mp)), dtype = 'float')
    # Time points
    time_points_mp = np.empty((len(flowline_ids_mp)), dtype = 'float')

    for idx, fid in enumerate(flowline_ids_mp):

        # Fill distance_points array
        distance_points_mp[idx] = summary_traveltimes_mp.loc[summary_traveltimes_mp.index == fid, ["xcoord","total_travel_time"]].sort_values(
                                    by = "total_travel_time", ascending = True).iloc[-1,0]
          
        # Fill time_points array
        time_points_mp[idx] = summary_traveltimes_mp.loc[summary_traveltimes_mp.index == fid, ["xcoord","total_travel_time"]].sort_values(
                                    by = "total_travel_time", ascending = True).iloc[-1,1]
    

    # Create travel time distribution plot
    fig, ax = plt.subplots(figsize= (10,10),dpi=300)

    # Plot time-radius plot
    plt.plot(distance_points_an, time_points_an, lw = 0., marker = '.')  # analytical
    plt.plot(distance_points_mp, time_points_mp, lw = 0.2)  # modpath
    plt.xlabel("Radial distance (m)")
    plt.ylabel("Travel time (days)")

    # travel distribution plot (analytic vs numeric)
    fpath_traveltime_plot = os.path.join(well1_mp.dstroot,well1_mp.schematisation_type + "_traveltimes_modpath.png")
    fig.savefig(fpath_traveltime_plot)


    # Create travel time plots
    fpath_scatter_times_log = os.path.join(well1_mp.dstroot,"log_travel_times_test.png")
    fpath_scatter_times = os.path.join(well1_mp.dstroot,"travel_times_test.png")
    # df particle
    df_particle = well1_mp.df_particle
    # time limits
    tmin, tmax = 0.1, 10000.
    # xcoord bounds
    xmin, xmax = 0., 50.

    # Create travel time plots (lognormal)
    well1_mp.plot_pathtimes(df_particle = df_particle, 
            vmin = tmin,vmax = tmax,
            fpathfig = fpath_scatter_times_log, figtext = None,x_text = 0,
            y_text = 0, lognorm = True, xmin = xmin, xmax = xmax,
            line_dist = 1, dpi = 192, trackingdirection = "forward",
            cmap = 'viridis_r')

    # Create travel time plots (linear)
    well1_mp.plot_pathtimes(df_particle = df_particle, 
            vmin = 0.,vmax = tmax,
            fpathfig = fpath_scatter_times, figtext = None,x_text = 0,
            y_text = 0, lognorm = False, xmin = xmin, xmax = xmax,
            line_dist = 1, dpi = 192, trackingdirection = "forward",
            cmap = 'viridis_r')


    # try:
    assert_frame_equal(summary_traveltimes_an.loc[:, summary_columns],
                        output_phreatic.loc[:, summary_columns],check_dtype=False, check_index_type = False)

    # except AssertionError:
    #     print("Assertion Exception Raised - TTD test")
    # else:
    #     print("Success, no error in TTD!")
9
# #%%
# def test_phreatic_defecation_withgravelpack(organism_name = "MS2"):
#     ''' Scenario with feces at surface level near pumping well.
#         Case: phreatic
#         organism_name: "MS2" 
#         Vadose zone: False
#         redox: anoxic
#         pH: neutral (7.5)
#         depth of a leak 'leak1' at 9.5 m depth
        
#         ''' 
#     ''' Phreatic scheme with gravelpack: modpath run.'''
    
#     # @ Steven 2022-7-11: tijdelijk niet invoegen --> eerst readthedocs werkend.
    
#     test_feces_contamination = AW.HydroChemicalSchematisation(schematisation_type='semiconfined',
#                                 computation_method = 'modpath',
#                                 what_to_export='all',
#                                 removal_function = 'mbo',
#                                 # biodegradation_sorbed_phase = False,
#                                 well_discharge=-1000.,
#                                 # vertical_resistance_shallow_aquifer=500,
#                                 porosity_vadose_zone=0.33,
#                                 porosity_shallow_aquifer=0.33,
#                                 porosity_target_aquifer=0.33,
#                                 porosity_gravelpack=0.33,
#                                 porosity_clayseal=0.33,
#                                 recharge_rate=0.001,
#                                 moisture_content_vadose_zone=0.15,
#                                 ground_surface = 0.0,
#                                 thickness_vadose_zone_at_boundary=0.,
#                                 thickness_shallow_aquifer=30.,
#                                 thickness_target_aquifer=20.,
#                                 hor_permeability_target_aquifer= 10.0,
#                                 hor_permeability_shallow_aquifer = 10.,
#                                 hor_permeability_gravelpack= 1000.,
#                                 hor_permeability_clayseal= 0.005,
#                                 vertical_anisotropy_shallow_aquifer = 5.,
#                                 vertical_anisotropy_target_aquifer = 5.,
#                                 vertical_anisotropy_gravelpack = 1.,
#                                 vertical_anisotropy_clayseal = 1.,
#                                 thickness_full_capillary_fringe=0.,
#                                 grainsize_vadose_zone= 0.00025,
#                                 grainsize_shallow_aquifer=0.00025,
#                                 grainsize_target_aquifer=0.00025,
#                                 redox_vadose_zone='anoxic', 
#                                 redox_shallow_aquifer='anoxic',
#                                 redox_target_aquifer='anoxic',
#                                 pH_vadose_zone=7.5,
#                                 pH_shallow_aquifer=7.5,
#                                 pH_target_aquifer=7.5,
#                                 temp_water=11.,
#                                 solid_density_vadose_zone= 2.650, 
#                                 solid_density_shallow_aquifer= 2.650, 
#                                 solid_density_target_aquifer= 2.650, 
#                                 diameter_borehole = 0.75,
#                                 name = organism_name,
#                                 diameter_filterscreen = 0.2,
#                                 diameter_gravelpack = 0.75,
#                                 diameter_clayseal = 0.75,
#                                 point_input_concentration = 1.,
#                                 diffuse_input_concentration=2.15E6,
#                                 discharge_point_contamination = 100.,#made up value
#                                 top_clayseal = 0,
#                                 bottom_clayseal = -1.,
#                                 top_gravelpack = -1.,
#                                 grainsize_gravelpack=0.001,
#                                 grainsize_clayseal=0.000001,

#                                 compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),
#                                 # Modpath grid parms
#                                 ncols_near_well = 20,
#                                 ncols_far_well = 80,
#                                 nlayers_shallow_aquifer = 60,
#                                 nlayers_target_aquifer = 40,
#                             )

#     test_feces_contamination.make_dictionary()

#     # Remove/empty point_parameters
#     test_feces_contamination.point_parameters = {}

#     # Add leakage at 9.5 m depth
#     leak_Q = -1.
#     leak_depth = -9.5
#     leak_depth_bot = leak_depth - 0.1
#     test_feces_contamination.well_parameters = {
#                 'leak1': {
#                     'well_discharge': leak_Q,
#                     'top': leak_depth,
#                     'bot': leak_depth_bot,
#                     'xmin': 0., # test_feces_contamination.diameter_filterscreen/2,  # next to filter screen [ibound=0]
#                     'xmax': test_feces_contamination.diameter_filterscreen/2.,
#                     },
#                 }
#     # Change inner_boundary ibound
#     test_feces_contamination.ibound_parameters.pop('inner_boundary_shallow_aquifer')
#     test_feces_contamination.ibound_parameters['inner_boundary_shallow_aquifer_above_leak'] = {'top': 0.0, 'bot': leak_depth,
#                                                                                     'xmin': 0, 'xmax': test_feces_contamination.diameter_filterscreen/2., 'ibound': 0}
#     test_feces_contamination.ibound_parameters['inner_boundary_shallow_aquifer_below_leak'] = {'top': leak_depth_bot, 'bot': test_feces_contamination.bottom_shallow_aquifer,
#                                                                                     'xmin': 0, 'xmax': test_feces_contamination.diameter_filterscreen/2., 'ibound': 0}
    
#     # # Add key to geo_parameters
#     # test_feces_contamination.geo_parameters['leak1'] = {
#     #                 'top': leak_depth,
#     #                 'bot': leak_depth_bot,
#     #                 'xmin': 0.,  # next to filter screen [ibound=0]
#     #                 'xmax': test_feces_contamination.diameter_filterscreen/2.,
#     #                 },

#     # Add leak as endpoint id
#     test_feces_contamination.endpoint_id = copy.deepcopy(test_feces_contamination.well_parameters)
#     test_feces_contamination.endpoint_id['leak1'].pop('well_discharge')

#     # print(test_phrea.__dict__)
#     modpath_phrea = mpw.ModPathWell(test_feces_contamination, #test_phrea,
#                             workspace = "test12_defecation_human",
#                             modelname = "semiconfined",
#                             bound_left = "xmin",
#                             bound_right = "xmax")
#     # print(modpath_phrea.__dict__)
#     # Run phreatic schematisation
#     modpath_phrea.run_model(run_mfmodel = True,
#                         run_mpmodel = True)


#     # microbial removal properties
#     # organism_name = 'MS2'
#     alpha0 = {"suboxic": 1.e-3, "anoxic": 1.e-5, "deeply_anoxic": 1.e-5}
#     pH0 = {"suboxic": 6.6, "anoxic": 6.8, "deeply_anoxic": 6.8}
#     organism_diam =  2.33e-8
#     mu1 = {"suboxic": 0.149,"anoxic": 0.023,"deeply_anoxic": 0.023}

#     removal_parameters = {organism_name: 
#                     {"organism_name": organism_name,
#                         "alpha0": alpha0,
#                         "pH0": pH0,
#                         "organism_diam": organism_diam,
#                         "mu1": mu1
#                     }
#                 }
#     # Removal parameters organism
#     rem_parms = removal_parameters[organism_name]

    # # Pollutant to define removal parameters for
    # organism = TR.MicrobialOrganism(organism_name=organism_name,
    #                         alpha0_suboxic = rem_parms["alpha0"]["suboxic"],
    #                         alpha0_anoxic = rem_parms["alpha0"]["anoxic"],
    #                         alpha0_deeply_anoxic =rem_parms["alpha0"]["deeply_anoxic"],
    #                         pH0_suboxic =rem_parms["pH0"]["suboxic"],
    #                         pH0_anoxic =rem_parms["pH0"]["anoxic"],
    #                         pH0_deeply_anoxic =rem_parms["pH0"]["deeply_anoxic"],
    #                         mu1_suboxic = rem_parms["mu1"]["suboxic"],
    #                         mu1_anoxic = rem_parms["mu1"]["anoxic"],
    #                         mu1_deeply_anoxic = rem_parms["mu1"]["deeply_anoxic"],
    #                         organism_diam = rem_parms["organism_diam"]
    # )

    # # Calculate advective microbial removal
    # modpath_removal = TR.SubstanceTransport(modpath_phrea,
    #                         pollutant = organism)
 
#     # Calculate advective microbial removal
#     # Final concentration per endpoint_id
#     C_final = {}
#     for endpoint_id in modpath_phrea.schematisation_dict.get("endpoint_id"):
#         df_particle, df_flowline, C_final[endpoint_id] = modpath_removal.calc_advective_microbial_removal(
#                                             modpath_phrea.df_particle, modpath_phrea.df_flowline, 
#                                             endpoint_id = endpoint_id,
#                                             trackingdirection = modpath_phrea.trackingdirection,
#                                             conc_start = 1., conc_gw = 0.)


#         # df_particle file name 
#         particle_fname = os.path.join(modpath_phrea.dstroot,modpath_phrea.schematisation_type + "_df_particle_microbial_removal.csv")
#         # Save df_particle 
#         df_particle.to_csv(particle_fname)
        
#         # df_flowline file name
#         flowline_fname = os.path.join(modpath_phrea.dstroot,modpath_phrea.schematisation_type + "_df_flowline_microbial_removal.csv")
#         # Save df_flowline
#         df_flowline.to_csv(flowline_fname)


    
#         # Create travel time plots
#         fpath_scatter_times_log = os.path.join(modpath_phrea.dstroot,"log_travel_times_" + endpoint_id + ".png")
#         fpath_scatter_times = os.path.join(modpath_phrea.dstroot,"travel_times_" + endpoint_id + ".png")
#         # # df particle
#         # df_particle = modpath_removal.df_particle
#         # time limits
#         tmin, tmax = 0.1, 50000.
#         # xcoord bounds
#         xmin, xmax = 0., 25.

#         # Create travel time plots (lognormal)
#         modpath_removal.plot_age_distribution(df_particle=df_particle,
#                 vmin = tmin,vmax = tmax,
#                 fpathfig = fpath_scatter_times_log, figtext = None,x_text = 0,
#                 y_text = 0, lognorm = True, xmin = xmin, xmax = xmax,
#                 line_dist = 1, dpi = 192, trackingdirection = "forward",
#                 cmap = 'viridis_r')

#         # Create travel time plots (linear)
#         modpath_removal.plot_age_distribution(df_particle=df_particle,
#                 vmin = 0.,vmax = tmax,
#                 fpathfig = fpath_scatter_times, figtext = None,x_text = 0,
#                 y_text = 0, lognorm = False, xmin = xmin, xmax = xmax,
#                 line_dist = 1, dpi = 192, trackingdirection = "forward",
#                 cmap = 'viridis_r')

#         # Create concentration plots
#         fpath_scatter_removal_log = os.path.join(modpath_phrea.dstroot,"log_removal_" + endpoint_id + ".png")

#         # relative conc limits
#         cmin, cmax = 1.e-21, 1.
#         # xcoord bounds
#         xmin, xmax = 0., 25.

#         # Create travel time plots (lognormal)
#         modpath_removal.plot_logremoval(df_particle=df_particle,
#                 df_flowline=df_flowline,
#                 vmin = cmin,vmax = cmax,
#                 fpathfig = fpath_scatter_removal_log, figtext = None,x_text = 0,
#                 y_text = 0, lognorm = True, xmin = xmin, xmax = xmax,
#                 line_dist = 1, dpi = 192, trackingdirection = "forward",
#                 cmap = 'viridis_r')          

#     assert modpath_phrea.success_mp

def test_omp_removal_analyticalwell_input(substance_name = 'AMPA'):
    ''' calculate the removal of a default substance using AnalyticalWell class.'''

    # Lets start with a simple example defining a HydroChemicalSchematisation object for a phreatic aquifer:

    phreatic_schematisation = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                                        computation_method='analytical',
                                                        removal_function = 'omp',
                                                        well_discharge=-7500, #m3/day
                                                        recharge_rate=0.0008, #m/day
                                                        thickness_vadose_zone_at_boundary=5, #m
                                                        thickness_shallow_aquifer=10,  #m
                                                        thickness_target_aquifer=40, #m
                                                        hor_permeability_target_aquifer=35, #m/day
                                                        redox_vadose_zone='anoxic',
                                                        redox_shallow_aquifer='anoxic',
                                                        redox_target_aquifer='deeply_anoxic',
                                                        pH_target_aquifer=7.,
                                                        temp_water=11.,
                                                        name='benzene',
                                                        diffuse_input_concentration = 100, #ug/L
                                                        )

    # The parameters from the HydroChemicalSchematisation class are added as attributes to
    # the class and can be accessed for example:
    phreatic_schematisation.schematisation_type
    phreatic_schematisation.well_discharge
    phreatic_schematisation.porosity_shallow_aquifer

    # If not defined, default values are used for the rest of the parameters. To view all parameters in the schematisation:
    phreatic_schematisation.__dict__

    # Step 2: Run the AnalyticalWell class
    # =====================================
    # Next we create an AnalyticalWell object for the HydroChemicalSchematisation object we just made.
    phreatic_well = AW.AnalyticalWell(phreatic_schematisation)

    # Then we calculate the travel time for each of the zones unsaturated, shallow aquifer and target aquifer zones
    # by running the .phreatic() function for the well object. 
    phreatic_well.phreatic()

    # The total travel time can be plotted as a function of radial distance from the well, or as a function
    # of the cumulative fraction of abstracted water: 
    radial_plot = phreatic_well.plot_travel_time_versus_radial_distance(xlim=[0, 2000], ylim=[1e3, 1e6])
    cumulative_plot = phreatic_well.plot_travel_time_versus_cumulative_abstracted_water(xlim=[0, 1], ylim=[1e3, 1e6])

    # .. image:: https://github.com/KWR-Water/sutra2/blob/main/docs/_images/travel_time_versus_radial_distance_phreatic.png?raw=true
    #   :width: 600
    #   :alt: travel_time_versus_radial_distance_phreatic.png

    # .. image:: https://github.com/KWR-Water/sutra2/blob/main/docs/_images/travel_time_versus_cumulative_abs_water_phreatic.png?raw=true
    #   :width: 600
    #   :alt: travel_time_versus_cumulative_abs_water_phreatic.png

    # From the AnalyticalWell class two other important outputs are:

    # * df_particle - Pandas dataframe with data about the different flowlines per zone (unsaturated/shallow/target)
    # * df_flowline - Pandas dataframe with data about the flowlines per flowline (eg. total travel time per flowline)

    # Step 3: View the Substance class (Optional)
    # ===========================================
    # You can retrieve the default substance parameters used to calculate the removal in the
    # SubstanceTransport class. The data are stored in a dictionary
    test_substance = TR.Substance(substance_name=substance_name)
    test_substance.substance_dict


    # Step 4: Run the SubstanceTransport class
    # ========================================
    # To calculate the removal and the steady-state concentration in each zone, create a concentration
    # object by running the SubstanceTransport class with the phreatic_well object and specifying
    # the OMP (or pathogen) of interest.

    # In this example we use benzene. First we create the object and view the substance properties:
    phreatic_concentration = TR.Transport(well = phreatic_well, pollutant = test_substance)
    phreatic_concentration.removal_parameters

    # .. Optional: You may specify a different value for the substance parameters, for example
    # .. a different half-life for the anoxic redox zone. This can be input in the HydroChemicalSchematisation
    # .. and this will be used in the calculation for the removal for the OMP. The AnalyticalWell and 
    # .. phreatic() functions must be rerun:
       
    # Optional: You may specify a different value for the substance parameters, for example
    # a different half-life for the anoxic redox zone. This can be input in the SubstanceTransport
    # and this will be used in the calculation for the removal for the OMP. The SubstanceTransportclass
    # must be reloaded with the new input.
    test_substance2 = TR.Substance(substance_name='benzene',
                                                partition_coefficient_water_organic_carbon=2,
                                                molar_mass = 78.1,
                                                dissociation_constant=1,
                                                halflife_suboxic=12, 
                                                halflife_anoxic=420, 
                                                halflife_deeply_anoxic=6000)
    
    phreatic_concentration = TR.Transport(well = phreatic_well, 
                                            pollutant = test_substance2)

    # If you have specified values for the substance (e.g. half-life, pKa, log_Koc),
    # the default value is overriden and used in the calculation of the removal. You can
    # view the updated removal parameters ('substance dictionary') from the concentration object:
    phreatic_concentration.removal_parameters

    # Then we can compute the removal by running the 'compute_omp_removal' function:
    # phreatic_concentration.compute_omp_removal()
    phreatic_concentration.compute_omp_removal()

    # Once the removal has been calculated, you can view the steady-state concentration
    # and breakthrough time per zone for the OMP in the df_particle:
    phreatic_concentration.df_particle[['flowline_id', 'zone', 'steady_state_concentration', 'travel_time']].head(4)

    # View the steady-state concentration of the flowline or the steady-state
    # contribution of the flowline to the concentration in the well
    phreatic_concentration.df_flowline[['flowline_id', 'breakthrough_concentration', 'total_breakthrough_travel_time']].head(5)

def test_omp_removal_modpath_input(substance_name = 'AMPA'):
    ''' calculate the removal of a default substance using ModPathWell class.'''

    # Lets start with a simple example defining a HydroChemicalSchematisation object for a phreatic aquifer:

    phreatic_schematisation = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                                        computation_method = 'modpath',
                                                        removal_function = 'omp',
                                                        well_discharge=-7500, #m3/day
                                                        recharge_rate=0.0008, #m/day
                                                        thickness_vadose_zone_at_boundary=5, #m
                                                        thickness_shallow_aquifer=10,  #m
                                                        thickness_target_aquifer=40, #m
                                                        hor_permeability_target_aquifer=35, #m/day
                                                        redox_vadose_zone='anoxic',
                                                        redox_shallow_aquifer='anoxic',
                                                        redox_target_aquifer='deeply_anoxic',
                                                        pH_target_aquifer=7.,
                                                        temp_water=11.,
                                                        name='benzene',
                                                        diffuse_input_concentration = 100, #ug/L
                                                        )

    # Then, we create a ModpathWell object for the HydroChemicalSchematisation object that we just made.
    # The ModpathWell object requires a dictionary of the subsurface schematisation and a set of boundary conditions
    # the numerical model has to abide by in calculating flow velocity and direction of flow.
    phreatic_schematisation.make_dictionary()

    # Step 2: Run the ModpathWell class
    # =====================================
    # Next we create an ModpathWell object for the HydroChemicalSchematisation object we just made.
    # The data files will be stored in location workspace using a given modelname.
    modpath_phrea = mpw.ModPathWell(phreatic_schematisation,
                                workspace = os.path.join(path,"test7_omp_removal"),
                                modelname = "phreatic")

    # Now we run the Modpath model, which numerically calculates the flow in the subsurface using the 
    # 'schematisation' dictionary stored in the HydroChemicalSchematisation object. By default the model will
    # calculate both the hydraulic head distribution (using modflow: 'run_mfmodel' = True) and
    # the particle pathlines [X,Y,Z,T-data] (using modpath: 'run_mpmodel' = True) with which OMP removal
    # or microbial organism ('mbo') removal is later calculated.
    modpath_phrea.run_model(run_mfmodel = True,
                        run_mpmodel = True)

    # Step 3: Collect removal parameters
    # ===========================================

    # Step 3a: View the Substance class (Optional)
    # ============================================
    # You can retrieve the default removal parameters used to calculate the removal of organic micropollutants [OMP] 
    # in the SubstanceTransport class. The data are stored in a dictionary
    test_substance = TR.Substance(substance_name='benzene')
    test_substance.substance_dict

    # Step 4: Run the SubstanceTransport class
    # ========================================
    # To calculate the removal and the steady-state concentration in each zone (analytical solution) or per particle node (modpath), create a concentration
    # object by running the SubstanceTransport class with the phreatic_well object and specifying
    # the OMP or microbial organism (mbo) of interest. 
    # The type of removal is defined using the option 'removal_function: 'omp' or 'mbo'
    # All required parameters for removal are stored as 'removal_parameters'.

    # Step 4b: Calculate the OMP removal
    # ========================================
    # As example, we take the default removal parameters for the substances 'AMPA'.
    # Note: For OMP you will have to specify values relevant for substances (e.g. half-life, pKa, log_Koc).
    # Any/all default values will be stored and used in the calculation of the removal. 
    # Note that by default the class expects the removal of microbial organisms copied from removal_function 
    # entered in modpath_phrea. We have to explicitly enter the removal_function below for removal op substances.
    # removal_function == 'omp'

    # substance (AMPA)
    substance_name = 'AMPA'
    substance_default = TR.Substance(substance_name=substance_name,
                                    partition_coefficient_water_organic_carbon=None,
                                    molar_mass = None,
                                    dissociation_constant=None,
                                    halflife_suboxic=None,
                                    halflife_anoxic=None,
                                    halflife_deeply_anoxic=None)
    # Calculate removal of organic micro-pollutants (removal_function = 'omp')
    modpath_removal = TR.Transport(well = modpath_phrea,
                                    pollutant = substance_default)

    # View the updated removal_parameters dictionary from the SubstanceTransport object
    modpath_removal.removal_parameters

    # We compute the removal by running the 'compute_omp_removal' function:
    # modpath_removal.compute_omp_removal()
    modpath_removal.compute_omp_removal()

    # Once the removal has been calculated, you can view the steady-state concentration
    # and breakthrough time per zone for the OMP in the df_particle:
    modpath_removal.df_particle.loc[:,['zone', 'steady_state_concentration', 'travel_time']].head(4)

    # View the steady-state concentration of the flowline or the steady-state
    # contribution of the flowline to the concentration in the well
    modpath_removal.df_flowline.loc[:,['breakthrough_concentration', 'total_breakthrough_travel_time']].head(5)

#=======

#%%
# if __name__ == "__main__":
#     test_modpath_run_phreatic_withgravelpack()