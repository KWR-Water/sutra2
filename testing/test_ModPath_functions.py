
#%%
import pytest
from pandas import read_csv
import datetime as dt
import numpy as np
import pandas as pd
import os
import ast  # abstract syntax trees
import sys
# path = os.getcwd()  # path of working directory
from pathlib import Path

from zmq import zmq_version_info

    
import greta.Analytical_Well as AW
import greta.ModPath_functions as MP
import greta.Substance_Transport as ST

# from Substance_Transport import *
# from ModPath_functions import ModPathWell

from pandas._testing import assert_frame_equal
# from greta.ModPath_functions import ModPathWell  
# =======
# from greta.Analytical_Well import *
# from greta.Substance_Transport import *
# from pandas.testing import assert_frame_equal

# get directory of this file
path = Path(__file__).parent #os.getcwd() #path of working directory

# Create test files dir
testfiles_dir = os.path.join(path,"test_files")
if not os.path.exists(testfiles_dir):
    os.makedirs(testfiles_dir)


# <<<<<<< HEAD
#%% 
def test_modflow_run_phreatic_withgravelpack():
    ''' Phreatic scheme with gravelpack: modflow run.'''
    test_phrea = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
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
                                      halflife_suboxic= 530,
                                      halflife_anoxic= 2120,
                                      halflife_deeply_anoxic= 2120,
                                      partition_coefficient_water_organic_carbon= 6.43,
                                      dissociation_constant= 99,
                                      diameter_filterscreen = 0.2,
                                      point_input_concentration = 100.,
                                      discharge_point_contamination = 100.,#made up value
                                      top_clayseal = 17,
                                      compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),
                                      # substance = 'benzene',
                                      # halflife_suboxic=600,
                                      # partition_coefficient_water_organic_carbon = 3.3,
                                      ncols_near_well = 20,
                                      ncols_far_well = 80,
                                    )

    test_phrea.make_dictionary()
    # print(test_phrea.__dict__)
    modpath_phrea = MP.ModPathWell(test_phrea,
                            workspace = "test_ws",
                            modelname = "phreatic",
                            bound_left = "xmin",
                            bound_right = "xmax")
    # print(modpath_phrea.__dict__)
    # Run phreatic schematisation
    modpath_phrea.run_model(run_mfmodel = True,
                        run_mpmodel = False)

    print(modpath_phrea.material)
    parms = ["material","hk","vani","porosity","recharge"]
    for iParm in parms:
        np.save(os.path.join(testfiles_dir, iParm + "_phrea.npy"), getattr(modpath_phrea,iParm))

    assert modpath_phrea.success_mf
#%%
def test_modpath_run_phreatic_nogravelpack():
    ''' Phreatic scheme without gravelpack: modpath run.'''
    test_phrea = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                    computation_method = 'modpath',
                                    what_to_export='all',
                                    removal_function = 'pathogen',
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
                                      substance = 'norovirus',
                                      halflife_suboxic= 530,
                                      halflife_anoxic= 2120,
                                      halflife_deeply_anoxic= 2120,
                                      partition_coefficient_water_organic_carbon= 6.43,
                                      dissociation_constant= 99,
                                      # diameter_filterscreen = 0.2,
                                      point_input_concentration = 100.,
                                      discharge_point_contamination = 100.,#made up value
                                      top_clayseal = 17,
                                      compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),

                                      # substance = 'benzene',
                                      # halflife_suboxic=600,
                                      # partition_coefficient_water_organic_carbon = 3.3,
                                    
                                      ncols_near_well = 20,
                                      ncols_far_well = 80,
                                    )

    test_phrea.make_dictionary()
    phreatic_dict_1 = { 'simulation_parameters' : test_phrea.simulation_parameters,
                        'endpoint_id': test_phrea.endpoint_id,
                        'mesh_refinement': test_phrea.mesh_refinement,
                        'geo_parameters' : test_phrea.geo_parameters,
                        'ibound_parameters' : test_phrea.ibound_parameters,
                        'recharge_parameters' : test_phrea.recharge_parameters,
                        'well_parameters' : test_phrea.well_parameters,
                        'point_parameters' : test_phrea.point_parameters,
                        'substance_parameters' : test_phrea.substance_parameters,
                        'bas_parameters' : test_phrea.bas_parameters,
                        }
    # Remove/empty point_parameters
    test_phrea.point_parameters = {}

    # print(test_phrea.__dict__)
    modpath_phrea = MP.ModPathWell(test_phrea,
                            workspace = "test_phrea_nogp",
                            modelname = "phreatic",
                            bound_left = "xmin",
                            bound_right = "xmax")
    # print(modpath_phrea.__dict__)
    # Run phreatic schematisation
    modpath_phrea.run_model(run_mfmodel = True,
                        run_mpmodel = True)

    # Calculate advective microbial removal
    modpath_removal = ST.SubstanceTransport(modpath_phrea, substance = 'norovirus') #, df_particle, df_flowline)
 
    # Calculate advective microbial removal
    # Final concentration per endpoint_id
    C_final = {}
    for endpoint_id in modpath_phrea.schematisation_dict.get("endpoint_id"):
        df_particle, df_flowline, C_final[endpoint_id] = modpath_removal.calc_advective_microbial_removal(
                                            modpath_phrea.df_particle, modpath_phrea.df_flowline, 
                                            endpoint_id = endpoint_id,
                                            trackingdirection = modpath_phrea.trackingdirection,
                                            mu1 = 0.023, grainsize = 0.00025, alpha0 = 1.E-5, pH0 = 6.8, const_BM = 1.38e-23,
                                            temp_water = 11., rho_water = 999.703, species_diam = 2.33e-8,
                                            conc_start = 1., conc_gw = 0.)

    # df_particle file name 
    particle_fname = os.path.join(modpath_phrea.dstroot,modpath_phrea.schematisation_type + "_df_particle_microbial_removal.csv")
    # Save df_particle 
    df_particle.to_csv(particle_fname)
    
    # df_flowline file name
    flowline_fname = os.path.join(modpath_phrea.dstroot,modpath_phrea.schematisation_type + "_df_flowline_microbial_removal.csv")
    # Save df_flowline
    df_flowline.to_csv(flowline_fname)
    
    assert modpath_phrea.success_mp

def test_modpath_run_horizontal_flow_points():
    ''' Horizontal flow test in target_aquifer: modpath run.'''
    # well discharge
    well_discharge = -1000.
    # distance to boundary
    distance_boundary = 550.
    # Center depth of target aquifer
    z_point = 0 - 0.1 - 10.
    x_point = distance_boundary - 1.
    # Phreatic scheme without gravelpack: modpath run.
    test_conf_hor = AW.HydroChemicalSchematisation(schematisation_type='semiconfined',
                                computation_method = 'modpath',
                                what_to_export='all',
                                removal_function = 'pathogen',
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
                                thickness_target_aquifer=20.0,
                                hor_permeability_target_aquifer=10.0,
                                hor_permeability_shallow_aquifer = 1.,
                                vertical_anisotropy_shallow_aquifer = 1000.,
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
                                temperature=12.,
                                solid_density_vadose_zone= 2.650, 
                                solid_density_shallow_aquifer= 2.650, 
                                solid_density_target_aquifer= 2.650, 
                                diameter_borehole = 0.2,
                                substance = 'norovirus',
                                halflife_suboxic= 530,
                                halflife_anoxic= 2120,
                                halflife_deeply_anoxic= 2120,
                                partition_coefficient_water_organic_carbon= 6.43,
                                dissociation_constant= 99,
                                diameter_filterscreen = 0.2,
                                top_clayseal = 0,
                                compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),
                                # substance = 'benzene',
                                # halflife_suboxic=600,
                                # partition_coefficient_water_organic_carbon = 3.3,
                                model_radius = distance_boundary,
                                ncols_near_well = 20,
                                ncols_far_well = 80,
                                # Point contamination
                                point_input_concentration = 1.,
                                discharge_point_contamination = abs(well_discharge),
                                distance_point_contamination_from_well = x_point,
                                depth_point_contamination = z_point,
                                )

    test_conf_hor.make_dictionary()
    ### Adjust ibound_parameters to add horizontal flow ###
    
    # Confined top boundary ; no recharge_parameters
    test_conf_hor.ibound_parameters.pop("top_boundary1")
    # Add outer boundary for horizontal flow test
    test_conf_hor.ibound_parameters["outer_boundary_target_aquifer"] = {
                            'top': test_conf_hor.bottom_shallow_aquifer,
                            'bot': test_conf_hor.bottom_target_aquifer,
                            'xmin': test_conf_hor.model_radius - 1.,
                            'xmax': test_conf_hor.model_radius,
                            'ibound': -1
                            }
    # Remove/empty diffuse_parameters
    test_conf_hor.diffuse_parameters = {}
    

    # Add point parameters                       
    startpoint_id = ["outer_boundary_target_aquifer"]
    

    # microbial removal properties
    substance_name = 'norovirus'
    alpha0 = {"suboxic": 1.e-3, "anoxic": 1.e-5, "deeply_anoxic": 1.e-5},
    reference_pH = {"suboxic": 6.6, "anoxic": 6.8, "deeply_anoxic": 6.8},
    species_diam =  2.33e-8,
    mu1 = {"suboxic": 0.149,"anoxic": 0.023,"deeply_anoxic": 0.023}

    substance_parameters = {substance_name: 
                    {"substance_name": substance_name,
                        "alpha0": alpha0,
                        "reference_pH": reference_pH,
                        "species_diam": species_diam,
                        "mu1": mu1
                    }
                }
    test_conf_hor.substance_parameters["norovirus"] = substance_parameters["norovirus"]


    phreatic_dict_2 = { 'simulation_parameters' : test_conf_hor.simulation_parameters,
                    'endpoint_id': test_conf_hor.endpoint_id,
                    'mesh_refinement': test_conf_hor.mesh_refinement,
                    'geo_parameters' : test_conf_hor.geo_parameters,
                    'ibound_parameters' : test_conf_hor.ibound_parameters,
                    'recharge_parameters' : test_conf_hor.recharge_parameters,
                    'well_parameters' : test_conf_hor.well_parameters,
                    'point_parameters' : test_conf_hor.point_parameters,
                    'substance_parameters' : test_conf_hor.substance_parameters,
                    'bas_parameters' : test_conf_hor.bas_parameters,
                    }

    # ModPath well
    modpath_hor = MP.ModPathWell(test_conf_hor, #test_phrea,
                            workspace = "test_ws_conf_hor_pnts",
                            modelname = "confined_hor",
                            bound_left = "xmin",
                            bound_right = "xmax")
    # print(modpath_phrea.__dict__)
    # Run phreatic schematisation
    modpath_hor.run_model(run_mfmodel = True,
                        run_mpmodel = True)


    # Calculate advective microbial removal
    modpath_removal = ST.SubstanceTransport(modpath_hor, substance = 'norovirus') #, df_particle, df_flowline)
 
    # Calculate advective microbial removal
    # Final concentration per endpoint_id
    C_final = {}
    for endpoint_id in modpath_hor.schematisation_dict.get("endpoint_id"):
        df_particle, df_flowline, C_final[endpoint_id] = modpath_removal.calc_advective_microbial_removal(
                                            modpath_hor.df_particle, modpath_hor.df_flowline, 
                                            endpoint_id = endpoint_id,
                                            trackingdirection = modpath_hor.trackingdirection,
                                            mu1 = 0.023, grainsize = 0.00025, alpha0 = 1.E-5, pH0 = 6.8, const_BM = 1.38e-23,
                                            temp_water = 11., rho_water = 999.703, species_diam = 2.33e-8,
                                            conc_start = 1., conc_gw = 0.)

    # df_particle file name 
    particle_fname = os.path.join(modpath_hor.dstroot,modpath_hor.schematisation_type + "_df_particle_microbial_removal.csv")
    # Save df_particle 
    df_particle.to_csv(particle_fname)
    
    # df_flowline file name
    flowline_fname = os.path.join(modpath_hor.dstroot,modpath_hor.schematisation_type + "_df_flowline_microbial_removal.csv")
    # Save df_flowline
    df_flowline.to_csv(flowline_fname)

    assert modpath_hor.success_mp

def test_modpath_run_horizontal_flow_diffuse():
    ''' Horizontal flow test in target_aquifer: modpath run.'''
    # well discharge
    well_discharge = -1000.
    # distance to boundary
    distance_boundary = 550.
    # Center depth of target aquifer
    z_point = 0 - 0.1 - 10.
    x_point = distance_boundary - 1.
    # Phreatic scheme without gravelpack: modpath run.
    test_conf_hor = AW.HydroChemicalSchematisation(schematisation_type='semiconfined',
                                computation_method = 'modpath',
                                what_to_export='all',
                                removal_function = 'pathogen',
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
                                thickness_target_aquifer=20.0,
                                hor_permeability_target_aquifer=10.0,
                                hor_permeability_shallow_aquifer = 1.,
                                vertical_anisotropy_shallow_aquifer = 1000.,
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
                                temperature=12.,
                                solid_density_vadose_zone= 2.650, 
                                solid_density_shallow_aquifer= 2.650, 
                                solid_density_target_aquifer= 2.650, 
                                diameter_borehole = 0.2,
                                substance = 'norovirus',
                                halflife_suboxic= 530,
                                halflife_anoxic= 2120,
                                halflife_deeply_anoxic= 2120,
                                partition_coefficient_water_organic_carbon= 6.43,
                                dissociation_constant= 99,
                                diameter_filterscreen = 0.2,
                                top_clayseal = 0,
                                compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),
                                # substance = 'benzene',
                                # halflife_suboxic=600,
                                # partition_coefficient_water_organic_carbon = 3.3,
                                model_radius = distance_boundary,
                                ncols_near_well = 20,
                                ncols_far_well = 80,
                                # Point contamination
                                point_input_concentration = 1.,
                                discharge_point_contamination = abs(well_discharge),
                                distance_point_contamination_from_well = x_point,
                                depth_point_contamination = z_point,
                                )

    test_conf_hor.make_dictionary()
    ### Adjust ibound_parameters to add horizontal flow ###
    
    # Confined top boundary ; no recharge_parameters
    test_conf_hor.ibound_parameters.pop("top_boundary1")
    # Add outer boundary for horizontal flow test
    test_conf_hor.ibound_parameters["outer_boundary_target_aquifer"] = {
                            'top': test_conf_hor.bottom_shallow_aquifer,
                            'bot': test_conf_hor.bottom_target_aquifer,
                            'xmin': test_conf_hor.model_radius - 1.,
                            'xmax': test_conf_hor.model_radius,
                            'ibound': -1
                            }
    # diffuse_parameters
    test_conf_hor.diffuse_parameters = {'source1': {
                'substance_name': "norovirus",
                'xmin': test_conf_hor.model_radius - 1.,
                'xmax': test_conf_hor.model_radius,
                'zmax': test_conf_hor.bottom_shallow_aquifer,
                'zmin': test_conf_hor.bottom_target_aquifer,
                'input_concentration': test_conf_hor.diffuse_input_concentration
                }
            # 'source2' :{}> surface water (BAR & RBF) #@MartinvdS come back to this when we start this module
            }
    # Remove/empty point_parameters
    test_conf_hor.point_parameters = {}

    # Add point parameters                       
    startpoint_id = ["outer_boundary_target_aquifer"]
    

    # microbial removal properties
    substance_name = 'norovirus'
    alpha0 = {"suboxic": 1.e-3, "anoxic": 1.e-5, "deeply_anoxic": 1.e-5},
    reference_pH = {"suboxic": 6.6, "anoxic": 6.8, "deeply_anoxic": 6.8},
    species_diam =  2.33e-8,
    mu1 = {"suboxic": 0.149,"anoxic": 0.023,"deeply_anoxic": 0.023}

    substance_parameters = {substance_name: 
                    {"substance_name": substance_name,
                        "alpha0": alpha0,
                        "reference_pH": reference_pH,
                        "species_diam": species_diam,
                        "mu1": mu1
                    }
                }
    test_conf_hor.substance_parameters["norovirus"] = substance_parameters["norovirus"]


    phreatic_dict_2 = { 'simulation_parameters' : test_conf_hor.simulation_parameters,
                    'endpoint_id': test_conf_hor.endpoint_id,
                    'mesh_refinement': test_conf_hor.mesh_refinement,
                    'geo_parameters' : test_conf_hor.geo_parameters,
                    'ibound_parameters' : test_conf_hor.ibound_parameters,
                    'recharge_parameters' : test_conf_hor.recharge_parameters,
                    'well_parameters' : test_conf_hor.well_parameters,
                    'point_parameters' : test_conf_hor.point_parameters,
                    'substance_parameters' : test_conf_hor.substance_parameters,
                    'bas_parameters' : test_conf_hor.bas_parameters,
                    }

    # ModPath well
    modpath_hor = MP.ModPathWell(test_conf_hor, #test_phrea,
                            workspace = "test_ws_conf_hor_diffuse",
                            modelname = "confined_hor",
                            bound_left = "xmin",
                            bound_right = "xmax")
    # print(modpath_phrea.__dict__)
    # Run phreatic schematisation
    modpath_hor.run_model(run_mfmodel = True,
                        run_mpmodel = True)


    # Calculate advective microbial removal
    modpath_removal = ST.SubstanceTransport(modpath_hor, substance = 'norovirus') #, df_particle, df_flowline)
 
    # Calculate advective microbial removal
    # Final concentration per endpoint_id
    C_final = {}
    for endpoint_id in modpath_hor.schematisation_dict.get("endpoint_id"):
        df_particle, df_flowline, C_final[endpoint_id] = modpath_removal.calc_advective_microbial_removal(
                                            modpath_hor.df_particle, modpath_hor.df_flowline, 
                                            endpoint_id = endpoint_id,
                                            trackingdirection = modpath_hor.trackingdirection,
                                            mu1 = 0.023, grainsize = 0.00025, alpha0 = 1.E-5, pH0 = 6.8, const_BM = 1.38e-23,
                                            temp_water = 11., rho_water = 999.703, species_diam = 2.33e-8,
                                            conc_start = 1., conc_gw = 0.)

    # df_particle file name 
    particle_fname = os.path.join(modpath_hor.dstroot,modpath_hor.schematisation_type + "_df_particle_microbial_removal.csv")
    # Save df_particle 
    df_particle.to_csv(particle_fname)
    
    # df_flowline file name
    flowline_fname = os.path.join(modpath_hor.dstroot,modpath_hor.schematisation_type + "_df_flowline_microbial_removal.csv")
    # Save df_flowline
    df_flowline.to_csv(flowline_fname)

    assert modpath_hor.success_mp

def test_modpath_run_phreatic_withgravelpack_traveltimes():


    ''' Phreatic scheme with gravelpack: modpath run.'''
    test_phrea_gp = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                computation_method = 'modpath',
                                what_to_export='all',
                                removal_function = 'pathogen',
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
                                substance = 'norovirus',
                                halflife_suboxic= 530,
                                halflife_anoxic= 2120,
                                halflife_deeply_anoxic= 2120,
                                partition_coefficient_water_organic_carbon= 6.43,
                                dissociation_constant= 99,
                                diameter_filterscreen = 0.2,
                                point_input_concentration = 100.,
                                discharge_point_contamination = 100.,#made up value
                                top_clayseal = 17,
                                compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),
                                # substance = 'benzene',
                                # halflife_suboxic=600,
                                # partition_coefficient_water_organic_carbon = 3.3,
                                ncols_near_well = 20,
                                ncols_far_well = 80,
                            )

    test_phrea_gp.make_dictionary()

    phreatic_dict_2 = { 'simulation_parameters' : test_phrea_gp.simulation_parameters,
                        'endpoint_id': test_phrea_gp.endpoint_id,
                        'mesh_refinement': test_phrea_gp.mesh_refinement,
                        'geo_parameters' : test_phrea_gp.geo_parameters,
                        'ibound_parameters' : test_phrea_gp.ibound_parameters,
                        'recharge_parameters' : test_phrea_gp.recharge_parameters,
                        'well_parameters' : test_phrea_gp.well_parameters,
                        'point_parameters' : test_phrea_gp.point_parameters,
                        'substance_parameters' : test_phrea_gp.substance_parameters,
                        'bas_parameters' : test_phrea_gp.bas_parameters,
                        }
    # Remove/empty point_parameters
    test_phrea_gp.point_parameters = {}

    # print(test_phrea.__dict__)
    modpath_phrea = MP.ModPathWell(test_phrea_gp, #test_phrea,
                            workspace = "test_ws_phrea_gp",
                            modelname = "phreatic",
                            bound_left = "xmin",
                            bound_right = "xmax")
    # print(modpath_phrea.__dict__)
    # Run phreatic schematisation
    modpath_phrea.run_model(run_mfmodel = True,
                        run_mpmodel = True)

    
    # Create travel time plots
    fpath_scatter_times_log = os.path.join(modpath_phrea.dstroot,"log_travel_times_test.png")
    fpath_scatter_times = os.path.join(modpath_phrea.dstroot,"travel_times_test.png")
    # df particle
    df_particle = modpath_phrea.df_particle
    # time limits
    tmin, tmax = 0.1, 10000.
    # xcoord bounds
    xmin, xmax = 0., 50.

    # Create travel time plots (lognormal)
    modpath_phrea.plot_pathtimes(df_particle = df_particle, 
            vmin = tmin,vmax = tmax,
            fpathfig = fpath_scatter_times_log, figtext = None,x_text = 0,
            y_text = 0, lognorm = True, xmin = xmin, xmax = xmax,
            line_dist = 1, dpi = 192, trackingdirection = "forward",
            cmap = 'viridis_r')

    # Create travel time plots (linear)
    modpath_phrea.plot_pathtimes(df_particle = df_particle, 
            vmin = 0.,vmax = tmax,
            fpathfig = fpath_scatter_times, figtext = None,x_text = 0,
            y_text = 0, lognorm = False, xmin = xmin, xmax = xmax,
            line_dist = 1, dpi = 192, trackingdirection = "forward",
            cmap = 'viridis_r')

    assert modpath_phrea.success_mp

#=======
#%%

def test_modpath_run_semiconfined_nogravelpack_traveltimes():


    ''' Phreatic scheme with gravelpack: modpath run.'''
    test_semiconf = AW.HydroChemicalSchematisation(schematisation_type='semiconfined',
                                    computation_method = 'modpath',
                                    what_to_export='all',
                                    removal_function = 'pathogen',
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
                                      substance = 'norovirus',
                                      halflife_suboxic= 530,
                                      halflife_anoxic= 2120,
                                      halflife_deeply_anoxic= 2120,
                                      partition_coefficient_water_organic_carbon= 6.43,
                                      dissociation_constant= 99,
                                      # diameter_filterscreen = 0.2,
                                      point_input_concentration = 100.,
                                      discharge_point_contamination = 100.,#made up value
                                      top_clayseal = 17,
                                      compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),
                                      
                                      # substance = 'benzene',
                                      # halflife_suboxic=600,
                                      # partition_coefficient_water_organic_carbon = 3.3,
                                      ncols_near_well = 20,
                                      ncols_far_well = 20,
                                      nlayers_target_aquifer = 20,
                                    #   # Pathogen removal parameters
                                    #   mu1 = 0.023,  # nog uitschrijven
                                    #   alpha0 = 

                                    )

    test_semiconf.make_dictionary()

    # Test diffuse_parameters (check line 796 in ModPath_functions)
    diffuse_parameters = test_semiconf.recharge_parameters
    # diffuse_parameters['xmin'] = 15.
    # diffuse_parameters['xmax'] = 16.
    semiconf_dict_1 = { 'simulation_parameters' : test_semiconf.simulation_parameters,
                        'endpoint_id': test_semiconf.endpoint_id,
                        'mesh_refinement': test_semiconf.mesh_refinement,
                        'geo_parameters' : test_semiconf.geo_parameters,
                        'ibound_parameters' : test_semiconf.ibound_parameters,
                        'recharge_parameters' : test_semiconf.recharge_parameters,
                        'well_parameters' : test_semiconf.well_parameters,
                        'diffuse_parameters': test_semiconf.recharge_parameters,
                        'point_parameters' : test_semiconf.point_parameters,
                        'substance_parameters' : test_semiconf.substance_parameters,
                        'bas_parameters' : test_semiconf.bas_parameters,
                        }
    # Remove/empty point_parameters
    test_semiconf.point_parameters = {}

    modpath_semiconf = MP.ModPathWell(test_semiconf, # semiconf_dict_1, #test_phrea,
                            workspace = "test1_ws_semiconf_nogp",
                            modelname = "semi_conf_nogp",
                            bound_left = "xmin",
                            bound_right = "xmax")

    # Run phreatic schematisation
    modpath_semiconf.run_model(run_mfmodel = True,
                        run_mpmodel = True)

    # Calculate advective microbial removal
    modpath_removal = ST.SubstanceTransport(modpath_semiconf, substance = 'norovirus') #, df_particle, df_flowline)
    # modpath_removal.compute_omp_removal()
    # Final concentration per endpoint_id
    C_final = {}
    for endpoint_id in modpath_semiconf.schematisation_dict.get("endpoint_id"):
        df_particle, df_flowline, C_final[endpoint_id] = modpath_removal.calc_advective_microbial_removal(
                                            modpath_semiconf.df_particle, modpath_semiconf.df_flowline, 
                                            endpoint_id = endpoint_id,
                                            trackingdirection = modpath_semiconf.trackingdirection,
                                            mu1 = 0.023, grainsize = 0.00025, alpha0 = 1.E-5, pH0 = 6.8, const_BM = 1.38e-23,
                                            temp_water = 11., rho_water = 999.703, species_diam = 2.33e-8,
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

#%%

def test_diffuse_modpath_run_semiconfined_nogravelpack_traveltimes():

    ''' Phreatic scheme without gravelpack: modpath run.'''
    test_semiconf = AW.HydroChemicalSchematisation(schematisation_type='semiconfined',
                                    computation_method = 'modpath',
                                    what_to_export='all',
                                    removal_function = 'pathogen',
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
                                      halflife_suboxic= 530,
                                      halflife_anoxic= 2120,
                                      halflife_deeply_anoxic= 2120,
                                      partition_coefficient_water_organic_carbon= 6.43,
                                      dissociation_constant= 99,
                                      # diameter_filterscreen = 0.2,
                                      point_input_concentration = 100.,
                                      discharge_point_contamination = 100.,#made up value
                                      top_clayseal = 17,
                                      compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),
                                      
                                      # substance = 'benzene',
                                      # halflife_suboxic=600,
                                      # partition_coefficient_water_organic_carbon = 3.3,
                                      ncols_near_well = 20,
                                      ncols_far_well = 20,
                                      nlayers_target_aquifer = 20,
                                    )

    test_semiconf.make_dictionary()

    semiconf_dict_1 = { 'simulation_parameters' : test_semiconf.simulation_parameters,
                        'endpoint_id': test_semiconf.endpoint_id,
                        'mesh_refinement': test_semiconf.mesh_refinement,
                        'geo_parameters' : test_semiconf.geo_parameters,
                        'ibound_parameters' : test_semiconf.ibound_parameters,
                        'recharge_parameters' : test_semiconf.recharge_parameters,
                        'well_parameters' : test_semiconf.well_parameters,
                        'point_parameters' : test_semiconf.point_parameters,
                        'substance_parameters' : test_semiconf.substance_parameters,
                        'bas_parameters' : test_semiconf.bas_parameters,
                        }
        # Remove/empty point_parameters
    test_semiconf.point_parameters = {}

    modpath_semiconf = MP.ModPathWell(test_semiconf, #test_phrea,
                            workspace = "test1_ws_semiconf",
                            modelname = "semi_conf_no_gp",
                            bound_left = "xmin",
                            bound_right = "xmax")

    # Run phreatic schematisation
    modpath_semiconf.run_model(run_mfmodel = True,
                        run_mpmodel = True)
    
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

#%%

def test_phreatic_scheme_withgravelpack_dictinput():
    ''' Check writing and reading of dictionary files: phreatic scheme.'''
    
    phreatic_scheme= AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                    computation_method = 'modpath',
                                    what_to_export='all',
                                    removal_function = 'pathogen',
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
                                      halflife_suboxic= 530,
                                      halflife_anoxic= 2120,
                                      halflife_deeply_anoxic= 2120,
                                      partition_coefficient_water_organic_carbon= 6.43,
                                      dissociation_constant= 99,
                                      diameter_filterscreen = 0.2,
                                      point_input_concentration = 100.,
                                      discharge_point_contamination = 100.,#made up value
                                      top_clayseal = 17,
                                      compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),
                                      # substance = 'benzene',
                                      # halflife_suboxic=600,
                                      # partition_coefficient_water_organic_carbon = 3.3,
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
        'substance_parameters' : phreatic_scheme.substance_parameters,
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

def test_travel_time_distribution_phreatic():
    output_phreatic = pd.read_csv(path / 'phreatic_test.csv')
    output_phreatic = output_phreatic.round(7) #round to 7 digits (or any digit), keep same as for the output for the model to compare

    test_ = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                        computation_method= 'analytical',
                                        what_to_export='all',
                                        removal_function = 'pathogen',
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

    well1 = AW.AnalyticalWell(test_)
    well1.phreatic()
    output = well1.df_output
    output = output[["total_travel_time", "travel_time_unsaturated",
                     "travel_time_shallow_aquifer", "travel_time_target_aquifer",
                     "radial_distance", ]]
    output = output.round(7)

    try:
        assert_frame_equal(output, output_phreatic,check_dtype=False)

    except AssertionError:
        print("Assertion Exception Raised - TTD test")
    else:
        print("Success, no error in TTD!")

# # %%
# if __name__ == "__main__":
#     test_modpath_run_phreatic_withgravelpack()
