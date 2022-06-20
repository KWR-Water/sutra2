
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
def test_modflow_run_phreatic_withgravelpack(organism_name = "MS2"):
    ''' Phreatic scheme with gravelpack: modflow run.'''
    test_phrea = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                    computation_method = 'modpath',
                                    what_to_export='all',
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
                                      name = organism_name,
                                      diameter_filterscreen = 0.2,
                                      point_input_concentration = 100.,
                                      discharge_point_contamination = 100.,#made up value
                                      top_clayseal = 17,
                                      compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),
                                      ncols_near_well = 20,
                                      ncols_far_well = 80,
                                    )

    test_phrea.make_dictionary()
    # print(test_phrea.__dict__)
    modpath_phrea = MP.ModPathWell(test_phrea,
                            workspace = "test1_phrea_mf",
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
def test_modpath_run_phreatic_nogravelpack(organism_name = "MS2"):
    ''' Phreatic scheme without gravelpack: modpath run.'''
    test_phrea = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
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
                                      temperature=11.,
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
                                      ncols_far_well = 80,
                                    )

    test_phrea.make_dictionary()
    # Remove/empty point_parameters
    test_phrea.point_parameters = {}

    # print(test_phrea.__dict__)
    modpath_phrea = MP.ModPathWell(test_phrea,
                            workspace = "test2_phrea_nogp",
                            modelname = "phreatic",
                            bound_left = "xmin",
                            bound_right = "xmax")
    # print(modpath_phrea.__dict__)
    # Run phreatic schematisation
    modpath_phrea.run_model(run_mfmodel = True,
                        run_mpmodel = True)

    # Calculate advective microbial removal
    modpath_removal = ST.SubstanceTransport(modpath_phrea,
                                            organism = organism_name)

    # Calculate advective microbial removal
    # Final concentration per endpoint_id
    C_final = {}
    for endpoint_id in modpath_phrea.schematisation_dict.get("endpoint_id"):
        df_particle, df_flowline, C_final[endpoint_id] = modpath_removal.calc_advective_microbial_removal(
                                            modpath_phrea.df_particle, modpath_phrea.df_flowline,
                                            endpoint_id = endpoint_id,
                                            trackingdirection = modpath_phrea.tracking_direction,
                                            mu1 = 0.023, grainsize = 0.00025, alpha0 = 1.E-5, reference_pH = 6.8, const_BM = 1.38e-23,
                                            temp_water = 11., rho_water = 999.703, organism_diam = 2.33e-8,
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

def test_modpath_run_horizontal_flow_points(organism_name = "MS2"):
    ''' Horizontal flow test in target_aquifer: modpath run.'''

    # microbial removal properties
    # organism_name = 'MS2'
    alpha0 = {"suboxic": 1.e-3, "anoxic": 1.e-5, "deeply_anoxic": 1.e-5}
    reference_pH = {"suboxic": 6.6, "anoxic": 6.8, "deeply_anoxic": 6.8}
    organism_diam =  2.33e-8
    mu1 = {"suboxic": 0.039,"anoxic": 0.023,"deeply_anoxic": 0.023}

    removal_parameters = {organism_name:
                    {"organism_name": organism_name,
                        "alpha0": alpha0,
                        "reference_pH": reference_pH,
                        "organism_diam": organism_diam,
                        "mu1": mu1
                    }
                }
    # Removal parameters organism
    rem_parms = removal_parameters[organism_name]

    # dataframes
    df_particle, df_flowline = {}, {}

    # Distances to boundary
    dist_boundary = list(range(25,100,5)) + list(range(100,560,10)) # [260] #  +
    # [10,20, 50, 100, 150, 200,250,300,350,400,450,500,550]

    df_conc = pd.DataFrame(index=dist_boundary, columns = ["Final_concentration"])
    df_conc.index.name = "Boundary_distance"

    for iDist in dist_boundary:


        # well discharge
        well_discharge = -1000.
        # distance to boundary
        distance_boundary = float(iDist)
        # Center depth of target aquifer
        z_point = 0 - 0.1 - 10.
        x_point = distance_boundary - 0.5
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
                                    name = organism_name,
                                    diameter_filterscreen = 0.2,
                                    top_clayseal = 0,
                                    compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),
                                    model_radius = distance_boundary,
                                    ncols_near_well = 20,
                                    ncols_far_well = int(distance_boundary/10),
                                    # Point contamination
                                    point_input_concentration = 1.,
                                    discharge_point_contamination = abs(well_discharge),
                                    distance_point_contamination_from_well = x_point,
                                    depth_point_contamination = z_point,
                                    )

        test_conf_hor.make_dictionary()


        # Remove "vadose" layer from geo_parameters
        test_conf_hor.geo_parameters.pop("vadose")

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
        # Remove/empty concentration_boundary_parameters
        test_conf_hor.concentration_boundary_parameters = {}

        # ModPath well
        modpath_hor = MP.ModPathWell(test_conf_hor, #test_phrea,
                                workspace = "test3_conf_hor_pnts",
                                modelname = "confined_hor",
                                bound_left = "xmin",
                                bound_right = "xmax")
        # print(modpath_phrea.__dict__)

        # Run phreatic schematisation
        modpath_hor.run_model(run_mfmodel = True,
                            run_mpmodel = True)


        # Calculate advective microbial removal
        modpath_removal = ST.SubstanceTransport(modpath_hor,
                                organism = organism_name,
                                alpha0_suboxic = rem_parms["alpha0"]["suboxic"],
                                alpha0_anoxic = rem_parms["alpha0"]["anoxic"],
                                alpha0_deeply_anoxic =rem_parms["alpha0"]["deeply_anoxic"],
                                reference_pH_suboxic =rem_parms["reference_pH"]["suboxic"],
                                reference_pH_anoxic =rem_parms["reference_pH"]["anoxic"],
                                reference_pH_deeply_anoxic =rem_parms["reference_pH"]["deeply_anoxic"],
                                mu1_suboxic = rem_parms["mu1"]["suboxic"],
                                mu1_anoxic = rem_parms["mu1"]["anoxic"],
                                mu1_deeply_anoxic = rem_parms["mu1"]["deeply_anoxic"],
                                organism_diam = rem_parms["organism_diam"]
                                )

        # Calculate advective microbial removal
        # Final concentration per endpoint_id
        C_final = {}
        for endpoint_id in modpath_hor.schematisation_dict.get("endpoint_id"):
            df_particle[iDist], df_flowline[iDist], C_final[endpoint_id] = modpath_removal.calc_advective_microbial_removal(
                                                modpath_hor.df_particle, modpath_hor.df_flowline,
                                                endpoint_id = endpoint_id,
                                                trackingdirection = modpath_hor.tracking_direction,
                                                conc_start = 1., conc_gw = 0.)

            # print("Final concentration at " + str(iDist) + " m is: " + str(round(C_final[endpoint_id],4)))
            # Add final concentration to summary dataframe
            df_conc.loc[iDist,"Final_concentration " + endpoint_id] = C_final[endpoint_id]

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

def test_modpath_run_horizontal_flow_diffuse(organism_name = "MS2"):

    ''' Horizontal flow test in target_aquifer: modpath run.'''
    # well discharge
    well_discharge = -1000.
    # distance to boundary
    distance_boundary = 50.
    # Center depth of target aquifer
    z_point = 0 - 0.1 - 10.
    x_point = distance_boundary - 0.5
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
                                name = organism_name,
                                diameter_filterscreen = 0.2,
                                top_clayseal = 0,
                                compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),
                                model_radius = distance_boundary,
                                ncols_near_well = 10,
                                ncols_far_well = int(distance_boundary/5),
                                # Point contamination
                                point_input_concentration = 1.,
                                discharge_point_contamination = abs(well_discharge),
                                distance_point_contamination_from_well = x_point,
                                depth_point_contamination = z_point,
                                )

    test_conf_hor.make_dictionary()


    # Remove "vadose" layer from geo_parameters
    test_conf_hor.geo_parameters.pop("vadose")

    ### Adjust ibound_parameters to add horizontal flow ###

    # Confined top boundary ; no recharge_parameters
    test_conf_hor.ibound_parameters.pop("top_boundary1")

    # Width (delr) outer boundary cells
    delr_outer = 1.

    # Add outer boundary for horizontal flow test
    test_conf_hor.ibound_parameters["outer_boundary_target_aquifer"] = {
                            'top': test_conf_hor.bottom_shallow_aquifer,
                            'bot': test_conf_hor.bottom_target_aquifer,
                            'xmin': test_conf_hor.model_radius - delr_outer,
                            'xmax': test_conf_hor.model_radius,
                            'ibound': -1
                            }
    # concentration_boundary_parameters
    test_conf_hor.concentration_boundary_parameters = {'source1': {
                'organism_name': organism_name,
                'xmin': test_conf_hor.model_radius - delr_outer * 0.5,
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


    # ModPath well
    modpath_hor = MP.ModPathWell(test_conf_hor, #test_phrea,
                            workspace = "test4_conf_hor_diffuse",
                            modelname = "confined_hor",
                            bound_left = "xmin",
                            bound_right = "xmax")
    # print(modpath_phrea.__dict__)
    # Run phreatic schematisation
    modpath_hor.run_model(run_mfmodel = True,
                        run_mpmodel = True)

    # microbial removal properties
    # organism_name = 'MS2'
    alpha0 = {"suboxic": 1.e-3, "anoxic": 1.e-5, "deeply_anoxic": 1.e-5}
    reference_pH = {"suboxic": 6.6, "anoxic": 6.8, "deeply_anoxic": 6.8}
    organism_diam =  2.33e-8
    mu1 = {"suboxic": 0.149,"anoxic": 0.023,"deeply_anoxic": 0.023}

    removal_parameters = {organism_name:
                    {"organism_name": organism_name,
                        "alpha0": alpha0,
                        "reference_pH": reference_pH,
                        "organism_diam": organism_diam,
                        "mu1": mu1
                    }
                }
    # Removal parameters organism
    rem_parms = removal_parameters[organism_name]

    # Calculate advective microbial removal
    modpath_removal = ST.SubstanceTransport(modpath_hor,
                            organism = organism_name,
                            alpha0_suboxic = rem_parms["alpha0"]["suboxic"],
                            alpha0_anoxic = rem_parms["alpha0"]["anoxic"],
                            alpha0_deeply_anoxic =rem_parms["alpha0"]["deeply_anoxic"],
                            reference_pH_suboxic =rem_parms["reference_pH"]["suboxic"],
                            reference_pH_anoxic =rem_parms["reference_pH"]["anoxic"],
                            reference_pH_deeply_anoxic =rem_parms["reference_pH"]["deeply_anoxic"],
                            mu1_suboxic = rem_parms["mu1"]["suboxic"],
                            mu1_anoxic = rem_parms["mu1"]["anoxic"],
                            mu1_deeply_anoxic = rem_parms["mu1"]["deeply_anoxic"],
                            organism_diam = rem_parms["organism_diam"]
                            )

    # Calculate advective microbial removal
    # Final concentration per endpoint_id
    C_final = {}
    for endpoint_id in modpath_hor.schematisation_dict.get("endpoint_id"):
        df_particle, df_flowline, C_final[endpoint_id] = modpath_removal.calc_advective_microbial_removal(
                                            modpath_hor.df_particle, modpath_hor.df_flowline,
                                            endpoint_id = endpoint_id,
                                            trackingdirection = modpath_hor.tracking_direction,
                                            conc_start = 1., conc_gw = 0.)

    # df_particle file name
    particle_fname = os.path.join(modpath_hor.dstroot,modpath_hor.schematisation_type + "_df_particle_microbial_removal.csv")
    # Save df_particle
    df_particle.to_csv(particle_fname)

    # df_flowline file name
    flowline_fname = os.path.join(modpath_hor.dstroot,modpath_hor.schematisation_type + "_df_flowline_microbial_removal.csv")
    # Save df_flowline
    df_flowline.to_csv(flowline_fname)


    # Create travel time plots
    fpath_scatter_times_log = os.path.join(modpath_hor.dstroot,"log_travel_times_" + endpoint_id + ".png")
    fpath_scatter_times = os.path.join(modpath_hor.dstroot,"travel_times_" + endpoint_id + ".png")
    # # df particle
    # df_particle = modpath_removal.df_particle
    # time limits
    tmin, tmax = 0.1, 10000.
    # xcoord bounds
    xmin, xmax = 0., min(50., x_point)

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
    fpath_scatter_removal_log = os.path.join(modpath_hor.dstroot,"log_removal_" + endpoint_id + ".png")

    # relative conc limits
    cmin, cmax = 1.e-21, 1.
    # xcoord bounds
    xmin, xmax = 0., min(50., x_point)

    # Create travel time plots (lognormal)
    modpath_removal.plot_logremoval(df_particle=df_particle,
            df_flowline=df_flowline,
            vmin = cmin,vmax = cmax,
            fpathfig = fpath_scatter_removal_log, figtext = None,x_text = 0,
            y_text = 0, lognorm = True, xmin = xmin, xmax = xmax,
            line_dist = 1, dpi = 192, trackingdirection = "forward",
            cmap = 'viridis_r')

    assert modpath_hor.success_mp

def test_modpath_run_phreatic_withgravelpack_traveltimes(organism_name = "MS2"):


    ''' Phreatic scheme with gravelpack: modpath run.'''
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
                                temperature=11.,
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
    modpath_phrea = MP.ModPathWell(test_phrea_gp, #test_phrea,
                            workspace = "test5_phrea_gp",
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

#%%
def test_modpath_run_phreatic_withgravelpack_removal(organism_name = "MS2"):


    ''' Phreatic scheme with gravelpack: modpath run.'''
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
                                temperature=11.,
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
    modpath_phrea = MP.ModPathWell(test_phrea_gp, #test_phrea,
                            workspace = "test6_phrea_gp_removal",
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
    reference_pH = {"suboxic": 6.6, "anoxic": 6.8, "deeply_anoxic": 6.8}
    organism_diam =  2.33e-8
    mu1 = {"suboxic": 0.149,"anoxic": 0.023,"deeply_anoxic": 0.023}

    removal_parameters = {organism_name:
                    {"organism_name": organism_name,
                        "alpha0": alpha0,
                        "reference_pH": reference_pH,
                        "organism_diam": organism_diam,
                        "mu1": mu1
                    }
                }
    # Removal parameters organism
    rem_parms = removal_parameters[organism_name]

    # Calculate advective microbial removal
    modpath_removal = ST.SubstanceTransport(modpath_phrea,
                            organism = organism_name,
                            alpha0_suboxic = rem_parms["alpha0"]["suboxic"],
                            alpha0_anoxic = rem_parms["alpha0"]["anoxic"],
                            alpha0_deeply_anoxic =rem_parms["alpha0"]["deeply_anoxic"],
                            reference_pH_suboxic =rem_parms["reference_pH"]["suboxic"],
                            reference_pH_anoxic =rem_parms["reference_pH"]["anoxic"],
                            reference_pH_deeply_anoxic =rem_parms["reference_pH"]["deeply_anoxic"],
                            mu1_suboxic = rem_parms["mu1"]["suboxic"],
                            mu1_anoxic = rem_parms["mu1"]["anoxic"],
                            mu1_deeply_anoxic = rem_parms["mu1"]["deeply_anoxic"],
                            organism_diam = rem_parms["organism_diam"]
                            )

    # Calculate advective microbial removal
    # Final concentration per endpoint_id
    C_final = {}
    for endpoint_id in modpath_phrea.schematisation_dict.get("endpoint_id"):
        df_particle, df_flowline, C_final[endpoint_id] = modpath_removal.calc_advective_microbial_removal(
                                            modpath_phrea.df_particle, modpath_phrea.df_flowline,
                                            endpoint_id = endpoint_id,
                                            trackingdirection = modpath_phrea.tracking_direction,
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

#=======
#%%

def test_modpath_run_semiconfined_nogravelpack_traveltimes(organism_name = "MS2"):


    ''' Phreatic scheme with gravelpack: modpath run.'''
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
                                      temperature=11.,
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

    modpath_semiconf = MP.ModPathWell(test_semiconf, # semiconf_dict_1, #test_phrea,
                            workspace = "test7_semiconf_nogp",
                            modelname = "semi_conf_nogp",
                            bound_left = "xmin",
                            bound_right = "xmax")

    # Run phreatic schematisation
    modpath_semiconf.run_model(run_mfmodel = True,
                        run_mpmodel = True)

    # Calculate advective microbial removal
    modpath_removal = ST.SubstanceTransport(modpath_semiconf,
                                            organism = organism_name,
                                            )
    # modpath_removal.compute_omp_removal()
    # Final concentration per endpoint_id
    C_final = {}
    for endpoint_id in modpath_semiconf.schematisation_dict.get("endpoint_id"):
        df_particle, df_flowline, C_final[endpoint_id] = modpath_removal.calc_advective_microbial_removal(
                                            modpath_semiconf.df_particle, modpath_semiconf.df_flowline,
                                            endpoint_id = endpoint_id,
                                            trackingdirection = modpath_semiconf.tracking_direction,
                                            mu1 = 0.023, grainsize = 0.00025, alpha0 = 1.E-5, reference_pH = 6.8, const_BM = 1.38e-23,
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

#%%

def test_diffuse_modpath_run_semiconfined_nogravelpack_traveltimes(organism_name = "MS2"):

    ''' Phreatic scheme without gravelpack: modpath run.'''
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
                                      temperature=11.,
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

    # Remove/empty point_parameters
    test_semiconf.point_parameters = {}

    modpath_semiconf = MP.ModPathWell(test_semiconf, #test_phrea,
                            workspace = "test8_semiconf",
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
                                      temperature=11.,
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

def test_travel_time_distribution_phreatic_analytical(organism_name = "MS2"):
    output_phreatic = pd.read_csv(path / 'phreatic_test.csv')
    output_phreatic = output_phreatic.round(7) #round to 7 digits (or any digit), keep same as for the output for the model to compare

    test_ = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                        computation_method= 'analytical',
                                        what_to_export='all',
                                        removal_function = 'mbo',
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
                                        name=organism_name
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

#%%

def test_travel_time_distribution_phreatic_analytical_plus_modpath(organism_name = "MS2"):
    output_phreatic = pd.read_csv(path / 'phreatic_test.csv')
    output_phreatic = output_phreatic.round(7) #round to 7 digits (or any digit), keep same as for the output for the model to compare

    test_phrea = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                        computation_method= 'analytical',
                                        what_to_export='all',
                                        removal_function = 'mbo',
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
                                        name=organism_name
                                        )

    # AnalyticalWell object
    well1_an = AW.AnalyticalWell(test_phrea)
    well1_an.phreatic()
    output_an = well1_an.df_particle
    output_an = output_an.round(7)

    # ModPath well object
    well1_mp = MP.ModPathWell(test_phrea, #test_phrea,
                            workspace = "test11_phrea",
                            modelname = "phreatic",
                            bound_left = "xmin",
                            bound_right = "xmax")

    # Run phreatic schematisation
    well1_mp.run_model(run_mfmodel = True,
                        run_mpmodel = True)

    output_mp = well1_mp.df_particle
    # output_mp = output_mp[["total_travel_time", "travel_time_unsaturated",
    #                  "travel_time_shallow_aquifer", "travel_time_target_aquifer",
    #                  "radial_distance", ]]
    output_mp = output_mp.round(7)

    # Export output of analytical and modpath calculations

    # df_particle file name (analytical well data)
    particle_fname_an = os.path.join(well1_mp.dstroot,well1_mp.schematisation_type + "_df_particle_analytical.csv")
    # Save df_particle
    output_an.to_csv(particle_fname_an)

    # df_particle file name
    particle_fname_mp = os.path.join(well1_mp.dstroot,well1_mp.schematisation_type + "_df_particle_modpath.csv")
    # Save df_particle
    output_mp.to_csv(particle_fname_mp)


    try:
        assert_frame_equal(output_an, output_phreatic,check_dtype=False)

    except AssertionError:
        print("Assertion Exception Raised - TTD test")
    else:
        print("Success, no error in TTD!")


#%%
# if __name__ == "__main__":
#     test_modpath_run_phreatic_withgravelpack()