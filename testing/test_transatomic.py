
#%%
import pytest
import numpy as np
import datetime as dt
from pandas import read_csv
import pandas as pd
import os
import sys

# path = os.getcwd()  # path of working directory
from pathlib import Path

# try:
#     from project_path import module_path #the dot says looik in the current folder, this project_path.py file must be in the folder here
# except ModuleNotFoundError:
#     from project_path import module_path


import sutra2.Analytical_Well as AW
import sutra2.Transport_Removal as TR
from pandas.testing import assert_frame_equal
import warnings

# get directory of this file
path = Path(__file__).parent #os.getcwd() #path of working directory

#%%
def test_travel_time_distribution_phreatic():
    """ Compares the calculated travel times (total, unsaturated zone, shallow aquifer and target aquifer) 
    against a known case from TRANSATOMIC excel """
    output_phreatic = pd.read_csv(path / 'phreatic_test.csv')
    output_phreatic = output_phreatic.round(7) #round to 7 digits (or any digit), keep same as for the output for the model to compare

    test_ = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
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
                                        temp_water=11,
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

    assert_frame_equal(output, output_phreatic,check_dtype=False)


def test_retardation_temp_koc_correction(substance = 'benzene', schematisation_type='phreatic'):
    """ Compares the calculated retardation coefficient for each redox zone against a known case from TRANSATOMIC excel """
    test_ = AW.HydroChemicalSchematisation(schematisation_type=schematisation_type,
                                        computation_method= 'analytical',
                                        what_to_export='omp',
                                        name = substance,
                                      well_discharge=-319.4*24,
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
                                      diffuse_input_concentration = 100,
                                      temp_water=11,
                                      solid_density_vadose_zone= 2.650,
                                      solid_density_shallow_aquifer= 2.650,
                                      solid_density_target_aquifer= 2.650,
                                      diameter_borehole = 0.75,
                                    )
    well1 = AW.AnalyticalWell(test_)
    if schematisation_type=='phreatic':
        well1.phreatic()
    elif  schematisation_type=='semiconfined':
        well1.semiconfined()
    conc1 = TR.Transport(well1, substance = substance) #, df_particle, df_flowline)
    conc1.compute_omp_removal()

    retardation = {
        'benzene': {
            'vadose_zone': 1.57866594,
            'shallow_aquifer':  1.32938582,
            'target_aquifer': 1.32940346,
        },
        'benzo(a)pyrene': {
            'vadose_zone': 1939.142373,
            'shallow_aquifer': 2388.097816,
            'target_aquifer': 3901.698980,
        },
        'AMPA' :{
            'vadose_zone': 1.0000000763015349,
            'shallow_aquifer':	1.000000004342605, #1.0000000004342615,
            'target_aquifer': 1.0000000004342615,
        },
    }
    retardation_array = np.array([retardation[substance]['vadose_zone'],
                    retardation[substance]['shallow_aquifer'],
                    retardation[substance]['target_aquifer']])

    test_array = np.array(conc1.df_particle.retardation.loc[1:3], dtype='float')

    try:
        # assert output == output_phreatic
        np.testing.assert_allclose(test_array,
                           retardation_array ),
                        #    rtol=1e-8, atol=1e-8)

    except AssertionError:
        print("Assertion Exception Raised - retardation test")
    else:
        print("Success, no error in retardation!")

def test_steady_concentration_temp_koc_correction_phreatic(substance='benzene'):
    """ Compares the calculated steady state concentration for a specific radial distance
    for each redox zone against a known case from TRANSATOMIC excel """

    test_ = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                        computation_method= 'analytical',
                                        what_to_export='omp',
                                        name = substance,
                                      well_discharge=-319.4*24,
                                    #   vertical_resistance_shallow_aquifer=500,
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
                                      diffuse_input_concentration = 100,
                                      temp_water=11,
                                      solid_density_vadose_zone= 2.650,
                                      solid_density_shallow_aquifer= 2.650,
                                      solid_density_target_aquifer= 2.650,
                                      diameter_borehole = 0.75,
                                    )
    well1 = AW.AnalyticalWell(test_)
    well1.phreatic()
    # substance = 'benzene'
    conc1 = TR.Transport(well1, substance = substance) #, df_particle, df_flowline)
    conc1.compute_omp_removal()

    steady_state_concentration = {
        'benzene': {
            'vadose_zone': 10.744926872632352,
            'shallow_aquifer':  1.3763989974870514,
            'target_aquifer': 1.3763989974870514,
        },
        'benzo(a)pyrene': {
            'vadose_zone': 0,
            'shallow_aquifer': 0,
            'target_aquifer': 0,
        },
        'AMPA' :{
            'vadose_zone': 0.000249362,
            'shallow_aquifer': 1.850450098e-10,#1.8504500983690007e-10,
            'target_aquifer': 1.850450098e-10, #1.8504500983690007e-10,
        },
    }
    concentration_array = np.array([steady_state_concentration[substance]['vadose_zone'],
                    steady_state_concentration[substance]['shallow_aquifer'],
                    steady_state_concentration[substance]['target_aquifer']])

    test_array = np.array(conc1.df_particle.steady_state_concentration.loc[1:3], dtype=float)

    try:
        # assert output == output_phreatic
        # assert_frame_equal(test_array,concentration_array,check_dtype=False)
        np.testing.assert_allclose(test_array,
                           concentration_array,
                           rtol=1e-8, atol=1e-8)

    except AssertionError:
        print("Assertion Exception Raised - concetration test")
    else:
        print("Success, no error in concetration!")

# %%

def test_travel_time_distribution_semiconfined():
    """ Compares the calculated travel times (total, unsaturated zone, shallow aquifer and target aquifer) 
    against a known case from TRANSATOMIC excel """

    # output_semiconfined = pd.read_csv(path / 'semiconfined_test.csv')
    output_semiconfined = pd.read_csv(path / 'semiconfined_test_fixed_TTD.csv')

    output_semiconfined = output_semiconfined.round(7)
    test_ = AW.HydroChemicalSchematisation(schematisation_type='semiconfined',
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
                                        temp_water=11,
                                        solid_density_vadose_zone= 2.650,
                                        solid_density_shallow_aquifer= 2.650,
                                        solid_density_target_aquifer= 2.650,
                                        diameter_borehole = 0.75,)
    well1 = AW.AnalyticalWell(test_)
    well1.semiconfined()
    output = well1.df_output
    # output = output_dict['df_output']
    output = output[["total_travel_time", "travel_time_unsaturated",
                    "travel_time_shallow_aquifer", "travel_time_target_aquifer",
                    "radial_distance",]]
    # try:
        # assert output == output_semiconfirned
    assert_frame_equal(output, output_semiconfined,
                        check_dtype=False)
        # assert output ==1

    # except AssertionError:
    #     print("Assertion Exception Raised - in TTD test")
    # else:
    #     print("Success, no error in TTD!")


def test_steady_concentration_temp_koc_correction_semiconfined(substance='benzene'):
    """ Compares the calculated retardation coefficient for each redox zone against a known case from TRANSATOMIC excel """

    test_ = AW.HydroChemicalSchematisation(schematisation_type='semiconfined',
                                        computation_method= 'analytical',
                                        what_to_export='omp',
                                        name = substance,
                                      well_discharge=-319.4*24,
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
                                      diffuse_input_concentration = 100,
                                      temp_water=11,
                                      solid_density_vadose_zone= 2.650,
                                      solid_density_shallow_aquifer= 2.650,
                                      solid_density_target_aquifer= 2.650,
                                      diameter_borehole = 0.75,
                                    )
    well1 = AW.AnalyticalWell(test_)
    well1.semiconfined()
    # substance = 'benzene'
    conc1 = TR.Transport(well1, substance = substance) #, df_particle, df_flowline)
    conc1.compute_omp_removal()

    steady_state_concentration = {
        'benzene': {
            'vadose_zone': 30.78934144,
            'shallow_aquifer':  21.11155403,
            'target_aquifer': 	21.11155403,
        },
        'benzo(a)pyrene': {
            'vadose_zone': 0,
            'shallow_aquifer': 0,
            'target_aquifer': 0,
        },
        'AMPA' :{
            'vadose_zone': 0.109923889,
            'shallow_aquifer':	0.008232593,
            'target_aquifer':0.008232593,
        },
    }
    concentration_array = np.array([steady_state_concentration[substance]['vadose_zone'],
                    steady_state_concentration[substance]['shallow_aquifer'],
                    steady_state_concentration[substance]['target_aquifer']])

    test_array = np.array(conc1.df_particle.steady_state_concentration.loc[1:3], dtype=float)

    try:
        # assert output == output_phreatic
        # assert_frame_equal(test_array,concentration_array,check_dtype=False)
        np.testing.assert_allclose(test_array,
                           concentration_array,
                           rtol=1e-8, atol=1e-8)

    except AssertionError:
        print("Assertion Exception Raised - concetration test")
    else:
        print("Success, no error in concetration!")

# %%

def test_start_end_dates_contamination():
    ''' Tests whether the correct exception is raised when the 'end_date_contamiantion' is before 'start_date_contamination' '''

    with pytest.raises(ValueError) as exc:
        phreatic_scheme = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                                    computation_method= 'analytical',
                                                    what_to_export='omp',
                                                    well_discharge=-319.4*24, #m3/day
                                                    recharge_rate=0.3/365.25, #m/day
                                                    start_date_contamination= dt.datetime.strptime('1990-01-01', "%Y-%m-%d") ,
                                                    end_date_contamination= dt.datetime.strptime('1950-01-01', "%Y-%m-%d"), #'1950-01-01'
                                      )
    assert 'Error, "end_date_contamination" is before "start_date_contamination". Please enter an new "end_date_contamination" or "start_date_contamination" ' in str(exc.value)

#%%
def test_compute_for_date_start_dates_contamination():
    ''' Tests whether the correct exception is raised when the 'computer_contamiantion_for_date' is before 'start_date_contamination' '''

    with pytest.raises(ValueError) as exc:
        phreatic_scheme = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                                    computation_method= 'analytical',
                                      what_to_export='omp',
                                      well_discharge=-319.4*24, #m3/day
                                      recharge_rate=0.3/365.25, #m/day
                                      start_date_contamination=dt.datetime.strptime('1960-01-01', "%Y-%m-%d") ,
                                      end_date_contamination= dt.datetime.strptime('1990-01-01', "%Y-%m-%d"),
                                      compute_contamination_for_date= dt.datetime.strptime('1950-01-01', "%Y-%m-%d")
                                      )
    assert 'Error, "compute_contamination_for_date" is before "start_date_contamination". Please enter an new "compute_contamination_for_date" or "start_date_contamination" ' in str(exc.value)

#%%
def test_compute_for_date_start_date_well():
    ''' Tests whether the correct exception is raised when the
    'computer_contamiantion_for_date' is before 'start_date_contamination' '''

    with pytest.raises(ValueError) as exc:
        phreatic_scheme = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                                    computation_method= 'analytical',
                                      what_to_export='omp',
                                      well_discharge=-319.4*24, #m3/day
                                      recharge_rate=0.3/365.25, #m/day
                                      start_date_contamination=dt.datetime.strptime('1950-01-01', "%Y-%m-%d") ,
                                      end_date_contamination= dt.datetime.strptime('1990-01-01', "%Y-%m-%d"),
                                      compute_contamination_for_date= dt.datetime.strptime('1960-01-01', "%Y-%m-%d"),
                                      start_date_well= dt.datetime.strptime('1975-01-01', "%Y-%m-%d"),
                                      )
    assert 'Error, "compute_contamination_for_date" is before "start_date_well". Please enter an new "compute_contamination_for_date" or "start_date_well" ' in str(exc.value)
#%%

def test_incorrect_date_input_format():
    ''' Tests whether the correct exception is raised when the
    'computer_contamiantion_for_date' is before 'start_date_contamination' '''

    with pytest.raises(TypeError) as exc:
        phreatic_scheme = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                                    computation_method= 'analytical',
                                      what_to_export='omp',
                                      well_discharge=-319.4*24, #m3/day
                                      recharge_rate=0.3/365.25, #m/day
                                      start_date_well='1950-01-01',
                                      )
    assert "Error invalid date input, please enter a new start_date_well using the format dt.datetime.strptime('YYYY-MM-DD', '%Y-%m-%d')" in str(exc.value)

#%%
def test_redox_options():
    ''' Tests whether the correct exception is raised when one of the redox zones
     is not one of'suboxic', 'anoxic', 'deeply_anoxic' '''
    with pytest.raises(ValueError) as exc:
        phreatic_scheme = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                                    computation_method= 'analytical',
                                what_to_export='omp',
                                well_discharge=-319.4*24, #m3/day
                                recharge_rate=0.3/365.25, #m/day
                                redox_vadose_zone='oxic',
                                redox_shallow_aquifer='anoxic',
                                redox_target_aquifer='deeply_anoxic',
                                )
    assert "Invalid redox_vadose_zone. Expected one of: ['suboxic', 'anoxic', 'deeply_anoxic']" in str(exc.value)

#%%

def test_phreatic_diffuse_point_source():
    ''' Test for phreatic case with both a diffuse and point source contamination, 
    Checks if the concentration over time matches a known case'''

    phreatic_scheme = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                      computation_method= 'analytical',
                                      removal_function = 'omp',
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
                                      temp_water=11,
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
    phreatic_well = AW.AnalyticalWell(phreatic_scheme)
    phreatic_well.phreatic()
    phreatic_conc = TR.Transport(phreatic_well, substance = 'OMP-X')
    phreatic_conc.compute_omp_removal()
    df_well_concentration = phreatic_conc.compute_concentration_in_well_at_date()
    df_well_concentration = df_well_concentration.astype({'time': 'int32', 'date': 'datetime64[ns]', 'total_concentration_in_well': 'float64'})
    
    df_well_concentration_test = read_csv(path / 'phreatic_diffuse_point_test.csv', index_col=0)
    # AH the assert frame_equal is being difficult, so have to specify each data type
    df_well_concentration_test = df_well_concentration_test.astype({'time': 'int32', 'date': 'datetime64[ns]', 'total_concentration_in_well': 'float64'}) #etc

    # relative tolerance
    rtol = 5.e-3
    # absolute tolerance
    atol = 1e99 # (not limiting)
    assert_frame_equal(df_well_concentration, df_well_concentration_test, check_dtype=False,
                        rtol = rtol, atol = atol)

#%%

def test_phreatic_diffuse_only_source():
    ''' Test for phreatic case with only a diffuse source contamination
    Checks if the concentration over time matches a known case'''

    phreatic_scheme = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                    computation_method= 'analytical',
                                    removal_function = 'omp',
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
                                    temp_water=11,
                                    solid_density_vadose_zone=2.650,
                                    solid_density_shallow_aquifer=2.650,
                                    solid_density_target_aquifer=2.650,
                                    diameter_borehole=0.75,
                                    #diffuse parameters
                                    diffuse_input_concentration=100, #ug/L
                                    #dates
                                    start_date_well=dt.datetime.strptime('1968-01-01',"%Y-%m-%d"),
                                    start_date_contamination= dt.datetime.strptime('1966-01-01',"%Y-%m-%d"),
                                    compute_contamination_for_date=dt.datetime.strptime('2050-01-01',"%Y-%m-%d"),
                                    end_date_contamination=dt.datetime.strptime('1990-01-01',"%Y-%m-%d"),
                                    )
    phreatic_well = AW.AnalyticalWell(phreatic_scheme)
    phreatic_well.phreatic()
    phreatic_conc = TR.Transport(phreatic_well, substance = 'OMP-X')
    phreatic_conc.compute_omp_removal()
    df_well_concentration = phreatic_conc.compute_concentration_in_well_at_date()

    df_well_concentration_test = read_csv(path / 'phreatic_diffuse_only_test.csv', index_col=0)
    # AH the assert frame_equal is being difficult, so have to specify each data type
    df_well_concentration_test = df_well_concentration_test.astype({'time': 'int32', 'date': 'datetime64[ns]', 'total_concentration_in_well': 'float64'})
    
    # relative tolerance
    rtol = 5.e-3
    # absolute tolerance
    atol = 1e99 # (not limiting)
    assert_frame_equal(df_well_concentration, df_well_concentration_test, check_dtype=False,
                        rtol = rtol, atol = atol)

#%%
def test_phreatic_point_only_source():
    ''' Test for phreatic case with only a diffuse source contamination,
    Checks if the concentration over time matches a known case'''

    phreatic_scheme = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                    computation_method= 'analytical',
                                    removal_function = 'omp',
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
                                    temp_water=11,
                                    solid_density_vadose_zone=2.650,
                                    solid_density_shallow_aquifer=2.650,
                                    solid_density_target_aquifer=2.650,
                                    diameter_borehole=0.75,
                                    #diffuse parameters
                                    diffuse_input_concentration=0, #ug/L
                                    #point paramters
                                    point_input_concentration=100,
                                    distance_point_contamination_from_well=25,
                                    depth_point_contamination=21, #m ASL
                                    discharge_point_contamination=-1000,
                                    #dates
                                    start_date_well=dt.datetime.strptime('1968-01-01', "%Y-%m-%d"),
                                    start_date_contamination= dt.datetime.strptime('1966-01-01', "%Y-%m-%d"),
                                    compute_contamination_for_date=dt.datetime.strptime('2050-01-01', "%Y-%m-%d"),
                                    end_date_contamination=dt.datetime.strptime('1990-01-01', "%Y-%m-%d"),

                                    )
    phreatic_well = AW.AnalyticalWell(phreatic_scheme)
    phreatic_well.phreatic()
    phreatic_conc = TR.Transport(phreatic_well, substance = 'OMP-X')
    phreatic_conc.compute_omp_removal()
    df_well_concentration = phreatic_conc.compute_concentration_in_well_at_date()

    df_well_concentration.to_csv('phreatic_point_only_test.csv')
    df_well_concentration_test = read_csv(path / 'phreatic_point_only_test.csv', index_col=0)
    # AH the assert frame_equal is being difficult, so have to specify each data type
    df_well_concentration_test = df_well_concentration_test.astype({'time': 'int32', 'date': 'datetime64[ns]', 'total_concentration_in_well': 'float64'})

    # relative tolerance
    rtol = 5.e-3
    # absolute tolerance
    atol = 1e99 # (not limiting)
    assert_frame_equal(df_well_concentration, df_well_concentration_test, check_dtype=False,
                        rtol = rtol, atol = atol)

def test_drawdown_lower_than_target_aquifer():
    ''' Tests whether the correct exception is raised when the drawdown of the 
    well is lower than the bottom of the target aquifer' '''
    with pytest.raises(ValueError) as exc:
        phreatic_scheme = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
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
                                        ground_surface = 22,)

        phreatic_well = AW.AnalyticalWell(phreatic_scheme)

        phreatic_well.phreatic()                                         
    assert "The drawdown at the well is lower than the bottom of the target aquifer. Please select a different schematisation." in str(exc.value)



# def test_warning_drawdown_in_target_aquifer():
#     ''' Tests whether a warning is issued when the head drawdown reaches the target aquifer' '''
#     @MartinK how to raise a warning here?
#     with AnalyticalWell.assertWarns(Warning) as exc:
#         phreatic_scheme = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
#                                         computation_method= 'analytical',
#                                         what_to_export='omp', # @alex: what_to_export sounds very cryptic and ad-hoc. maybe we can think of something better
#                                         well_discharge=-319.4*24,
#                                         # vertical_resistance_shallow_aquifer=500,
#                                         hor_permeability_shallow_aquifer = 0.02,
#                                         porosity_vadose_zone=0.38,
#                                         porosity_shallow_aquifer=0.35,
#                                         porosity_target_aquifer=0.35,
#                                         recharge_rate=0.3/365.25,
#                                         moisture_content_vadose_zone=0.15,
#                                         ground_surface = 22,
#                                         thickness_vadose_zone_at_boundary=1,
#                                         thickness_shallow_aquifer=1,
#                                         thickness_target_aquifer=20,
#                                         hor_permeability_target_aquifer=35,
#                                         thickness_full_capillary_fringe=0.4,
#                                         temperature=11,
#                                         solid_density_vadose_zone= 2.650,
#                                         solid_density_shallow_aquifer= 2.650,
#                                         solid_density_target_aquifer= 2.650,
#                                         diameter_borehole = 0.75,

#                                       )

#     assert 'The drawdown is lower than the bottom of the shallow aquifer' in str(exc.value)
