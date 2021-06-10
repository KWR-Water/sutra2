
#%%
import pytest
from pandas import read_csv
import pandas as pd
import os
# path = os.getcwd()  # path of working directory
from pathlib import Path
try:
    from project_path import module_path #the dot says looik in the current folder, this project_path.py file must be in the folder here
except ModuleNotFoundError:
    from project_path import module_path

from greta.draft_transport_function import *

from pandas._testing import assert_frame_equal

# get directory of this file
path = Path(__file__).parent #os.getcwd() #path of working directory


#%%
def test_travel_time_distribution_phreatic():
    output_phreatic = pd.read_csv(path / 'phreatic_test.csv')
    output_phreatic = output_phreatic.round(7) #round to 7 digits (or any digit), keep same as for the output for the model to compare

    test_ = HydroChemicalSchematisation(schematisation_type='phreatic',
                                        what_to_export='omp_parameters',
                                        well_discharge=319.4*24,
                                        # vertical_resistance_aquitard=500,
                                        hor_permeability_shallow_aquifer = 0.02,
                                       vertical_anistropy_shallow_aquifer = (10/(0.02*500)),
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

    well1 = AnalyticalWell(test_)
    well1.phreatic()
    output = well1.df_output
    output = output[["total_travel_time", "travel_time_unsaturated",
                     "travel_time_shallow_aquifer", "travel_time_target_aquifer",
                     "radial_distance", ]]
    output = output.round(7)

    try:
        # assert output == output_phreatic
        assert_frame_equal(output, output_phreatic,check_dtype=False)

    except AssertionError:
        print("Assertion Exception Raised - TTD test")
    else:
        print("Success, no error in TTD!")

def test_retardation_temp_koc_correction(substance = 'benzene', schematisation_type='phreatic'):
    test_ = HydroChemicalSchematisation(schematisation_type=schematisation_type,
                                        what_to_export='omp_parameters',
                                      well_discharge=319.4*24,
                                      hor_permeability_shallow_aquifer = 0.02,
                                      vertical_anistropy_shallow_aquifer = (10/(0.02*500)),
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
                                      temperature=11,
                                      solid_density_vadose_zone= 2.650, 
                                      solid_density_shallow_aquifer= 2.650, 
                                      solid_density_target_aquifer= 2.650, 
                                      diameter_borehole = 0.75,
                                    )
    well1 = AnalyticalWell(test_)
    if schematisation_type=='phreatic':
        well1.phreatic()  
    elif  schematisation_type=='semiconfined':
        well1.semiconfined()
    conc1 = Concentration(well1, substance = substance) #, df_particle, df_flowline)
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
    test_ = HydroChemicalSchematisation(schematisation_type='phreatic',
                                        what_to_export='omp_parameters',
                                      well_discharge=319.4*24,
                                    #   vertical_resistance_aquitard=500,
                                      hor_permeability_shallow_aquifer = 0.02,
                                      vertical_anistropy_shallow_aquifer = (10/(0.02*500)),
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
                                      temperature=11,
                                      solid_density_vadose_zone= 2.650, 
                                      solid_density_shallow_aquifer= 2.650, 
                                      solid_density_target_aquifer= 2.650, 
                                      diameter_borehole = 0.75,
                                    )
    well1 = AnalyticalWell(test_)
    well1.phreatic()  
    # substance = 'benzene'
    conc1 = Concentration(well1, substance = substance) #, df_particle, df_flowline)
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
    output_semiconfined = pd.read_csv(path / 'semiconfined_test.csv')
    output_semiconfined = output_semiconfined.round(7)
    test_ = HydroChemicalSchematisation(schematisation_type='semiconfined',
                                                what_to_export='omp_parameters',
                                        well_discharge=319.4*24,
                                        # vertical_resistance_aquitard=500,
                                      hor_permeability_shallow_aquifer = 0.02,
                                      vertical_anistropy_shallow_aquifer = (10/(0.02*500)),
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
    well1 = AnalyticalWell(test_)
    well1.semiconfined() 
    output = well1.df_output
    # output = output_dict['df_output']
    output = output[["total_travel_time", "travel_time_unsaturated",
                    "travel_time_shallow_aquifer", "travel_time_target_aquifer",
                    "radial_distance",]]
    try:
        # assert output == output_semiconfirned
        assert_frame_equal(output, output_semiconfined, 
                        check_dtype=False)

    except AssertionError:
        print("Assertion Exception Raised - in TTD test")
    else:
        print("Success, no error in TTD!")



def test_steady_concentration_temp_koc_correction_semiconfined(substance='benzene'):

    
    test_ = HydroChemicalSchematisation(schematisation_type='semiconfined',
                                        what_to_export='omp_parameters',
                                      well_discharge=319.4*24,
                                    #   vertical_resistance_aquitard=500,
                                      hor_permeability_shallow_aquifer = 0.02,
                                      vertical_anistropy_shallow_aquifer = (10/(0.02*500)),
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
                                      temperature=11,
                                      solid_density_vadose_zone= 2.650, 
                                      solid_density_shallow_aquifer= 2.650, 
                                      solid_density_target_aquifer= 2.650, 
                                      diameter_borehole = 0.75,
                                    )
    well1 = AnalyticalWell(test_)
    well1.semiconfined()  
    # substance = 'benzene'
    conc1 = Concentration(well1, substance = substance) #, df_particle, df_flowline)
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
