
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
                                        well_discharge=319.4*24,
                                        vertical_resistance_aquitard=500,
                                        porosity_vadose_zone=0.38,
                                        porosity_shallow_aquifer=0.35,
                                        porosity_target_aquifer=0.35,
                                        recharge_rate=0.3/365.25,
                                        soil_moisture_content_vadose_zone=0.15,
                                      ground_surface = 22,
                                    #   groundwater_level_ASL = 17,
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
    output = well1.phreatic()
    output = output[["total_travel_time", "travel_time_unsaturated",
                     "travel_time_shallow_aquifer", "travel_time_target_aquifer",
                     "radial_distance", ]]
    output = output.round(7)

    try:
        # assert output == output_phreatic
        assert_frame_equal(output, output_phreatic,check_dtype=False)

    except AssertionError:
        print("Assertion Exception Raised.")
    else:
        print("Success, no error!")

# %%


def test_travel_time_distribution_semiconfined():
    output_semiconfined = pd.read_csv(path / 'semiconfined_test.csv')
    output_semiconfined = output_semiconfined.round(7)
    test_ = HydroChemicalSchematisation(schematisation_type='semi-confined',
                                                well_discharge=319.4,
                                                vertical_resistance_aquitard=500,
                                                porosity_vadose_zone=0.38,
                                                porosity_shallow_aquifer=0.35,
                                                porosity_target_aquifer=0.35,
                                                recharge_rate=0.3,
                                                soil_moisture_content_vadose_zone=0.15,
                                                groundwater_level_ASL = 17,
                                                thickness_vadose_zone=5,
                                                thickness_shallow_aquifer=10,
                                                thickness_target_aquifer=40,
                                                KD=1400,
                                                thickness_full_capillary_fringe=0.4,)
    well1 = AnalyticalWell(test_)
    output = well1.semiconfined() 
    # output = output_dict['df_output']
    output = output[["total_travel_time", "travel_time_unsaturated",
                    "travel_time_shallow_aquifer", "travel_time_target_aquifer",
                    "radial_distance",]]
    try:
        # assert output == output_semiconfirned
        assert_frame_equal(output, output_semiconfined, 
                        check_dtype=False)

    except AssertionError:
        print("Assertion Exception Raised.")
    else:
        print("Success, no error!")

# %%
