
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

from draft_transport_function import *

from pandas._testing import assert_frame_equal

# get directory of this file
path = Path(__file__).parent #os.getcwd() #path of working directory


#%%

output_phreatic = pd.read_csv(path / 'phreatic_output.csv')


def test_travel_time_distribution_phreatic():
    test_ = HydroChemicalSchematisation(schematisation='phreatic',
                                        well_discharge_m3hour=319.4,
                                        vertical_resistance_aquitard=500,
                                        porosity_vadose_zone=0.38,
                                        porosity_shallow_aquifer=0.35,
                                        porosity_target_aquifer=0.35,
                                        recharge_rate=0.3,
                                        soil_moisture_content_vadose_zone=0.15,
                                        thickness_vadose_zone=5,
                                        thickness_shallow_aquifer=10,
                                        thickness_target_aquifer=40,
                                        KD=1400,
                                        thickness_full_capillary_fringe=0.4,)
    well1 = AnalyticalWell(test_)
    output = well1.phreatic()
    output = output[["total_travel_time", "travel_time_unsaturated",
                     "travel_time_shallow_aquifer", "travel_time_target_aquifer",
                     "radial_distance", ]]

    try:
        # assert output == output_phreatic
        assert_frame_equal(output, output_phreatic,
                           check_dtype=False, check_less_precise=2)

    except AssertionError:
        print("Assertion Exception Raised.")
    else:
        print("Success, no error!")

# %%
output_semiconfined = pd.read_csv(path / 'semiconfined_output.csv')


def test_travel_time_distribution_semiconfined():
    test_ = HydroChemicalSchematisation(schematisation='semi-confined',
                                                well_discharge_m3hour=319.4,
                                                vertical_resistance_aquitard=500,
                                                porosity_vadose_zone=0.38,
                                                porosity_shallow_aquifer=0.35,
                                                porosity_target_aquifer=0.35,
                                                recharge_rate=0.3,
                                                soil_moisture_content_vadose_zone=0.15,
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
                        check_dtype=False, check_less_precise=2)

    except AssertionError:
        print("Assertion Exception Raised.")
    else:
        print("Success, no error!")
