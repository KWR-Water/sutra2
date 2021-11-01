
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
# try:
#     from project_path import module_path #the dot says looik in the current folder, this project_path.py file must be in the folder here
# except ModuleNotFoundError:
#     from project_path import module_path

# <<<<<<< HEAD
# Import schematisation functions
try:
    from greta.ModPath_functions import ModPathWell   
    import greta.Analytical_Well as AW
    from greta.Substance_Transport import *
except ModuleNotFoundError as e:
    print(e, ": second try.")
    module_path = os.path.join("..","greta")
    if module_path not in sys.path:
        sys.path.insert(0,module_path)

    from ModPath_functions import ModPathWell   
    import Analytical_Well as AW
    from Substance_Transport import *

    print("Second try to import modules succeeded.")

# from greta.Analytical_Well import *
# from greta.Substance_Transport import *
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
def test_modflow_run_phreatic():
    test_phrea = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                        computation_method= 'analytical', 
                                        what_to_export='omp',
                                        well_discharge=319.4*24,
                                        # vertical_resistance_shallow_aquifer=500,
                                        hor_permeability_shallow_aquifer = 0.02,
                                        vertical_anisotropy_shallow_aquifer = (10/(0.02*500)),
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
    test_phrea.make_dictionary()
    # print(test_phrea.__dict__)
    modpath_phrea = ModPathWell(test_phrea,
                            workspace = "test_ws",
                            modelname = "phreatic",
                            bound_left = "xmin",
                            bound_right = "xmax")
    # print(modpath_phrea.__dict__)
    # Run phreatic schematisation
    modpath_phrea.run_model(run_mfmodel = True,
                        run_mpmodel = False)
    assert modpath_phrea.success_mf
#%%
def test_modpath_run_phreatic():
    test_phrea = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                        computation_method= 'analytical', 
                                        what_to_export='omp',
                                        well_discharge=319.4*24,
                                        # vertical_resistance_shallow_aquifer=500,
                                        hor_permeability_shallow_aquifer = 0.02,
                                        vertical_anisotropy_shallow_aquifer = (10/(0.02*500)),
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
    test_phrea.make_dictionary()
    # print(test_phrea.__dict__)
    modpath_phrea = ModPathWell(test_phrea,
                            workspace = "test_ws",
                            modelname = "phreatic",
                            bound_left = "xmin",
                            bound_right = "xmax")
    # print(modpath_phrea.__dict__)
    # Run phreatic schematisation
    modpath_phrea.run_model(run_mfmodel = True,
                        run_mpmodel = True)
    print(modpath_phrea.material)
    parms = ["material","hk","vani","porosity","recharge"]
    for iParm in parms:
        np.save(os.path.join(testfiles_dir, iParm + "_phrea.npy"), getattr(modpath_phrea,iParm))
                     
    assert modpath_phrea.success_mp


#=======
#%%

def test_phreatic_scheme_dictinput():
    phreatic_scheme= AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                    computation_method = 'modpath',
                                    what_to_export='omp',
                                    # biodegradation_sorbed_phase = False,
                                      well_discharge=319.4*24,
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
                                      # diameter_filterscreen = 0.2,
                                      point_input_concentration = 100.,
                                      discharge_point_contamination = 100.,#made up value
                                      top_clayseal = 17,
                                      compute_contamination_for_date=dt.datetime.strptime('2020-01-01',"%Y-%m-%d"),

                                      # substance = 'benzene',
                                      # halflife_suboxic=600,
                                      # partition_coefficient_water_organic_carbon = 3.3,
                                    )

    phreatic_scheme.make_dictionary()  
    
    phreatic_dict_1 = { #'simulation_parameters' : phreatic_scheme.simulation_parameters,
            'endpoint_id': {phreatic_scheme.endpoint_id["name"]: phreatic_scheme.endpoint_id},
            'mesh_refinement': phreatic_scheme.mesh_refinement,
            'geo_parameters' : phreatic_scheme.geo_parameters,
            'ibound_parameters' : phreatic_scheme.ibound_parameters,
            'recharge_parameters' : phreatic_scheme.recharge_parameters,
            'well_parameters' : phreatic_scheme.well_parameters,
            # 'point_parameters' : phreatic_scheme.point_parameters,
            'substance_parameters' : phreatic_scheme.substance_parameters,
            'bas_parameters' : phreatic_scheme.bas_parameters,
    }
    ''' Hier gebleven: 1-11-21 --> uncomment 'simulation_parameters' en 'point_parameters'. '''
    fpath_research = os.path.abspath(os.path.join(path,os.pardir,"research"))
    fpath_phreatic_dict_check1 = os.path.join(fpath_research,"phreatic_dict_nogravel_test.txt")
    with open (fpath_phreatic_dict_check1, "w") as dict_file:
        dict_file.write(str(phreatic_dict_1))

    with open(fpath_phreatic_dict_check1,"r") as dict_file:
        dict_raw = dict_file.read()
        phreatic_dict_check1 = ast.literal_eval(dict_raw)  # ast --> abstract syntax trees
        # pd.read_csv(fpath_phreatic_dict_check1, delimiter=" ", header = None)
    
    assert True
    # assert phreatic_dict_1 == phreatic_dict_check1

#%%

def test_travel_time_distribution_phreatic():
    output_phreatic = pd.read_csv(path / 'phreatic_test.csv')
    output_phreatic = output_phreatic.round(7) #round to 7 digits (or any digit), keep same as for the output for the model to compare

    test_ = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                        computation_method= 'analytical',
                                        what_to_export='omp',
                                        well_discharge=319.4*24,
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

# %%
if __name__ == "__main__":
    test_modpath_run_phreatic()
