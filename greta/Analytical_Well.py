#%% ----------------------------------------------------------------------------
# A. Hockin, March 2021
# KWR BO 402045-247
# ZZS verwijdering bodempassage
# AquaPriori - Transport Model
# With Martin Korevaar, Martin vd Schans, Steven Ros, Bas vd Grift
#
# Based on Stuyfzand, P. J. (2020). Predicting organic micropollutant behavior 
#               for 4 public supply well field types, with TRANSATOMIC Lite+
#               (Vol. 2). Nieuwegein, Netherlands.
# ------------------------------------------------------------------------------

#### CHANGE LOG ####

# things which must be checked indicated in comments with AH
# specific questions flagged for;
# @MartinvdS // @steven //@martinK

#@MartinK alternative to self.schematisaiton.schematisation... this becomes cumbersome the more classes we have

# LEFT OFF: MAY 27
# add a dict with user input substance aprams which checks and updates dict if needed
# still need to do the RST file, point injection and others....

# May 19
# moved calculation of travel time unsaturated zone (and linked functions/params)
# to the hydrochemicalSchematisation class (from the AnalyticalWell class). Working.
# added the updated list of parameters to the hydrochemical schematisation and 
# dicitonaries from excel sheet

#To Do
# RST file/explanation of how to use the model
# comments ah_todo
# point injection
# need to incorporate the point/diffuse if/else statements
# @MartinK the testing functionality went away... worked for a few days just fine then stopped working today

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
import datetime
from datetime import timedelta  

path = os.getcwd()  # path of working directory


#%%
class HydroChemicalSchematisation:

    """ Converts input parameters of AquaPriori GUI to a complete parameterisation
    
    Parameters
    ----------
    schematisation_type: string 
        #AH add definition here this is the kind of well yuo use...
        'phreatic', 'semiconfined', 'riverbankfiltration', 'basinfiltration'
    computation_method: string
        'analytical', 'modpath'
    removal_function: string
        'omp' -> get parameters for OMP
        'microbiology' -> get parameters for microbiology
    what_to_export: string
        'all', 'omp', 'microbiology'
    temp_correction_Koc: Bool
        Default: True
    temp_correction_halflife: Bool
    biodegradation_sorbed_phase: Bool
    compute_thickness_vadose_zone: Bool
    ground_surface: float
        meters above sea level (ASL)
    thickness_vadose_zone_at_boundary: float
        meters
    bottom_vadose_zone_at_boundary: float
        meters
    head_boundary: float
        mASL
    thickness_shallow_aquifer: float 
        m
    bottom_shallow_aquifer: float 
        mASL
    thickness_target_aquifer: float 
        m
    bottom_target_aquifer: float 
        mASL
    thickness_full_capillary_fringe: float 
        m
    porosity_vadose_zone: float 
    porosity_shallow_aquifer: float 
    porosity_target_aquifer: float 
    moisture_content_vadose_zone: float 
    solid_density_vadose_zone: float 
        kg/L
    solid_density_shallow_aquifer: float 
        kg/L
    solid_density_target_aquifer: float 
        kg/L
    fraction_organic_carbon_vadose_zone: float
    fraction_organic_carbon_shallow_aquifer: float
    fraction_organic_carbon_target_aquifer: float
    redox_vadose_zone: string
        'oxic','anoxic', 'deeply_anoxic;
    redox_shallow_aquifer: string     
        'oxic','anoxic', 'deeply_anoxic;
    redox_target_aquifer: string
            'oxic','anoxic', 'deeply_anoxic;
    dissolved_organic_carbon_vadose_zone: float 
        mg/L
    dissolved_organic_carbon_shallow_aquifer: float     
        mg/L
    dissolved_organic_carbon_target_aquifer: float 
        mg/L
    dissolved_organic_carbon_infiltration_water: float
        mg/L
    total_organic_carbon_infiltration_water: float 
        mg/L
    pH_vadose_zone: float 
    pH_shallow_aquifer: float
    pH_target_aquifer: float
    temperature: float 
        C
    temperature_vadose_zone: float 
        C
    temperature_shallow_aquifer: float 
        C
    temperature_target_aquifer: float 
        C
    recharge_rate: float 
        m/d
    """
    # initialize non-defaults to None to start...
    def __init__(self,
                schematisation_type= None, #'phreatic',

                computation_method= None, #'analytical',
                removal_function= 'omp',
                what_to_export = 'all',
                temp_correction_Koc=True,
                temp_correction_halflife=True,
                biodegradation_sorbed_phase=True,
                compute_thickness_vadose_zone=True,

                ground_surface=0.0,
                thickness_vadose_zone_at_boundary=1,
                bottom_vadose_zone_at_boundary=None,
                head_boundary=None,

                thickness_shallow_aquifer=10.0,
                thickness_target_aquifer=10.0,
                thickness_full_capillary_fringe=0.0,
                porosity_vadose_zone=0.35,
                porosity_shallow_aquifer=0.35,
                porosity_target_aquifer=0.35,
                moisture_content_vadose_zone=0.2,
                solid_density_vadose_zone=2.65, # kg/L
                solid_density_shallow_aquifer=2.65,
                solid_density_target_aquifer=2.65,
                fraction_organic_carbon_vadose_zone=0.001,
                fraction_organic_carbon_shallow_aquifer=0.0005,
                fraction_organic_carbon_target_aquifer=0.0005,

                redox_vadose_zone='oxic',
                redox_shallow_aquifer='anoxic',
                redox_target_aquifer='deeply_anoxic',
                dissolved_organic_carbon_vadose_zone=0.0,
                dissolved_organic_carbon_shallow_aquifer=0.0,
                dissolved_organic_carbon_target_aquifer=0.0,
                dissolved_organic_carbon_infiltration_water=0.0,
                total_organic_carbon_infiltration_water=0.0,
                pH_vadose_zone=7.0,
                pH_shallow_aquifer=7.0,
                pH_target_aquifer=7.0,
                temperature=10.0,
                temperature_vadose_zone=None,
                temperature_shallow_aquifer=None,
                temperature_target_aquifer=None,

                recharge_rate=0.001,
                well_discharge=-1000.0,
                basin_length=None, # BAR params \/
                basin_width=None,
                basin_xmin=None,
                basin_xmax=None,
                bottom_basin=None,
                bottom_sludge=None,
                hor_permeability_sludge=None,
                vertical_anistropy_sludge=None,
                head_basin=None,
                basin_infiltration_rate=None,
                travel_time_h20_shallow_aquifer=None,
                minimum_travel_time_h20_shallow_aquifer=None,
                travel_time_h20_deeper_aquifer=None,
                minimum_travel_time_h20_target_aquifer=None, # BAR params ^
                diameter_borehole=0.75,
                top_filterscreen=None,
                bottom_filterscreen=None,
                diameter_filterscreen=0.75,
                inner_diameter_filterscreen=None,
                top_gravelpack=None,
                bottom_gravelpack=None,
                diameter_gravelpack=None,
                inner_diameter_gravelpack=None,
                top_clayseal=None,
                bottom_clayseal=None,
                diameter_clayseal=None,
                hor_permeability_shallow_aquifer=1.0,
                hor_permeability_target_aquifer=10.0,
                hor_permebility_gravelpack=None,
                hor_permeability_clayseal=None,
                vertical_anistropy_shallow_aquifer=1.0,
                vertical_anistropy_target_aquifer=1.0,
                vertical_anistropy_gravelpack=None,
                vertical_anistropy_clayseal=None,

                substance=None,
                partition_coefficient_water_organic_carbon=None,
                dissociation_constant=None,
                halflife_oxic=None,
                halflife_anoxic=None,
                halflife_deeply_anoxic=None,

                diffuse_input_concentration=1,
                start_date_well='1950-01-01', #("Enter date in YYYY-MM-DD format")
                start_date_contamination=None,
                end_date_contamination='2000-01-01', #("Enter date in YYYY-MM-DD format")None, #AH my own default
                compute_contamination_for_date=None,

                concentration_point_contamination=None,
                distance_point_contamination_from_well=np.array([1]),
                depth_point_contamination=None,
                discharge_point_contamination=None, #AH_todo add this to the calculation-> what calculation???? how?

                source_area_radius=None,
                number_of_spreading_distance=None,
                model_radius=None,
                model_width=1, #AH @MartinvdS for the ycoord? cell width, see def _export_to_df(self, ):
                relative_position_starting_points_radial=None,
                relative_position_starting_points_in_basin=None,
                relative_position_starting_points_outside_basin=None,

                 # AH to be removed when if/else for point/diffuse added, see below
                #  recharge_concentration=None,
                #  input_concentration=None,
                 particle_release_date=None,
                
                #AH modpath params
                 ncols_near_well = 20,
                 ncols_far_well = 30,
                 nlayers_shallow_aquifer = None, 
                 nlayers_target_aquifer = None,

                 ):
        
        ''' Assign the parameters to be attributes of the class'''
        #AH_todo reorder these to have all the default/calcs done at the same place? or 
        # first just have them assigned then do something with them?

        # System
        self.schematisation_type = schematisation_type 

        # Settings
        self.computation_method = computation_method
        self.removal_function = removal_function
        self.what_to_export = what_to_export
        self.temp_correction_Koc = temp_correction_Koc
        self.temp_correction_halflife = temp_correction_halflife
        self.biodegradation_sorbed_phase = biodegradation_sorbed_phase
        self.compute_thickness_vadose_zone = compute_thickness_vadose_zone

        # Porous Medium
        self.ground_surface = ground_surface #as meters above sea level (m ASL)
        self.thickness_vadose_zone_at_boundary = thickness_vadose_zone_at_boundary
        self.bottom_vadose_zone_at_boundary = ground_surface - thickness_vadose_zone_at_boundary 
        self.head_boundary = head_boundary

        self.thickness_shallow_aquifer = thickness_shallow_aquifer
        self.bottom_shallow_aquifer = ground_surface - thickness_vadose_zone_at_boundary - thickness_shallow_aquifer
        self.thickness_target_aquifer = thickness_target_aquifer
        self.bottom_target_aquifer = self.bottom_shallow_aquifer - thickness_target_aquifer
        self.thickness_full_capillary_fringe = thickness_full_capillary_fringe
        self.porosity_vadose_zone = porosity_vadose_zone
        self.porosity_shallow_aquifer = porosity_shallow_aquifer
        self.porosity_target_aquifer = porosity_target_aquifer
        self.moisture_content_vadose_zone = moisture_content_vadose_zone
        self.solid_density_vadose_zone = solid_density_vadose_zone
        self.solid_density_shallow_aquifer = solid_density_shallow_aquifer
        self.solid_density_target_aquifer = solid_density_target_aquifer
        self.fraction_organic_carbon_vadose_zone = fraction_organic_carbon_vadose_zone
        self.fraction_organic_carbon_shallow_aquifer = fraction_organic_carbon_shallow_aquifer
        self.fraction_organic_carbon_target_aquifer = fraction_organic_carbon_target_aquifer

        # Hydrochemistry
        self.redox_vadose_zone = redox_vadose_zone
        self.redox_shallow_aquifer = redox_shallow_aquifer
        self.redox_target_aquifer = redox_target_aquifer
        self.dissolved_organic_carbon_vadose_zone = dissolved_organic_carbon_vadose_zone
        self.dissolved_organic_carbon_shallow_aquifer = dissolved_organic_carbon_shallow_aquifer
        self.dissolved_organic_carbon_target_aquifer = dissolved_organic_carbon_target_aquifer
        self.dissolved_organic_carbon_infiltration_water = dissolved_organic_carbon_infiltration_water
        self.total_organic_carbon_infiltration_water = total_organic_carbon_infiltration_water
        self.pH_vadose_zone = pH_vadose_zone
        self.pH_shallow_aquifer = pH_shallow_aquifer
        self.pH_target_aquifer = pH_target_aquifer
        self.temperature = temperature
        self.temperature_vadose_zone = temperature_vadose_zone
        self.temperature_shallow_aquifer = temperature_shallow_aquifer
        self.temperature_target_aquifer = temperature_target_aquifer

        # Hydrology
        self.recharge_rate = recharge_rate
        self.well_discharge = well_discharge
        self.basin_length = basin_length
        self.basin_width = basin_width
        self.basin_xmin = basin_xmin
        self.basin_xmax = basin_xmax
        self.bottom_basin = bottom_basin
        self.bottom_sludge = bottom_sludge
        self.hor_permeability_sludge = hor_permeability_sludge
        self.vertical_anistropy_sludge = vertical_anistropy_sludge
        self.head_basin = head_basin
        self.basin_infiltration_rate = basin_infiltration_rate
        self.travel_time_h20_shallow_aquifer = travel_time_h20_shallow_aquifer
        self.minimum_travel_time_h20_shallow_aquifer = minimum_travel_time_h20_shallow_aquifer
        self.travel_time_h20_deeper_aquifer = travel_time_h20_deeper_aquifer
        self.minimum_travel_time_h20_target_aquifer = minimum_travel_time_h20_target_aquifer
        self.diameter_borehole = diameter_borehole
        self.top_filterscreen = top_filterscreen
        self.bottom_filterscreen = bottom_filterscreen
        self.diameter_filterscreen = diameter_filterscreen
        self.inner_diameter_filterscreen = inner_diameter_filterscreen
        self.top_gravelpack = top_gravelpack
        self.bottom_gravelpack = bottom_gravelpack
        self.diameter_gravelpack = diameter_gravelpack
        self.inner_diameter_gravelpack = inner_diameter_gravelpack
        self.top_clayseal = top_clayseal
        self.bottom_clayseal = bottom_clayseal
        self.diameter_clayseal = diameter_clayseal
        self.hor_permeability_shallow_aquifer = hor_permeability_shallow_aquifer
        self.hor_permeability_target_aquifer = hor_permeability_target_aquifer
        self.hor_permebility_gravelpack = hor_permebility_gravelpack
        self.hor_permeability_clayseal = hor_permeability_clayseal
        self.vertical_anistropy_shallow_aquifer = vertical_anistropy_shallow_aquifer
        self.vertical_anistropy_target_aquifer = vertical_anistropy_target_aquifer
        self.vertical_anistropy_gravelpack = vertical_anistropy_gravelpack
        self.vertical_anistropy_clayseal = vertical_anistropy_clayseal

        self.vertical_resistance_aquitard = thickness_shallow_aquifer / (hor_permeability_shallow_aquifer *vertical_anistropy_shallow_aquifer)

        # Substance
        self.substance = substance
        self.partition_coefficient_water_organic_carbon = partition_coefficient_water_organic_carbon
        self.dissociation_constant = dissociation_constant
        self.halflife_oxic = halflife_oxic
        self.halflife_anoxic = halflife_anoxic
        self.halflife_deeply_anoxic = halflife_deeply_anoxic

        # Diffuse contamination override if point contamination specified
        # if concentration_point_contamination is None:
        self.diffuse_input_concentration = diffuse_input_concentration
            # self.input_concentration = diffuse_input_concentration

        # else: 
            # self.diffuse_input_concentration = 0 
            # self.input_concentration = concentration_point_contamination
        self.concentration_point_contamination = concentration_point_contamination

        def date_to_datetime(date):
            '''convert str input of date to datetime'''
            year, month, day = map(int, date.split('-'))
            date = datetime.date(year, month, day)
            return date

        self.start_date_well = date_to_datetime(start_date_well)
        self.end_date_contamination = date_to_datetime(end_date_contamination)

        # Contamination
        if start_date_contamination is None: 
            self.start_date_contamination = self.start_date_well
        else: 
            self.start_date_contamination = date_to_datetime(start_date_contamination)
        if compute_contamination_for_date is None: 
            self.compute_contamination_for_date = self.start_date_well + timedelta(days=100)
        else:
            self.compute_contamination_for_date = date_to_datetime(compute_contamination_for_date)
        if depth_point_contamination is None: 
            self.depth_point_contamination = self.ground_surface
        elif depth_point_contamination > self.ground_surface: 
            self.depth_point_contamination = self.ground_surface
        else:
            self.depth_point_contamination = depth_point_contamination


        # Point Contamination
        self.concentration_point_contamination = concentration_point_contamination
        self.distance_point_contamination_from_well = np.array([distance_point_contamination_from_well])
        self.discharge_point_contamination = discharge_point_contamination

        # Model size
        self.source_area_radius = source_area_radius
        self.number_of_spreading_distance = number_of_spreading_distance
        self.model_radius = model_radius
        self.model_width = model_width
        self.relative_position_starting_points_radial = relative_position_starting_points_radial
        self.relative_position_starting_points_in_basin = relative_position_starting_points_in_basin
        self.relative_position_starting_points_outside_basin = relative_position_starting_points_outside_basin

        ''' Default? calculations, calcs which should be done for all '''
        self.KD = hor_permeability_target_aquifer*thickness_target_aquifer
        self.groundwater_level =self.ground_surface-self.thickness_vadose_zone_at_boundary 
        
        # Temperature 
        if self.temperature_vadose_zone is None:
            self.temperature_vadose_zone = self.temperature
        if self.temperature_shallow_aquifer is None:
            self.temperature_shallow_aquifer = self.temperature
        if self.temperature_target_aquifer is None:
            self.temperature_target_aquifer = self.temperature


        # LEFT OFF HERE, NEED TO ADD SOMETHING (IF/ELSE) TO DECIDE WHICH CONCENTRATION
        # DIFFUSE/POINT TO ASSIGN AS THE CONCENTRATION
        # self.recharge_concentration = recharge_concentration
        # self.input_concentration = input_concentration
        if particle_release_date is None:
            self.particle_release_date = particle_release_date
        else: 
            self.particle_release_date =  date_to_datetime(particle_release_date)

        #Modpath params
        self.ncols_near_well = ncols_near_well
        self.ncols_far_well = ncols_far_well
        if nlayers_shallow_aquifer is None:
            self.nlayers_shallow_aquifer = int(self.thickness_shallow_aquifer)
        else: 
            self.nlayers_shallow_aquifer =nlayers_shallow_aquifer

        if nlayers_target_aquifer is None:
            self.nlayers_target_aquifer = int(self.thickness_target_aquifer)
        else: 
            self.nlayers_target_aquifer = nlayers_target_aquifer

        # if computation_method == 'modpath':
          
        top_shallow_aquifer = self.bottom_vadose_zone_at_boundary
        top_target_aquifer = self.bottom_shallow_aquifer

        if top_filterscreen	is None:
            self.top_filterscreen = top_target_aquifer
        if bottom_filterscreen is None: 
            self.bottom_filterscreen = self.bottom_target_aquifer
        if top_gravelpack is None: 
            self.top_gravelpack = top_target_aquifer
        if bottom_gravelpack is None: 
            self.bottom_gravelpack = self.bottom_target_aquifer
        if top_clayseal is None: 
            self.top_clayseal = self.ground_surface
        if bottom_clayseal is None: 
            self.bottom_clayseal = top_target_aquifer
        
        # DIAMETERS
        if diameter_gravelpack is None:
            self.diameter_gravelpack = self.diameter_borehole
        if inner_diameter_gravelpack is None:
            self.inner_diameter_gravelpack = diameter_filterscreen
        if inner_diameter_filterscreen is None:
            self.inner_diameter_filterscreen = self.diameter_filterscreen
        if diameter_clayseal is None:
            self.diameter_clayseal = diameter_borehole
            self.inner_diameter_clayseal = self.diameter_filterscreen

        #other default params
        if hor_permebility_gravelpack is None: 
            self.hor_permebility_gravelpack = 1000
        if hor_permeability_clayseal is None: 
            self.hor_permeability_clayseal = 0.001
        if vertical_anistropy_gravelpack is None: 
            self.vertical_anistropy_gravelpack = 1
        if vertical_anistropy_clayseal is None: 
            self.vertical_anistropy_clayseal = 1

        if model_radius is None: 
            if self.schematisation_type == 'phreatic':
                self.model_radius = (math.sqrt(self.well_discharge
                                                / (math.pi * self.recharge_rate ))) #AH SA*recharge = well_discharge
            elif self.schematisation_type == 'semiconfined':
                self.model_radius = math.sqrt(self.vertical_resistance_aquitard * self.KD * 3) # spreading_distance*3

        # check when initialization
        # thickness <=0 check, for all but the vadose zone
        # vadose zone can be 0, but the others should not


    def make_dictionary(self,):
        ''' Dicitonaries of the different parameters '''
        
        # @MartinvdS -> is this how to implement these for the different 
        # schematisations?
        if self.schematisation_type == 'phreatic':
            compute_thickness_vadose_zone = True # @MartinvdS what is this for?
            
            # additional meter added to the model radius for the fixed head boundary
            # only for the phreatic case, not for the semiconfined case
            self.model_radius_computed = self.model_radius + 1

            # only outer_boundary for phreatic 
            ibound_parameters = {
                'inner_boundary_shallow_aquifer': {
                    'top': self.bottom_vadose_zone_at_boundary,
                    'bot': self.bottom_shallow_aquifer,
                    'xmin': 0, 
                    'xmax': self.diameter_filterscreen/2, 
                    'ibound':0,
                    },

                'inner_boundary_target_aquifer': {
                    'head': self.bottom_vadose_zone_at_boundary,
                    'top': self.bottom_shallow_aquifer,
                    'bot': self.bottom_target_aquifer,
                    'xmin': 0, 
                    'xmax': self.diameter_filterscreen/2,
                    'ibound':-1,
                    }
                }

        elif self.schematisation_type == 'semiconfined':
            compute_thickness_vadose_zone = False # @MartinvdS what is this for?

            # ibound at the model radius (no additional meter added)
            self.model_radius_computed = self.model_radius

            # only top_boundary for semiconfined 
            ibound_parameters = {
                'top_boundary1': {
                    'head': self.bottom_vadose_zone_at_boundary,
                    'top': self.bottom_vadose_zone_at_boundary + 0.1,# 10 cm ficticous thickness to allow head boundary
                    'bottom': self.bottom_vadose_zone_at_boundary,
                    'xmin': self.diameter_gravelpack/2,
                    'xmax': self.model_radius_computed,
                    'ibound': -1, 
                        },
                    'inner_boundary_shallow_aquifer': {
                    'top': self.bottom_vadose_zone_at_boundary,
                    'bot':self.bottom_shallow_aquifer,
                    'xmin':0,
                    'xmax':self.diameter_filterscreen/2, #check this
                    'ibound':0, 
                    }
                    }

        # make dictionaries of each grouping of parameters
        simulation_parameters = {
            'schematisation_type': self.schematisation_type,
            'computation_method': self.computation_method,
            'temp_correction_Koc': self.temp_correction_Koc,
            'temp_correction_halflife': self.temp_correction_halflife,
            'biodegradation_sorbed_phase': self.biodegradation_sorbed_phase,
            'compute_thickness_vadose_zone': self.compute_thickness_vadose_zone,
            'start_date_well': self.start_date_well,
            'start_date_contamination': self.start_date_contamination,
            'compute_contamination_for_date': self.compute_contamination_for_date,
            }

        # Aquifer parameters dcitionary
        geo_parameters  = {
            'vadose': {
                'vadose': True,
                'top': self.ground_surface,
                'bot': self.bottom_vadose_zone_at_boundary,
                'xmin': self.diameter_gravelpack/2, 
                'xmax': self.model_radius,
                'porosity': self.porosity_vadose_zone,
                'moisture_content': self.moisture_content_vadose_zone,
                'solid_density': self.solid_density_vadose_zone,
                'f_oc': self.fraction_organic_carbon_vadose_zone,
                'redox': self.redox_vadose_zone,
                'DOC': self.dissolved_organic_carbon_vadose_zone,
                'pH': self.pH_vadose_zone,
                'T': self.temperature_vadose_zone,
                },
            'layer1': {
                'top': self.bottom_vadose_zone_at_boundary,
                'bot': self.bottom_shallow_aquifer,
                'xmin': self.diameter_gravelpack/2, 
                'xmax': self.model_radius_computed,
                'porosity': self.porosity_shallow_aquifer,
                'solid_density': self.solid_density_shallow_aquifer,
                'f_oc': self.fraction_organic_carbon_shallow_aquifer,
                'redox': self.redox_shallow_aquifer,
                'DOC': self.dissolved_organic_carbon_shallow_aquifer,
                'pH': self.pH_shallow_aquifer,
                'T': self.temperature_shallow_aquifer,
                'hk': self.hor_permeability_shallow_aquifer,
                'vani': self.vertical_anistropy_shallow_aquifer,
                'nlayers': self.nlayers_shallow_aquifer, 
                },
            'layer2': {
                'top': self.bottom_shallow_aquifer,
                'bot': self.bottom_target_aquifer,
                'xmin': self.diameter_gravelpack/2,
                'xmax': self.model_radius_computed,
                'porosity': self.porosity_target_aquifer,
                'solid_density': self.solid_density_target_aquifer,
                'f_oc': self.fraction_organic_carbon_target_aquifer,
                'redox': self.redox_target_aquifer,
                'DOC': self.dissolved_organic_carbon_target_aquifer,
                'pH': self.pH_target_aquifer,
                'T': self.temperature_target_aquifer,
                'hk': self.hor_permeability_target_aquifer,
                'vani': self.vertical_anistropy_target_aquifer,
                'nlayers': self.nlayers_target_aquifer, #AH maybe this will change to vertical_resolution
                },
            'gravelpack1': {
                'top': self.top_gravelpack,
                'bot': self.bottom_gravelpack,
                'xmin': self.inner_diameter_gravelpack/2,
                'xmax': self.diameter_gravelpack/2,
                'hk': self.hor_permebility_gravelpack,
                'vani': self.vertical_anistropy_gravelpack,
                },
            'clayseal1':{
                'top': self.top_clayseal,
                'bot': self.bottom_clayseal,
                'xmin': self.inner_diameter_clayseal/2, #@MartinvdS correct?
                'xmax': self.diameter_clayseal/2, 
                'hk': self.hor_permeability_clayseal,
                'vani': self.vertical_anistropy_clayseal,
                },
            # From well to borehole

            # From borehole to D target aquifer
            'mesh_refinement1': {
                'xmin': self.diameter_borehole/2, 
                'xmax': self.thickness_target_aquifer,
                'ncols': self.ncols_near_well, #indicates the number of columns close to the well
                },
            #From D target aquifer to model radiu #horizontal resolution
            'mesh_refinement2': {
                'xmin': self.thickness_target_aquifer, #@Martin from email... correct? self.diameter_gravelpack, 
                'xmax': self.model_radius, #mesh boundary at the model raidus, must line up AH
                'ncols': self.ncols_far_well #indicates the number of columns far from the well
                }, 
            }

        if self.schematisation_type == 'semiconfined':
            well_parameters = {
                'well1': {
                    'Q': self.well_discharge,
                    'top': self.top_filterscreen,
                    'bot': self.bottom_filterscreen,
                    'xmin': 0.0, 
                    'xmax': self.diameter_filterscreen/2,
                    },
                } 
        # elif will come here for the BAR, RBF cases AH
        else: 
            # no well_paramters dict needed for the phreatic case
            well_parameters = {
            }  

        recharge_parameters = {
            'source1': { # source1 -> recharge & diffuse sources
                'substance_name': self.substance,
                'recharge': self.recharge_rate,
                'xmin': self.diameter_gravelpack/2,
                'xmax': self.model_radius,
                'DOC': self.dissolved_organic_carbon_infiltration_water,
                'TOC': self.total_organic_carbon_infiltration_water,
                'c_in': self.diffuse_input_concentration, 
                },
            # 'source2' :{}> surface water (BAR & RBF) #@MartinvdS come back to this when we start this module
        }

        # Create point diciontary if point source concentration specified, 
        # otherwise pass empty dictionary
        if self.concentration_point_contamination is None:
            point_parameters= {}
        else:
            point_parameters= {
                'point1': { 
                    'substance_name': self.substance,
                    'c_in': self.concentration_point_contamination, 
                    'x_start': self.distance_point_contamination_from_well[0][0], #AH_todo check this works...
                    'z_start': self.depth_point_contamination,
                    'q_point': self.discharge_point_contamination,
                    },
                # 'point2': {} #AH_todo if there is more than one point source, 
                # input a dictionary directly, @maritnK, @martinvdS we will discuss this
                # and make a decision
                }

        #AH eventially to be computed by QSAR"
        substance_parameters = {
                'substance_name': self.substance,
                'log_Koc': self.partition_coefficient_water_organic_carbon,
                'pKa': self.dissociation_constant,
                'omp_half_life': {
                    'suboxic': self.halflife_oxic,
                    'anoxic': self.halflife_anoxic,
                    'deeply_anoxic': self.halflife_deeply_anoxic,
                    },
            }

        # @MartinvdS, the rest of the params are in the Modpath class, 
        # so make dictionary here or not?
        # bas = basic package modflow
        bas_parameters = {
            }

        # returned as attribute of function
        self.simulation_parameters = simulation_parameters
        self.geo_parameters = geo_parameters
        self.ibound_parameters = ibound_parameters
        self.recharge_parameters = recharge_parameters
        self.well_parameters = well_parameters
        self.point_parameters = point_parameters
        self.substance_parameters = substance_parameters
        self.bas_parameters = bas_parameters

    # functions to calculate the travel time through vadose zone, shared functions for
    # Analytical and Modflow models
    def _create_radial_distance_array(self):

        ''' Calculate the radial distance from a well,
        specifying the aquifer type (semiconfined, phreatic)
        '''
        #ah_todo change this to single array of 0.001 to 100
        fraction_flux = np.array([0.00001, 0.0001, 0.001, 0.005])
        fraction_flux = np.append(fraction_flux, np.arange(0.01, 1, 0.01))
        self.fraction_flux = np.append(fraction_flux, [0.995, 0.9999])

        radial_distance = self.radial_distance_recharge * \
            np.sqrt(self.fraction_flux)
            
        if self.schematisation_type == 'semiconfined':

            radial_distance = np.append(radial_distance,
                                    [(radial_distance[-1] + ((self.spreading_distance * 3) - radial_distance[-1]) / 3),
                                    (radial_distance[-1] + 2 *
                                    ((self.spreading_distance * 3) - radial_distance[-1]) / 3),
                                        (self.spreading_distance * 3)])
        
        self.radial_distance = radial_distance

    def calculate_hydraulic_head_phreatic(self, distance):
        head = (self.groundwater_level - self.well_discharge
                    / (2 * math.pi * self.KD) 
                    * np.log(self.radial_distance_recharge / distance))
        return head


    #ah_todo if time, split this into making the thickness, head and travel time (3 functions)
    def calculate_unsaturated_zone_travel_time_phreatic (self, 
                                                        distance, 
                                                        depth_point_contamination=None):
        ''''''

        head = self.calculate_hydraulic_head_phreatic(distance=distance)

        if depth_point_contamination is None:
            thickness_vadose_zone_drawdown = (self.groundwater_level 
                                                + self.thickness_vadose_zone_at_boundary) - head

            travel_time_unsaturated = (((self.groundwater_level + thickness_vadose_zone_drawdown
                                    - self.groundwater_level
                                    - self.thickness_full_capillary_fringe)
                                    * self.moisture_content_vadose_zone
                                    + self.thickness_full_capillary_fringe
                                    * self.porosity_vadose_zone)
                                / self.recharge_rate)

        elif depth_point_contamination >= head:
            thickness_vadose_zone_drawdown = depth_point_contamination - head
            if thickness_vadose_zone_drawdown < 0:
                travel_time_unsaturated =  np.array([0])
            else:
                travel_time_unsaturated = (((self.groundwater_level + thickness_vadose_zone_drawdown
                                    - self.groundwater_level
                                    - self.thickness_full_capillary_fringe)
                                    * self.moisture_content_vadose_zone
                                    + self.thickness_full_capillary_fringe
                                    * self.porosity_vadose_zone)
                                / self.recharge_rate)
        else:
            travel_time_unsaturated = np.array([0])
            thickness_vadose_zone_drawdown = 0 #AH_todo if time, replace this with the travel distance, not thickness_vadose because this is a stand in for the travel distance

        return travel_time_unsaturated, thickness_vadose_zone_drawdown, head 


    def _calculate_travel_time_unsaturated_zone(self, 
                                                distance = None, 
                                                depth_point_contamination=None):

        ''' Calculating the unsaturated zone travel time'''

        self.spreading_distance = math.sqrt(self.vertical_resistance_aquitard * self.KD)

        # AH do not change to model_radius, since the radial distance for recharge is based on the phreatic value for BOTH cases
        self.radial_distance_recharge =  (math.sqrt(self.well_discharge
                                                    / (math.pi * self.recharge_rate )))
        
        # Diffuse source or regular travel time calculation use radial distance array
        # Point source use specified distance or array? or distances? #AH_todo if we use multiple point sources
        if distance is None:
            self._create_radial_distance_array()
            distance= self.radial_distance
        else:
            distance = distance

        '''Equation A.11 in report '''

        if self.schematisation_type =='phreatic':
            travel_time_unsaturated, self.thickness_vadose_zone_drawdown, self.head = self.calculate_unsaturated_zone_travel_time_phreatic (distance= distance,depth_point_contamination=depth_point_contamination)
                                    
        elif self.schematisation_type == 'semiconfined':
            if depth_point_contamination is None:
                travel_distance = self.ground_surface - self.groundwater_level - self.thickness_full_capillary_fringe
            else: 
                #if point contamination at depth, assign ground surface to depth
                travel_distance =  depth_point_contamination - self.groundwater_level - self.thickness_full_capillary_fringe

            travel_time_unsaturated =(((travel_distance)
                                        * self.moisture_content_vadose_zone
                                        + self.thickness_full_capillary_fringe
                                        * self.porosity_vadose_zone)
                                    / self.recharge_rate)

            if travel_distance < 0:
                travel_time_unsaturated = 0

            # travel time in semiconfined is one value, make it array by repeating the value
            if isinstance(distance, float):
                pass
            else:
                travel_time_unsaturated = [travel_time_unsaturated] * (len(distance))
        
        self.travel_time_unsaturated = travel_time_unsaturated
    
class AnalyticalWell():
    """ Compute travel time distribution using analytical well functions."""

    def __init__(self, schematisation: HydroChemicalSchematisation): #change schematisation_instance to schematisation
        """ 'unpack/parse' all the variables from the hydrogeochemical schematizization """

        '''Parameters
        ----------
        df_flowline: pandas.DataFrame
            Column 'flowline_id': Integer
            Column 'discharge': Float
                Discharge associated with the flowline (m3/d)
            Column 'particle_release_date': Float
            Column 'input_concentration'
            Column 'endpoint_id': Integer
                ID of Well (or drain) where the flowline ends

        df_particle: pandas.DataFrame
            Column 'flowline_id'
            Column 'travel_time'
            Column 'xcoord'
            Column 'ycoord'
            Column 'zcoord'
            Column 'redox_zone'
            Column 'temperature'
        '''

        # get the non-default parameters
        # self.test_variable = None #AH test variable here to see if errors are caught....

        self.schematisation = schematisation
        
        #Make dictionaries
        # @MartinK if we make the dictionaries for all cases, do we need it as a seprate function?
        self.schematisation.make_dictionary()

    def _check_required_variables(self,required_variables):
        for req_var in required_variables:
            value = getattr(self.schematisation, req_var)
            if value is None:
                raise KeyError(f'Error, required variable {req_var} is not defined.')

    def _check_init_phreatic(self):
        '''check the variables that we need for the individual aquifer types are not NONE aka set by the user'''
        #AH_todo update these
        required_variables = ["schematisation_type", #repeat for all
                              "thickness_vadose_zone_at_boundary",
                              "thickness_shallow_aquifer",
                              "thickness_target_aquifer",
                              "porosity_vadose_zone",
                              "porosity_shallow_aquifer",
                              "porosity_target_aquifer",
                              "well_discharge",
                              "recharge_rate",
                              "vertical_resistance_aquitard",
                              "moisture_content_vadose_zone",
                              "KD",
                              "groundwater_level",
                              "thickness_full_capillary_fringe",
                            #   self.test_variable,
                            ]
        self._check_required_variables(required_variables)

    def _calculate_travel_time_shallow_aquifer_phreatic(self, 
                                                        head, 
                                                        depth_point_contamination=None,
                                                        ):

        if depth_point_contamination is None:
            travel_distance_shallow_aquifer = self.schematisation.thickness_shallow_aquifer - (self.schematisation.groundwater_level - head)
            travel_time_shallow_aquifer = ((travel_distance_shallow_aquifer)
                            * self.schematisation.porosity_shallow_aquifer / self.schematisation.recharge_rate)

        else:
            # travel_distance_shallow_aquifer = 
            if depth_point_contamination <= head:
                travel_distance_shallow_aquifer = depth_point_contamination - self.schematisation.bottom_shallow_aquifer
                if travel_distance_shallow_aquifer < 0:
                    travel_time_shallow_aquifer = np.array([0])
                else: 
                    travel_time_shallow_aquifer = np.array([(travel_distance_shallow_aquifer
                            * self.schematisation.porosity_shallow_aquifer / self.schematisation.recharge_rate)])
            else:
                travel_distance_shallow_aquifer = self.schematisation.thickness_shallow_aquifer - (self.schematisation.groundwater_level - head)
                travel_time_shallow_aquifer = ((travel_distance_shallow_aquifer)
                            * self.schematisation.porosity_shallow_aquifer / self.schematisation.recharge_rate)

        return travel_time_shallow_aquifer

    
    def _calculate_travel_time_target_aquifer_phreatic(self, 
                                                        fraction_flux = None, 
                                                        distance= None):
        ''' Use a fraction_flux array to calculate the travel time for all flowlines or specify
        a specific distance from the well to calculate the travel time '''

        if fraction_flux is None:
            '''point source calculation'''
            travel_time_target_aquifer = (self.schematisation.porosity_target_aquifer * self.schematisation.thickness_target_aquifer 
                                        / self.schematisation.recharge_rate
                                        * math.log(self.schematisation.well_discharge 
                                        / (self.schematisation.well_discharge - math.pi * self.schematisation.recharge_rate
                                         * distance ** 2 ) )) 
            travel_time_target_aquifer = np.array([travel_time_target_aquifer])

        elif distance is None:
            ''' diffuse calculation or regular TTD'''
            travel_time_target_aquifer = (self.schematisation.porosity_target_aquifer * self.schematisation.thickness_target_aquifer 
                                / self.schematisation.recharge_rate
                                * np.log(1 / (1 -fraction_flux)))
    
        return travel_time_target_aquifer

    def _check_init_confined(self):
        #AH_todo update these

        '''check the variables that we need for the individual aquifer types are not NONE aka set by the user'''
       
        required_variables = ["schematisation_type", #repeat for all
                              "thickness_vadose_zone",
                              "thickness_shallow_aquifer",
                              "thickness_target_aquifer",
                              "porosity_vadose_zone",
                              "porosity_shallow_aquifer",
                              "porosity_target_aquifer",
                              "well_discharge",
                              "recharge_rate",
                              "vertical_resistance_aquitard",
                              "moisture_content_vadose_zone",
                              "KD",
                              "groundwater_level",
                              "thickness_full_capillary_fringe",
                            #   self.test_variable,
                            ]
        self._check_required_variables(required_variables)

    def _calculate_travel_time_aquitard_semiconfined(self, distance,depth_point_contamination):
        ''' Calculate the travel time in shallow aquifer (aquitard in semiconfined case)
        using the the Peters (1985) solution( eq. 8.8 in Peters \cite{Peters1985})

        Equation A.12 in report BUT now implemented with n' (fraction of aquitard 
        contacted, to account for gaps in aquitard [-] = porosity of the shallow 
        aquifer (aquitard))
        '''          
        if depth_point_contamination is None:
            travel_distance_shallow_aquifer  = self.schematisation.thickness_shallow_aquifer
        elif depth_point_contamination > self.schematisation.bottom_shallow_aquifer:
            travel_distance_shallow_aquifer  = self.schematisation.thickness_shallow_aquifer

        else:
            travel_distance_shallow_aquifer  = depth_point_contamination - self.schematisation.bottom_shallow_aquifer 

        self.travel_time_shallow_aquifer = (self.schematisation.porosity_shallow_aquifer 
                                            * (2 * math.pi * self.schematisation.KD * self.schematisation.vertical_resistance_aquitard
                                            / (self.schematisation.well_discharge)
                                            * (travel_distance_shallow_aquifer 
                                            / besselk(0, distance
                                            / math.sqrt(self.schematisation.KD * self.schematisation.vertical_resistance_aquitard)))
                            ))
        if travel_distance_shallow_aquifer < 0:
            self.travel_time_shallow_aquifer =  np.array([0])

        return self.travel_time_shallow_aquifer


    def _calculate_travel_time_target_aquifer_semiconfined(self, distance):
        '''Calculate the travel time in zone 2 (aquifer with the production well)
        using the the Peters (1985) solution
        Equation A.13/A.14 in report

        '''
        # porosity_target_aquifer, AH #we keep this number for comparing to excel, but will change later to be the user defined porosities
        porosity_target_aquifer=0.32  
        thickness_target_aquifer=95

        self.travel_time_target_aquifer = (2 * math.pi * self.spreading_distance ** 2 / (self.schematisation.well_discharge)
                            * porosity_target_aquifer * thickness_target_aquifer
                            * (1.0872 * (distance / self.spreading_distance) ** 3
                                - 1.7689 * (distance /
                                            self.spreading_distance) ** 2
                                + 1.5842 * (distance / self.spreading_distance) - 0.2544)
                            )
        self.travel_time_target_aquifer[self.travel_time_target_aquifer < 0] = 0

        return self.travel_time_target_aquifer


    def _calculate_hydraulic_head(self, distance):
        '''Calculate the hydraulic head in meters above sea level for the 
        semiconfined case (AH confirm?)'''

        self.head = (-self.schematisation.well_discharge / (2 * math.pi * self.schematisation.KD)
                * besselk(0, distance / self.schematisation.spreading_distance)
                + self.schematisation.groundwater_level)

        return self.head


    def _calculate_flux_fraction(self,
                                radial_distance,
                                spreading_distance,
                                ):
        '''Calculates the fraction of the flux at distance x

        [-], Qr/Qw
        Qr = flux at distance x
        Qw = well flux '''

        flux_fraction = (radial_distance / spreading_distance
                        * besselk(1, radial_distance / spreading_distance))

        return flux_fraction

    def _create_output_dataframe(self, 
                                total_travel_time, 
                                travel_time_unsaturated, 
                                travel_time_shallow_aquifer, 
                                travel_time_target_aquifer,
                                distance, 
                                head, 
                                cumulative_fraction_abstracted_water,  
                                # flowline_discharge
                                ):
        # fraction flux ( array of 1,,23,....,99,99.5,99.9), row 7 in TTD phreatic, difference between cells divided by the well discharge
        # self.schematisation = schematisation
        column_names = ["total_travel_time", "travel_time_unsaturated", 
                        "travel_time_shallow_aquifer", "travel_time_target_aquifer",
                        "radial_distance", "head", "cumulative_fraction_abstracted_water", 
                        "flowline_discharge"
                        ]

        # if isinstance(cumulative_fraction_abstracted_water, float):
        # if isinstance(distance, float):
        #AH_here
        if len(distance) == 1:
            #point source
            flowline_discharge = cumulative_fraction_abstracted_water * self.schematisation.well_discharge

        else:
            #diffuse source
            flowline_discharge = (np.diff(np.insert(cumulative_fraction_abstracted_water,0,0., axis=0)))*self.schematisation.well_discharge

        data = [total_travel_time, 
                travel_time_unsaturated, 
                travel_time_shallow_aquifer, 
                travel_time_target_aquifer,
                distance, 
                head, 
                cumulative_fraction_abstracted_water,  
                flowline_discharge,
                ]
        df_output = pd.DataFrame (data = np.transpose(data), columns=column_names, dtype="object")
            
        return df_output
        

    def _export_to_df(self,
        df_output, 
        distance,
        total_travel_time, 
        travel_time_unsaturated, 
        travel_time_shallow_aquifer, 
        travel_time_target_aquifer,
        discharge_point_contamination=None,
        # input_concentration=self.schematisation.diffuse_input_concentration,
        ):
        """ Export to dataframe....

            Parameters
            ----------
            what_to_export: String
                options: 'all', 'omp_parameters', 'microbial_parameters'
            """
        what_to_export = self.schematisation.what_to_export
        #------------------------------
        # Make df_particle
        #------------------------------
        def fill_df_particle (df, 
                            distance, 
                            travel_time_unsaturated, 
                            travel_time_shallow_aquifer, 
                            travel_time_target_aquifer,
                            total_travel_time):
            
            df.loc[0] = [flowline_id, 
                        "surface", 
                        0, 
                        0, 
                        distance,
                        self.schematisation.model_width,
                         self.schematisation.ground_surface, 
                         None, 
                         self.schematisation.temperature, 
                         None, 
                         None, 
                         None, 
                         None, 
                         None,
                         None]

            df.loc[1] = [flowline_id, 
                         "vadose_zone",
                         travel_time_unsaturated,
                         travel_time_unsaturated,
                         distance,
                         self.schematisation.model_width,
                         self.schematisation.bottom_vadose_zone_at_boundary, # @MartinvdS should this be the thickness_vadose_zone_drawdown??
                         self.schematisation.redox_vadose_zone, 
                         self.schematisation.temperature,
                         self.schematisation.thickness_vadose_zone_at_boundary,
                         self.schematisation.porosity_vadose_zone,
                         self.schematisation.dissolved_organic_carbon_vadose_zone,
                         self.schematisation.pH_vadose_zone,
                         self.schematisation.fraction_organic_carbon_vadose_zone,
                         self.schematisation.solid_density_vadose_zone,
                         ]

            df.loc[2] = [flowline_id, 
                        "shallow_aquifer", 
                         travel_time_shallow_aquifer,
                         travel_time_unsaturated + travel_time_shallow_aquifer,
                         distance,
                         self.schematisation.model_width,
                         self.schematisation.bottom_shallow_aquifer,
                         self.schematisation.redox_shallow_aquifer, 
                         self.schematisation.temperature,
                         self.schematisation.thickness_shallow_aquifer, # @MartinvdS does this need to account for the vadose zone drawdown? 
                         self.schematisation.porosity_shallow_aquifer,
                         self.schematisation.dissolved_organic_carbon_shallow_aquifer,
                         self.schematisation.pH_shallow_aquifer,
                         self.schematisation.fraction_organic_carbon_shallow_aquifer,
                         self.schematisation.solid_density_shallow_aquifer, 
                         ]

            df.loc[3] = [flowline_id, 
                        "target_aquifer",  
                         travel_time_target_aquifer,
                         total_travel_time,
                         self.schematisation.diameter_borehole/2, #at the well
                         self.schematisation.model_width,
                         self.schematisation.bottom_target_aquifer,
                         self.schematisation.redox_target_aquifer, 
                         self.schematisation.temperature,
                         self.schematisation.thickness_target_aquifer,
                         self.schematisation.porosity_target_aquifer,
                         self.schematisation.dissolved_organic_carbon_target_aquifer,
                         self.schematisation.pH_target_aquifer,
                         self.schematisation.fraction_organic_carbon_target_aquifer,
                         self.schematisation.solid_density_target_aquifer,                         
                         ]
            return df

        df_particle = pd.DataFrame(columns=['flowline_id', 
                                                'zone', 
                                                'travel_time_zone', 
                                                'total_travel_time', 
                                                 'xcoord', #= radial_distcance,
                                                 'ycoord', #= the width of the cell .. default = 1 m 
                                                 'zcoord',
                                                 'redox_zone', 
                                                 'temperature', 
                                                 'thickness_layer', 
                                                 'porosity_layer', 
                                                 'dissolved_organic_carbon', 
                                                 'pH', 
                                                 'fraction_organic_carbon', 
                                                 'solid_density_layer'])

        df = df_particle.copy()

        for i in range(len(df_output)):
            flowline_id = i+1
            df = fill_df_particle (df= df, 
                                distance = distance[i], 
                                travel_time_unsaturated = travel_time_unsaturated[i], 
                                travel_time_shallow_aquifer=travel_time_shallow_aquifer[i], 
                                travel_time_target_aquifer = travel_time_target_aquifer[i],
                                total_travel_time= total_travel_time[i])

            df_particle = df_particle.append(df, ignore_index=True)
            df_particle['redox_zone'] = df_particle['redox_zone'].fillna('').astype(str)
            df_particle['flowline_id'] = df_particle['flowline_id'].astype(int)

        #------------------------------
        # Make df_flowline
        #------------------------------
        df_flowline = pd.DataFrame(columns=['flowline_id',
                                            'flowline_type', #diffuse_source or point_source
                                            'discharge',
                                            'particle_release_date',
                                            # 'input_concentration',
                                            'endpoint_id', ])
        df_flowline['discharge'] = df_output['flowline_discharge']
        df_flowline['flowline_type'] = 'diffuse_source'


        df_flowline['flowline_id'] =  df_flowline.index + 1
        df_flowline['particle_release_date'] = self.schematisation.particle_release_date


        # df_flowline['input_concentration'] = self.schematisation.input_concentration
        # #AH  Temporary solution here, @MartinK @MartinvdS, do we want to automatically 
        # make the dictionaries for modpath when running the analytical model? 
        # so far this is the only use for hte dicitonaries so they are made here
        self.schematisation.make_dictionary()
        endpoint_id = list(self.schematisation.well_parameters.items())[0][0]
        df_flowline['endpoint_id'] = endpoint_id 
        
        # AH which parameters for the 'microbial_parameters' option? @MartinvdS or @steven

        if what_to_export == 'all' or what_to_export== 'omp_parameters':

            df_flowline['well_discharge'] = self.schematisation.well_discharge
            df_flowline['recharge_rate'] = self.schematisation.recharge_rate
            df_flowline['vertical_resistance_aquitard'] = self.schematisation.vertical_resistance_aquitard
            df_flowline['KD'] = self.schematisation.KD
            df_flowline['thickness_full_capillary_fringe'] = self.schematisation.thickness_full_capillary_fringe
            # df_flowline['recharge_concentration'] = self.schematisation.recharge_concentration
            df_flowline['substance'] = self.schematisation.substance
            # df_flowline['input_concentration'] = input_concentration
            df_flowline['particle_release_date'] = self.schematisation.particle_release_date
            df_flowline['moisture_content_vadose_zone'] = self.schematisation.moisture_content_vadose_zone
            df_flowline['diameter_borehole'] = self.schematisation.diameter_borehole
            df_flowline['removal_function'] = self.schematisation.removal_function
            df_flowline['solid_density_vadose_zone'] = self.schematisation.solid_density_vadose_zone
            df_flowline['solid_density_shallow_aquifer'] = self.schematisation.solid_density_shallow_aquifer
            df_flowline['solid_density_target_aquifer'] = self.schematisation.solid_density_target_aquifer

        # if what_to_export == 'all':
        #     pass
            #AH come back to this and fill in with the final parameters of interest
            # df_flowline['borehole_diameter'] = self.schematisation.borehole_diameter
            # df_flowline['k_hor_aquifer'] = self.schematisation.k_hor_aquifer
            # df_flowline['vani_aquifer'] = self.schematisation.vani_aquifer
            # df_flowline['k_hor_confining'] = self.schematisation.k_hor_confining
            # df_flowline['vani_confining'] = self.schematisation.vani_confining
            # df_flowline['k_hor_gravelpack'] = self.schematisation.k_hor_gravelpack
            # df_flowline['vani_gravelpack'] = self.schematisation.vani_gravelpack
            # df_flowline['k_hor_clayseal'] = self.schematisation.k_hor_clayseal
            # df_flowline['vani_clayseal'] = self.schematisation.vani_clayseal
            # df_flowline['dz_well'] = self.schematisation.dz_well

        # self.df_flowline = df_flowline
        return df_flowline, df_particle

    def phreatic(self, 
                distance=None, 
                depth_point_contamination=None,
                cumulative_fraction_abstracted_water=None,
                ):
                
        '''
        Function to create array of travel time distributionss
        for distances from well for phreatic aquifer scenario

        Parameter - Input
        ---------
        well_discharge                           # [m3/day]
        spreading_distance               # [m]?, sqtr(K*D*c)
        vertical_resistance_aquitard   # [d], c_V
        porosity_vadose_zone           # [-]
        porosity_shallow_aquifer       # [-]
        porosity_target_aquifer        # [-]
        KD                             # [m2/d], permeability * tickness of vertical layer
        permeability                   # [m2]

        thickness_vadose_zone           # [m], thickness of unsaturated zone/layer
        thickness_shallow_aquifer       # [m], thickness of zone 1, if semiconfined aquitard,
                                        if phreatic aquifer, aquifer zone above top of well screens.
        thickness_target_aquifer        # [m], thickness of zone 2

        thickness_full_capillary_fringe #[m], cF

        recharge_rate                        # [m/d], recharge of well area
        moisture_content_vadose_zone           # [m3/m3], θ
        travel_time_H2O                 # [d],  travel time of water along flowline, in zone

        Output
        ------
        radial_distance_recharge        # [m], radial distance to well in order to recharge the well
        radial_distance                 # [m], radial distance to well (field),
                                        from site X within and any site on the groundwater divide
        '''
        self._check_init_phreatic()

        # travel time unsaturated now calculated in the HydrochemicalSchematisation class

        if distance is None:
            self.schematisation._calculate_travel_time_unsaturated_zone()
            self.radial_distance = self.schematisation.radial_distance
            fraction_flux=self.schematisation.fraction_flux
        else: 
            self.schematisation._calculate_travel_time_unsaturated_zone(distance=distance, depth_point_contamination=depth_point_contamination)
            self.radial_distance = distance
            fraction_flux=None
            
        self.travel_time_unsaturated = self.schematisation.travel_time_unsaturated
        self.head = self.schematisation.head

        self.travel_time_shallow_aquifer = self._calculate_travel_time_shallow_aquifer_phreatic(head=self.head, 
                                                                                                depth_point_contamination=depth_point_contamination)
        
        self.travel_time_target_aquifer = self._calculate_travel_time_target_aquifer_phreatic(fraction_flux=fraction_flux, distance=distance)

        self.total_travel_time = (self.travel_time_unsaturated + self.travel_time_shallow_aquifer
                            + self.travel_time_target_aquifer)
        
        if cumulative_fraction_abstracted_water is None:
            self.cumulative_fraction_abstracted_water = self.schematisation.fraction_flux 
        else:
            self.cumulative_fraction_abstracted_water = cumulative_fraction_abstracted_water

        self.df_output = self._create_output_dataframe(total_travel_time=self.total_travel_time, 
                    travel_time_unsaturated = self.travel_time_unsaturated, 
                    travel_time_shallow_aquifer=self.travel_time_shallow_aquifer, 
                    travel_time_target_aquifer=self.travel_time_target_aquifer,
                    distance=self.radial_distance,
                    head=self.head, 
                    cumulative_fraction_abstracted_water = self.cumulative_fraction_abstracted_water,  
                    )
            

        self.df_flowline, self.df_particle = self._export_to_df(df_output=self.df_output, 
                    distance=self.radial_distance,
                    total_travel_time=self.total_travel_time, 
                    travel_time_unsaturated = self.travel_time_unsaturated, 
                    travel_time_shallow_aquifer=self.travel_time_shallow_aquifer, 
                    travel_time_target_aquifer=self.travel_time_target_aquifer,
                    discharge_point_contamination = self.schematisation.discharge_point_contamination)

    def add_phreatic_point_sources(self, 
                distance=None, 
                depth_point_contamination=None,
                cumulative_fraction_abstracted_water=None,
                ):
                
        '''
        Function to create array of travel time distributionss
        for distances from well for phreatic aquifer scenario

        Parameter - Input
        ---------
        well_discharge                           # [m3/day]
        spreading_distance               # [m]?, sqtr(K*D*c)
        vertical_resistance_aquitard   # [d], c_V
        porosity_vadose_zone           # [-]
        porosity_shallow_aquifer       # [-]
        porosity_target_aquifer        # [-]
        KD                             # [m2/d], permeability * tickness of vertical layer
        permeability                   # [m2]

        thickness_vadose_zone           # [m], thickness of unsaturated zone/layer
        thickness_shallow_aquifer       # [m], thickness of zone 1, if semiconfined aquitard,
                                        if phreatic aquifer, aquifer zone above top of well screens.
        thickness_target_aquifer        # [m], thickness of zone 2

        thickness_full_capillary_fringe #[m], cF

        recharge_rate                        # [m/d], recharge of well area
        moisture_content_vadose_zone           # [m3/m3], θ
        travel_time_H2O                 # [d],  travel time of water along flowline, in zone

        Output
        ------
        radial_distance_recharge        # [m], radial distance to well in order to recharge the well
        radial_distance                 # [m], radial distance to well (field),
                                        from site X within and any site on the groundwater divide
        '''
        self._check_init_phreatic()

        # travel time unsaturated now calculated in the HydrochemicalSchematisation class

        if distance is None:
            self.schematisation._calculate_travel_time_unsaturated_zone()
            radial_distance = self.schematisation.radial_distance
            fraction_flux=self.schematisation.fraction_flux
        else: 
            self.schematisation._calculate_travel_time_unsaturated_zone(distance=distance, depth_point_contamination=depth_point_contamination)
            radial_distance = distance
            fraction_flux=None
            
        travel_time_unsaturated = self.schematisation.travel_time_unsaturated
        head = self.schematisation.head

        travel_time_shallow_aquifer = self._calculate_travel_time_shallow_aquifer_phreatic(head=head, 
                                                                                        depth_point_contamination=depth_point_contamination)
        
        travel_time_target_aquifer = self._calculate_travel_time_target_aquifer_phreatic(fraction_flux=fraction_flux,
                                                                                        distance=distance)

        total_travel_time = (travel_time_unsaturated + travel_time_shallow_aquifer
                            + travel_time_target_aquifer)
        
        if cumulative_fraction_abstracted_water is None:
            cumulative_fraction_abstracted_water = self.schematisation.fraction_flux 
        else:
            cumulative_fraction_abstracted_water = cumulative_fraction_abstracted_water

        df_output = self._create_output_dataframe(total_travel_time=total_travel_time, 
                    travel_time_unsaturated = travel_time_unsaturated, 
                    travel_time_shallow_aquifer=travel_time_shallow_aquifer, 
                    travel_time_target_aquifer=travel_time_target_aquifer,
                    distance=radial_distance,
                    head=head, 
                    cumulative_fraction_abstracted_water = cumulative_fraction_abstracted_water,  
                    )
            
        df_flowline, df_particle=self._export_to_df(df_output=df_output, 
                    distance=radial_distance,
                    total_travel_time=total_travel_time, 
                    travel_time_unsaturated = travel_time_unsaturated, 
                    travel_time_shallow_aquifer=travel_time_shallow_aquifer, 
                    travel_time_target_aquifer=travel_time_target_aquifer,
                    discharge_point_contamination = self.schematisation.discharge_point_contamination)
        
        return df_flowline, df_particle

    def semiconfined(self,
                    distance = None, 
                    depth_point_contamination=None,  ):
        # self._check_init_confined() #AH_todo this is not implemented correctly fix!

        '''
        "The transit time distribution (TTD), also called hydrological response
        curve (HRC), is defined as the cumulative frequency distribution of
        travel times, in our case from land surface or open water course to
        a well (field)." cite(Stuyfzand2020)

        Function to create array of travel time distributionss
        for distances from well.

        Parameter - Input
        ---------
        well_discharge                           # [m3/d]
        spreading_distance               # [m]?, sqtr(K*D*c)
        vertical_resistance_aquitard   # [d], c_V
        porosity_shallow_aquifer                 # [-]
        KD                             # [m2/d], permeability * tickness of vertical layer
        permeability                   # [m2]

        thickness_vadose_zone      # [m], thickness of unsaturated zone/layer
        thickness_shallow_aquifer                 # [m], thickness of zone 1, if semiconfined aquitard,
                                        if phreatic aquifer, aquifer zone above top of well screens.
        thickness_target_aquifer                 # [m], thickness of zone 2

        thickness_full_capillary_fringe #[m], cF

        recharge_rate                        # [m/d], recharge_rate of well area
        moisture_content_vadose_zone           # [m3/m3], θ
        travel_time_H2O                 # [d],  travel time of water along flowline, in zone

        Output
        ------
        radial_distance_recharge        # [m], radial distance to well in order to recharge_rate the well
        radial_distance                 # [m], radial distance to well (field),
                                        from site X within and any site on the groundwater divide

        '''
        # travel time unsaturated now calculated in the HydrochemicalSchematisation class
        # self.schematisation._calculate_travel_time_unsaturated_zone()
        # self.travel_time_unsaturated = self.schematisation.travel_time_unsaturated
        # self.radial_distance = self.schematisation.radial_distance
        # self.spreading_distance = self.schematisation.spreading_distance

        if distance is None:
            self.schematisation._calculate_travel_time_unsaturated_zone()
            self.radial_distance = self.schematisation.radial_distance
        else: 
            self.schematisation._calculate_travel_time_unsaturated_zone(distance=distance,
                                                                        depth_point_contamination = depth_point_contamination
                                                                        )
            self.radial_distance = distance
        
        self.travel_time_unsaturated = self.schematisation.travel_time_unsaturated
        self.spreading_distance = self.schematisation.spreading_distance

        # travel time in semiconfined is one value, make it array by repeating the value
        # self.travel_time_unsaturated = [self.travel_time_unsaturated] * (len(self.radial_distance))

        #left off here, need to make sure this is up to date? distances????
        self.travel_time_shallow_aquifer = self._calculate_travel_time_aquitard_semiconfined(distance=self.radial_distance,
                                                                                                depth_point_contamination=depth_point_contamination)

        self.travel_time_target_aquifer = self._calculate_travel_time_target_aquifer_semiconfined(distance=self.radial_distance) 

        self.total_travel_time = self.travel_time_unsaturated + self.travel_time_shallow_aquifer + self.travel_time_target_aquifer

        self.head = self._calculate_hydraulic_head(distance=self.radial_distance)

        # percent of the max head? AH
        self.head_minus_max = 100*(self.head - self.schematisation.groundwater_level) / (self.head[0] - self.schematisation.groundwater_level)

        self.flux_fraction = self._calculate_flux_fraction(radial_distance=self.radial_distance,
                                spreading_distance = self.spreading_distance,
                                )

        ''' Calculate the cumulative_fraction_abstracted_water
        1.1369 comes form pg. 52 in report, describes cutting off the 
        recharge_rate distance at 3 labda, need to increase the fraction abstracted from ~87%
        to 99.9% so multiply by 1.1369 to get to that
        Equation A.16 in report'''
        # AH, may want to change this, to eg. 6 labda or something else, adjust this number
        self.cumulative_fraction_abstracted_water = 1.1369 * (1 - self.flux_fraction)


        self.df_output = self._create_output_dataframe(total_travel_time=self.total_travel_time, 
                    travel_time_unsaturated = self.travel_time_unsaturated, 
                    travel_time_shallow_aquifer=self.travel_time_shallow_aquifer, 
                    travel_time_target_aquifer=self.travel_time_target_aquifer,
                    distance=self.radial_distance,
                    head=self.head, 
                    cumulative_fraction_abstracted_water = self.cumulative_fraction_abstracted_water,  
                    )

        self.df_flowline, self.df_particle=self._export_to_df(df_output=self.df_output, 
                    distance=self.radial_distance,
                    total_travel_time=self.total_travel_time, 
                    travel_time_unsaturated = self.travel_time_unsaturated, 
                    travel_time_shallow_aquifer=self.travel_time_shallow_aquifer, 
                    travel_time_target_aquifer=self.travel_time_target_aquifer,
                    discharge_point_contamination = self.schematisation.discharge_point_contamination)

    def add_semiconfined_point_sources(self, 
                                    distance = None, 
                                    depth_point_contamination=None,  ):
        ''' Same as the semiconfined except that the attributes are not updated, 
        this ensures that the attributes of the well class remain as the ones for 
        the flowlines from the well, not the point sources'''

        if distance is None:
            self.schematisation._calculate_travel_time_unsaturated_zone()
            radial_distance = self.schematisation.radial_distance
        else: 
            self.schematisation._calculate_travel_time_unsaturated_zone(distance=distance,
                                                                        depth_point_contamination = depth_point_contamination
                                                                        )
            radial_distance = distance
        
        travel_time_unsaturated = self.schematisation.travel_time_unsaturated
        spreading_distance = self.schematisation.spreading_distance

        # travel time in semiconfined is one value, make it array by repeating the value
        # self.travel_time_unsaturated = [self.travel_time_unsaturated] * (len(self.radial_distance))

        #left off here, need to make sure this is up to date? distances????
        travel_time_shallow_aquifer = self._calculate_travel_time_aquitard_semiconfined(distance=radial_distance,
                                                                                    depth_point_contamination=depth_point_contamination)

        travel_time_target_aquifer = self._calculate_travel_time_target_aquifer_semiconfined(distance=radial_distance) 

        total_travel_time = travel_time_unsaturated + travel_time_shallow_aquifer + travel_time_target_aquifer

        head = self._calculate_hydraulic_head(distance=radial_distance)

        # percent of the max head? AH
        head_minus_max = 100*(head - self.schematisation.groundwater_level) / (head[0] - self.schematisation.groundwater_level)

        flux_fraction = self._calculate_flux_fraction(radial_distance=radial_distance,
                                spreading_distance = spreading_distance,)

        ''' Calculate the cumulative_fraction_abstracted_water
        1.1369 comes form pg. 52 in report, describes cutting off the 
        recharge_rate distance at 3 labda, need to increase the fraction abstracted from ~87%
        to 99.9% so multiply by 1.1369 to get to that
        Equation A.16 in report'''
        # AH, may want to change this, to eg. 6 labda or something else, adjust this number
        cumulative_fraction_abstracted_water = 1.1369 * (1 - flux_fraction)


        df_output = self._create_output_dataframe(total_travel_time=total_travel_time, 
                    travel_time_unsaturated = travel_time_unsaturated, 
                    travel_time_shallow_aquifer=travel_time_shallow_aquifer, 
                    travel_time_target_aquifer=travel_time_target_aquifer,
                    distance=radial_distance,
                    head=head, 
                    cumulative_fraction_abstracted_water = cumulative_fraction_abstracted_water,  
                    )

        df_flowline, df_particle=self._export_to_df(df_output=df_output, 
                    distance=radial_distance,
                    total_travel_time=total_travel_time, 
                    travel_time_unsaturated = travel_time_unsaturated, 
                    travel_time_shallow_aquifer=travel_time_shallow_aquifer, 
                    travel_time_target_aquifer=travel_time_target_aquifer,
                    discharge_point_contamination = self.schematisation.discharge_point_contamination)
        
        return df_flowline, df_particle


    def plot_travel_time_versus_radial_distance(self,
                                                xlim=[0, 4000],
                                                ylim=[1, 5000]):
        ''' Plot the travel time versus the radial distance '''

        fig = plt.figure(figsize=[10, 5])
        plt.plot(self.radial_distance, self.total_travel_time, 'r', label=self.schematisation.schematisation_type)
        
        plt.xlim(xlim)
        plt.ylim(ylim) 
        plt.yscale('log')
        plt.xlabel('Radial distance to well(s) (m)')
        plt.ylabel('Total travel time (days)')
        plt.title('Aquifer type: ' + self.schematisation.schematisation_type)
        plt.grid()
        plt.legend()
        plt.savefig('travel_time_versus_radial_distance_'+self.schematisation.schematisation_type+'.png', dpi=300, bbox_inches='tight')  # save_results_to + '/


    def plot_travel_time_versus_cumulative_abstracted_water(self,
                                                            xlim,
                                                            ylim=[1, 5000]):
        
        ''' Plot the travel time versus the cumulative abstracted water '''

        fig = plt.figure(figsize=[10, 5])
        plt.plot(self.cumulative_fraction_abstracted_water, self.total_travel_time, 'r', label=self.schematisation.schematisation_type)
        plt.xlim(xlim)  # [0.01,10000])
        plt.ylim(ylim) 
        plt.yscale('log')
        plt.xlabel('Cumulative fraction of abstracted water')
        plt.ylabel('Total travel time (days)')
        plt.title('Aquifer type: ' + self.schematisation.schematisation_type)
        plt.grid()
        plt.legend()
        plt.savefig('travel_time_versus_cumulative_fraction_abstracted_water_'+self.schematisation.schematisation_type+'.png', dpi=300, bbox_inches='tight')  # save_results_to + '/

