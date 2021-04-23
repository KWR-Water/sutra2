#%% ----------------------------------------------------------------------------
# A. Hockin, March 2021
# KWR BO 402045-247
# ZZS verwijdering bodempassage
# AquaPriori - Transport Model
# With Martin Korevaar, Martin vd Schans, Steven Ros
#
# Based on Stuyfzand, P. J. (2020). Predicting organic micropollutant behavior 
#               for 4 public supply well field types, with TRANSATOMIC Lite+
#               (Vol. 2). Nieuwegein, Netherlands.
# ------------------------------------------------------------------------------

#### CHANGE LOG ####
# things which must be checked indicated in comm ents with AH

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

path = os.getcwd()  # path of working directory

#%% 

#%% 

# ------------------------------------------------------------------------------
# Questions
# ------------------------------------------------------------------------------

# 1. How to structure the input parameters?

# ------------------------------------------------------------------------------
# Phreatic and Semi-Confined Aquifer Functions
# ------------------------------------------------------------------------------

"""
Parameters
----------
schematisation: string
    'phreatic', 'semi-confined', 'riverbankfiltration', 'basinfiltration'
removal_function: string
    'omp' -> get parameters for OMP
    'microbiology' -> get parameters for microbiology
"""


#%%
# AH in general need more guidance on how to use this...

# AH does a long list of all the paramters come here?
class HydroChemicalSchematisation:

    """ Converts input parameters of AquaPriori GUI to a complete parameterisation."""

    #everything initialized to None to start...
    def __init__(self,
                 schematisation=None,
                 thickness_vadose_zone=None,
                 thickness_shallow_aquifer=None,
                 thickness_target_aquifer=None,
                 porosity_vadose_zone=None,
                 porosity_shallow_aquifer=None,
                 porosity_target_aquifer=None,
                 well_discharge_m3hour=None,
                 recharge_rate=None,
                 vertical_resistance_aquitard=None,
                 soil_moisture_content_vadose_zone=None,
                 KD=None,  # AH premeability*thicnkeness.
                 thickness_full_capillary_fringe=None,
                 dissolved_organic_carbon_vadose_zone=None,
                 dissolved_organic_carbon_shallow_aquifer=None,
                 dissolved_organic_carbon_target_aquifer=None,
                 recharge_concentration=None,
                 substance=None,
                 redox_vadose_zone=None,
                 redox_shallow_aquifer=None,
                 redox_target_aquifer=None,
                 input_concentration=None,
                 particle_release_date=None,
                 
                 fraction_organic_carbon_vadose_zone=None,
                 fraction_organic_carbon_shallow_aquifer=None,
                 fraction_organic_carbon_target_aquifer=None,
                 pH_vadose_zone=None,
                 pH_shallow_aquifer=None,
                 pH_target_aquifer=None,):
        #DEFAULT VALUES
        # AH can/should we (user) be able to alter these?
        # aka do we have to include these as arguements for the class?
        # AH do I add self. in from of each of these?
        self.well_screen_diameter = .75    # m
        self.borehole_diameter = .75        # m
        self.temperature = 11         # celcius
        self.k_hor_aquifer = 10             # m/d
        self.vani_aquifer = 1.             # -
        self.k_hor_confining = .001         # m/d
        self.vani_confining = 1.          # -
        self.k_hor_gravelpack = 100         # m/d
        self.vani_gravelpack = 1.           # -
        self.k_hor_clayseal = .001          # m/d
        self.vani_clayseal = 1.             # -
        self.density_aquifer = 2650.        # kg/m3
        self.removal_function = 'omp'
        self.dz_well = 0.5

        # non-default parameters
        self.schematisation = schematisation

        self.thickness_vadose_zone = thickness_vadose_zone
        self.thickness_shallow_aquifer = thickness_shallow_aquifer
        self.thickness_target_aquifer = thickness_target_aquifer

        self.porosity_vadose_zone = porosity_vadose_zone
        self.porosity_shallow_aquifer = porosity_shallow_aquifer
        self.porosity_target_aquifer = porosity_target_aquifer

        self.dissolved_organic_carbon_vadose_zone = dissolved_organic_carbon_vadose_zone
        self.dissolved_organic_carbon_shallow_aquifer = dissolved_organic_carbon_shallow_aquifer
        self.dissolved_organic_carbon_target_aquifer = dissolved_organic_carbon_target_aquifer

        # 1 = (sub)oxic; 2 = anoxic; 3 = deeply anoxic
        self.redox_vadose_zone = redox_vadose_zone
        self.redox_shallow_aquifer = redox_shallow_aquifer
        self.redox_target_aquifer = redox_target_aquifer

        self.fraction_organic_carbon_vadose_zone = fraction_organic_carbon_vadose_zone
        self.fraction_organic_carbon_shallow_aquifer = fraction_organic_carbon_shallow_aquifer
        self.fraction_organic_carbon_target_aquifer = fraction_organic_carbon_target_aquifer

        self.pH_vadose_zone = pH_vadose_zone
        self.pH_shallow_aquifer = pH_shallow_aquifer
        self.pH_target_aquifer = pH_target_aquifer


        self.vertical_resistance_aquitard = vertical_resistance_aquitard
        self.soil_moisture_content_vadose_zone = soil_moisture_content_vadose_zone

        self.recharge_rate = recharge_rate

        self.well_discharge_m3hour = well_discharge_m3hour

        # AH input these seperately? # [m2/d], permeability * thickness of vertical layer
        self.KD = KD
        self.groundwater_level = 0-self.thickness_vadose_zone
        self.thickness_full_capillary_fringe = thickness_full_capillary_fringe

        self.recharge_concentration = recharge_concentration

        self.substance = substance

        # AH are the input concentration and release dates single values? or do they vary over time/distance?
        self.input_concentration=input_concentration
        self.particle_release_date=particle_release_date

    def make_dictionary(self,):
    # make dictionaries of each grouping of parameters
    
        recharge = {
            'recharge1': {
                'recharge_rate': self.recharge_rate
                # 'compound': ....
                # 'concentration': ....
            }}

        geo_parameters = {
            'geolayer1': {
                'vadose': True,  # not part of modflow
                'top': 0,  # flow parameters
                'bot': 0 - self.thickness_vadose_zone,
                # 'K_hor': self.k_hor_aquifer, # AH what goes here?
                # 'K_vert': self.K_vert, # AH what goes here?
                'porosity': self.porosity_vadose_zone,
                'redoxzone': self.redox_vadose_zone,  # concentration parameters
                'dissolved_organic_carbon': self.dissolved_organic_carbon_vadose_zone,
                'soil_moisture_content_vadose_zone': self.soil_moisture_content_vadose_zone,
                'thickness_full_capillary_fringe': self.thickness_full_capillary_fringe,
            },
            'geolayer2': {
                'vadose': False,  # not part of modflow
                'top': 0 - self.thickness_vadose_zone,  # flow parameters
                'bot': 0 - self.thickness_vadose_zone - self.thickness_shallow_aquifer,
                'K_hor': self.k_hor_aquifer,
                # 'K_vert': self.K_vert, # AH what goes here?
                'porosity': self.porosity_shallow_aquifer,
                'redoxzone': self.redox_shallow_aquifer,  # concentration parameters
                'dissolved_organic_carbon': self.dissolved_organic_carbon_shallow_aquifer,
                'vertical_resistance_aquitard': self.vertical_resistance_aquitard, # AH with the 2nd layer or 3rd layer?
            },
            'geolayer3': {
                'vadose': False,  # not part of modflow
                'top': 0 - self.thickness_vadose_zone - self.thickness_shallow_aquifer,  # flow parameters
                'bot': (0 - self.thickness_vadose_zone
                        - self.thickness_shallow_aquifer
                        - self.thickness_target_aquifer),
                'K_hor': self.k_hor_aquifer,
                # 'K_vert': self.K_vert, # AH what goes here?
                'porosity': self.porosity_target_aquifer,
                'redoxzone': self.redox_target_aquifer,  # concentration parameters
                'dissolved_organic_carbon': self.dissolved_organic_carbon_target_aquifer,
                'KD': self.KD,
            }
        }

        well_parameters = {
            'well1': {
                'top': 0 - self.thickness_vadose_zone - self.thickness_shallow_aquifer,
                'bot': (0 - self.thickness_vadose_zone - self.thickness_shallow_aquifer
                        - self.thickness_target_aquifer),
                'discharge': - self.well_discharge_m3hour / 24, # AH should this be '* 24' to make  m3/day?
                # 'screenradius': self.well_screen_diameter / 2,
                'dz_well': self.dz_well,
            },
            'well2': {
                'IBOUND': 0,  # 0=inactive, -1 -> fixed head
                'head': 0,
                'top': 0,
                # 'screenradius': self.well_screen_diameter / 2,
                'dz_well': self.dz_well,
                'bot': 0 - self.thickness_vadose_zone - self.thickness_shallow_aquifer,
            },
        }

        # AH ask Martin vdS what to put in these dicts (not used in the analytical functions)
        # backfill_parameters = {
        #     'backfill1': {
        #         'boreholeradius': self.borehole_diameter / 2,
        #         # 'K_hor': self.k_hor_clayseal # AH which clayseal? or gravelpack?
        #     }
        # }

        # drain_parameters = {
        #     'depth':
        #     'discharge':
        #     }

        # AH what goes here?
        # surfacewater_parameters = {
        # 'water_level': ....
        #     }

        # AH how do we want the dictionary made?
        # it is not very convenient for the analytical well models to have in one dictionary, 
        # but outputting many dictionaries is also inconvenient
        hydrochemical_dictionary = {'well_parameters': well_parameters,
                                    'geo_parameters': geo_parameters,
                                    'recharge': recharge, }
        return hydrochemical_dictionary


#%%
class AnalyticalWell():
    # AH do we pass the HydroChemicalSchematisation class to this class?
    # or the individual dictionaries??
    """ Compute travel time distribution using analytical well functions."""


    def __init__(self, schematisation_instance):
        """ 'unpack/parse' all the variables from the hydrogeochemical schematizization """
        # AH can i pass the instance of the schematisation? Then I can use the
        # instance.parameter to get all the values I want...
        #
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
        # AH can I replace this whole thing, by inheriting the HydroChemicalSchematisation class?

        self.schematisation = schematisation_instance.schematisation
        self.thickness_vadose_zone = schematisation_instance.thickness_vadose_zone
        self.thickness_shallow_aquifer = schematisation_instance.thickness_shallow_aquifer
        self.thickness_target_aquifer = schematisation_instance.thickness_target_aquifer
        self.porosity_vadose_zone = schematisation_instance.porosity_vadose_zone
        self.porosity_shallow_aquifer = schematisation_instance.porosity_shallow_aquifer
        self.porosity_target_aquifer = schematisation_instance.porosity_target_aquifer
        self.dissolved_organic_carbon_vadose_zone = schematisation_instance.dissolved_organic_carbon_vadose_zone
        self.dissolved_organic_carbon_shallow_aquifer = schematisation_instance.dissolved_organic_carbon_shallow_aquifer
        self.dissolved_organic_carbon_target_aquifer = schematisation_instance.dissolved_organic_carbon_target_aquifer
        self.redox_vadose_zone = schematisation_instance.redox_vadose_zone
        self.redox_shallow_aquifer = schematisation_instance.redox_shallow_aquifer
        self.redox_target_aquifer = schematisation_instance.redox_target_aquifer
        self.vertical_resistance_aquitard = schematisation_instance.vertical_resistance_aquitard
        self.soil_moisture_content_vadose_zone = schematisation_instance.soil_moisture_content_vadose_zone
        self.recharge_rate = schematisation_instance.recharge_rate
        self.well_discharge_m3hour = schematisation_instance.well_discharge_m3hour
        self.KD = schematisation_instance.KD
        self.groundwater_level=schematisation_instance.groundwater_level
        self.thickness_full_capillary_fringe = schematisation_instance.thickness_full_capillary_fringe
        self.recharge_concentration = schematisation_instance.recharge_concentration
        self.substance = schematisation_instance.substance

        self.fraction_organic_carbon_vadose_zone= schematisation_instance.fraction_organic_carbon_vadose_zone
        self.fraction_organic_carbon_shallow_aquifer= schematisation_instance.fraction_organic_carbon_shallow_aquifer
        self.fraction_organic_carbon_target_aquifer= schematisation_instance.fraction_organic_carbon_target_aquifer
        self.pH_vadose_zone= schematisation_instance.pH_vadose_zone
        self.pH_shallow_aquifer= schematisation_instance.pH_shallow_aquifer
        self.pH_target_aquifer= schematisation_instance.pH_target_aquifer


        self.particle_release_date = schematisation_instance.particle_release_date
        self.input_concentration = schematisation_instance.input_concentration

        # get the default parameters
        self.well_screen_diameter = schematisation_instance.well_screen_diameter
        self.borehole_diameter = schematisation_instance.borehole_diameter
        self.temperature = schematisation_instance.temperature
        self.k_hor_aquifer = schematisation_instance.k_hor_aquifer
        self.vani_aquifer = schematisation_instance.vani_aquifer
        self.k_hor_confining = schematisation_instance.k_hor_confining
        self.vani_confining = schematisation_instance.vani_confining
        self.k_hor_gravelpack = schematisation_instance.k_hor_gravelpack
        self.vani_gravelpack = schematisation_instance.vani_gravelpack
        self.k_hor_clayseal = schematisation_instance.k_hor_clayseal
        self.vani_clayseal = schematisation_instance.vani_clayseal
        self.density_aquifer = schematisation_instance.density_aquifer
        self.removal_function = schematisation_instance.removal_function
        self.dz_well = schematisation_instance.dz_well

        # left off here, need to figure out how to extract these values for each of the dictionaries...
        # d = input_schematisation_dict['geo_parameters']['geolayer1']
        #     for key,val in d.items():
        #     exec(key + '=val')

    def _check_init_phreatic(self):
        # check the variables that we need for the individual aquifer types are not NONE aka set by the user

        required_variables = [self.schematisation,
                              self.thickness_vadose_zone,
                              self.thickness_shallow_aquifer,
                              self.thickness_target_aquifer,
                              self.porosity_vadose_zone,
                              self.porosity_shallow_aquifer,
                              self.porosity_target_aquifer,
                              self.well_discharge_m3hour,
                              self.recharge_rate,
                              self.vertical_resistance_aquitard,
                              self.soil_moisture_content_vadose_zone,
                              self.KD,
                              self.groundwater_level,
                              self.thickness_full_capillary_fringe,
                            #   self.test_variable,
                            ]

        #AH how to check this and raise error correctly??
        for req_var in required_variables:
            if req_var is None:
                raise KeyError('Error, somehow print the name of the variable and the value here!!')

    def _check_init_confined(self):
        # check the variables that we need for the individual aquifer types are not NONE aka set by the user

        required_variables = [self.schematisation,
                              self.thickness_vadose_zone,
                              self.thickness_shallow_aquifer,
                              self.thickness_target_aquifer,
                              self.porosity_vadose_zone,
                              self.porosity_shallow_aquifer,
                              self.porosity_target_aquifer,
                              self.well_discharge_m3hour,
                              self.recharge_rate,
                              self.vertical_resistance_aquitard,
                              self.soil_moisture_content_vadose_zone,
                              self.KD,
                              self.groundwater_level,
                              self.thickness_full_capillary_fringe,
                            #   self.test_variable,
                            ]
        #AH how to check this and raise error correctly??
        for req_var in required_variables:
            if req_var is None:
                print(req_var, ', !!somehow print the name of the variable and the value here!!')

    def _create_radial_distance_array(self):

        ''' Calculate the radial distance from a well,
        specifying the aquifer type (semi-confined, phreatic)
        '''
        self.percent_flux = np.array([0.001, 0.01, 0.1, 0.5])
        self.percent_flux = np.append(self.percent_flux, np.arange(1, 100, 1))
        self.percent_flux = np.append(self.percent_flux, [99.5, 99.99])

        self.radial_distance = self.radial_distance_recharge * \
            np.sqrt(self.percent_flux / 100)

        if self.schematisation == 'semi-confined':

            self.radial_distance = np.append(self.radial_distance,
                                    [(self.radial_distance[-1] + ((self.spreading_length * 3) - self.radial_distance[-1]) / 3),
                                    (self.radial_distance[-1] + 2 *
                                    ((self.spreading_length * 3) - self.radial_distance[-1]) / 3),
                                        (self.spreading_length * 3)])

        return self.radial_distance, self.percent_flux


    def _calculate_travel_time_unsaturated_zone(self):
        '''Equation A.11 in report '''

        # AH
        # (groundwater_level +thickness_vadose_zone) -
        # > before was the top of the unsaturated zone, 
        # confirm if this is the same?
        # aka is the top of the unsaturated zone always
        # the thickness_vadose_zone + groundwater_level
        self.travel_time_unsaturated = ((((self.groundwater_level + self.thickness_vadose_zone)
                                    - self.groundwater_level
                                    - self.thickness_full_capillary_fringe)
                                    * self.soil_moisture_content_vadose_zone
                                    + self.thickness_full_capillary_fringe
                                    * self.porosity_vadose_zone)
                                / self.recharge_rate)

        return self.travel_time_unsaturated


    def _calculate_travel_time_aquitard_zone1(self):
        ''' Calculate the travel time in zone 1 (aquitard in semi-confined case)
        using the the Peters (1985) solution
        Equation A.12 in report
        '''                                   

        self.travel_time_shallow_aquifer = (2 * math.pi * self.KD * self.vertical_resistance_aquitard
                            / (self.well_discharge_m3hour * 24)
                            * (self.thickness_shallow_aquifer
                                / besselk(0, self.radial_distance
                                        / math.sqrt(self.KD * self.vertical_resistance_aquitard)))
                            / 365.25)

        return self.travel_time_shallow_aquifer


    def _calculate_travel_time_target_aquifer(self):
        '''Calculate the travel time in zone 2 (aquifer with the production well)
        using the the Peters (1985) solution
        Equation A.13 in report

        '''
        # porosity_target_aquifer, AH #we keep this number for camparing to excel, but will change later to be the user defined porosities
        porosity_target_aquifer=0.32  
        thickness_target_aquifer=95

        self.travel_time_target_aquifer = (2 * math.pi * self.spreading_length ** 2 / (self.well_discharge_m3hour * 24)
                            * porosity_target_aquifer * thickness_target_aquifer
                            * (1.0872 * (self.radial_distance / self.spreading_length) ** 3
                                - 1.7689 * (self.radial_distance /
                                            self.spreading_length) ** 2
                                + 1.5842 * (self.radial_distance / self.spreading_length) - 0.2544)
                            / 365.25)

        self.travel_time_target_aquifer[self.travel_time_target_aquifer < 0] = 0

        return self.travel_time_target_aquifer


    def _calculuate_hydraulic_head(self):
        '''Calculate the hydraulic head in meters above sea level '''

        self.head = (-self.well_discharge_m3hour * 24 / (2 * math.pi * self.KD)
                * besselk(0, self.radial_distance / self.spreading_length)
                + self.groundwater_level)

        return self.head


    def _calculate_flux_fraction(self):
        '''Calculates the fraction of the flux at distance x
        [-], Qr/Qw
        Qr = flux at distance x
        Qw = well flux '''

        self.flux_fraction = (self.radial_distance / self.spreading_length
                        * besselk(1, self.radial_distance / self.spreading_length))

        return self.flux_fraction

    def _create_output_dataframe(self):
        # AH, confirm this is how to calculate the discharge of each flowline
        self.flowline_discharge = (np.diff( np.insert(self.cumulative_percent_abstracted_water,0,0., axis=0))/100)*self.well_discharge_m3hour

        column_names = ["total_travel_time", "travel_time_unsaturated", 
                        "travel_time_shallow_aquifer", "travel_time_target_aquifer",
                        "radial_distance", "head", "cumulative_percent_abstracted_water", 
                        "flowline_discharge"
                        ]

        # AH check the well_discharge_m3hour calculations for the streamline... Pr*Q? 
        self.data = [self.total_travel_time, self.travel_time_unsaturated, 
                    self.travel_time_shallow_aquifer, self.travel_time_target_aquifer,
                    self.radial_distance, self.head, self.cumulative_percent_abstracted_water, 
                    self.flowline_discharge
                    ]

        self.df_output = pd.DataFrame (data = np.transpose(self.data), columns=column_names)
        
        return self.df_output

    def phreatic(self):
        self._check_init_phreatic()
        
        # def travel_time_distribution_phreatic(self):

        '''
        Function to create array of travel time distributionss
        for distances from well for phreatic aquifer scenario

        Parameter - Input
        ---------
        well_discharge_m3hour                           # [m3/hr]
        spreading_length               # [m]?, sqtr(K*D*c)
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

        recharge_rate                        # [m/yr], recharge of well area
        soil_moisture_content_vadose_zone           # [m3/m3], θ
        travel_time_H2O                 # [d],  travel time of water along flowline, in zone

        Output
        ------
        radial_distance_recharge        # [m], radial distance to well in order to recharge the well
        radial_distance                 # [m], radial distance to well (field),
                                        from site X within and any site on the groundwater divide
        '''

        self.spreading_length = math.sqrt(self.vertical_resistance_aquitard * self.KD)

        self.radial_distance_recharge = (math.sqrt(self.well_discharge_m3hour
                                                    / (math.pi * self.recharge_rate / (365.25 * 24))))

        self.radial_distance, self.percent_flux = self._create_radial_distance_array()

        self.head = (self.groundwater_level - self.well_discharge_m3hour * 24 / (2 * math.pi * self.KD) 
                * np.log(self.radial_distance_recharge / self.radial_distance))

        # AH
        # (groundwater_level +thickness_vadose_zone) -
        # > before was the top of the unsaturated zone, 
        # confirm if this is the same?
        # aka is the top of the unsaturated zone always
        # the thickness_vadose_zone + groundwater_level
        self.depth_bottom_vadose_aquifer = self.thickness_vadose_zone #no drawdown
        self.depth_bottom_shallow_aquifer = self.thickness_vadose_zone + self.thickness_shallow_aquifer
        self.depth_bottom_target_aquifer = self.thickness_vadose_zone + self.thickness_shallow_aquifer + self.thickness_target_aquifer

        self.thickness_vadose_zone = (self.groundwater_level + self.thickness_vadose_zone) - self.head

        self.travel_time_unsaturated = self._calculate_travel_time_unsaturated_zone()

        self.travel_time_shallow_aquifer = ((self.thickness_shallow_aquifer - (self.groundwater_level - self.head))
                            * self.porosity_shallow_aquifer / self.recharge_rate)

        self.travel_time_target_aquifer = (self.porosity_target_aquifer * self.thickness_target_aquifer / self.recharge_rate
                            * np.log(1 / (1 - 0.01 * self.percent_flux)))

        self.total_travel_time = (self.travel_time_unsaturated + self.travel_time_shallow_aquifer
                            + self.travel_time_target_aquifer)
        
        self.cumulative_percent_abstracted_water = self.percent_flux

        self.df_output = self._create_output_dataframe()

        return self.df_output

    def semiconfined(self):
        self._check_init_confined

        '''
        "The transit time distribution (TTD), also called hydrological response
        curve (HRC), is defined as the cumulative frequency distribution of
        travel times, in our case from land surface or open water course to
        a well (field)." cite(Stuyfzand2020)

        Function to create array of travel time distributionss
        for distances from well.

        Parameter - Input
        ---------
        well_discharge_m3hour                           # [m3/hr]
        spreading_length               # [m]?, sqtr(K*D*c)
        vertical_resistance_aquitard   # [d], c_V
        porosity_shallow_aquifer                 # [-]
        KD                             # [m2/d], permeability * tickness of vertical layer
        permeability                   # [m2]

        thickness_vadose_zone      # [m], thickness of unsaturated zone/layer
        thickness_shallow_aquifer                 # [m], thickness of zone 1, if semiconfined aquitard,
                                        if phreatic aquifer, aquifer zone above top of well screens.
        thickness_target_aquifer                 # [m], thickness of zone 2

        thickness_full_capillary_fringe #[m], cF

        recharge_rate                        # [m/yr], recharge_rate of well area
        soil_moisture_content_vadose_zone           # [m3/m3], θ
        travel_time_H2O                 # [d],  travel time of water along flowline, in zone

        Output
        ------
        radial_distance_recharge        # [m], radial distance to well in order to recharge_rate the well
        radial_distance                 # [m], radial distance to well (field),
                                        from site X within and any site on the groundwater divide

        '''

        self.spreading_length = math.sqrt(self.vertical_resistance_aquitard * self.KD)

        self.radial_distance_recharge_semiconfined = self.spreading_length * 3

        self.radial_distance_recharge = (math.sqrt(self.well_discharge_m3hour
                                                    / (math.pi * self.recharge_rate / (365.25 * 24))))

        self.radial_distance, self.percent_flux = self._create_radial_distance_array()
        
        self.depth_bottom_shallow_aquifer = self.thickness_vadose_zone + self.thickness_shallow_aquifer
        self.depth_bottom_target_aquifer = self.thickness_vadose_zone + self.thickness_shallow_aquifer + self.thickness_target_aquifer


        # calculate the travel times in unsaturated zone, shallow aquifer and target aquifer in years
        self.travel_time_unsaturated = self._calculate_travel_time_unsaturated_zone()

        # travel time in semi-confined is one value, make it array by repeating the value
        self.travel_time_unsaturated = [self.travel_time_unsaturated] * (len(self.radial_distance))

        self.travel_time_shallow_aquifer = self._calculate_travel_time_aquitard_zone1()

        self.travel_time_target_aquifer = self._calculate_travel_time_target_aquifer() 

        self.total_travel_time = self.travel_time_unsaturated + \
            self.travel_time_shallow_aquifer + self.travel_time_target_aquifer

        self.head = self._calculuate_hydraulic_head()

        # percent of the max head? AH
        self.head_minus_max = 100*(self.head - self.groundwater_level) / (self.head[0] - self.groundwater_level)

        self.flux_fraction = self._calculate_flux_fraction()

        # [%], AH 1.1369 comes form pg. 52 in report, describes cutting off the 
        # recharge_rate distance at 3 labda, need to increase the percent abstracted from ~87%
        # to 99.9% so multiply by 1.1369 to get to that
        # AH, may want to change this, to eg. 6 labda or something else, adjust this number
        # equation A.16 in report
        self.cumulative_percent_abstracted_water = 1.1369 * 100 * (1 - self.flux_fraction)

        self.df_output = self._create_output_dataframe()

        return self.df_output

    def export_to_df(self, what_to_export='all'):
        """ Export to dataframe....

            Parameters
            ----------
            what_to_export: String
                options: 'all', 'omp_parameters', 'microbial_parameters'
            """
        
        #------------------------------
        # Make df_particle
        #------------------------------
        self.df_particle = pd.DataFrame(columns=['flowline_id', 'travel_time', 
                                                #  'xcoord','ycoord', # AH convert from radial distance? how?
                                                 "radial_distance",
                                                 'zcoord',
                                                 'redox_zone', 'temperature', 
                                                 'thickness_layer', 'porosity_layer', 
                                                 'dissolved_organic_carbon', 
                                                 'pH', 'fraction_organic_carbon'])

        df = self.df_particle.copy()

        for i in range(len(self.df_output)):
            flowline_id = i+1
            df.loc[0] = [flowline_id, 0, self.radial_distance[i],
                         0, None, self.temperature, None, None, None, 
                         None, None,]

            if self.schematisation == 'phreatic':
                vadose_zone = self.thickness_vadose_zone[i]
                bottom_vadose = self.depth_bottom_vadose_aquifer
            else: 
                vadose_zone = self.thickness_vadose_zone
                bottom_vadose = self.thickness_vadose_zone

            df.loc[1] = [flowline_id, self.travel_time_unsaturated[i],
                         self.radial_distance[i],
                         -1*vadose_zone,
                         self.redox_vadose_zone, self.temperature,
                         bottom_vadose,
                         self.porosity_vadose_zone,
                         self.dissolved_organic_carbon_vadose_zone,
                         self.pH_vadose_zone,
                         self.fraction_organic_carbon_vadose_zone
                         ]

            df.loc[2] = [flowline_id,
                         self.travel_time_unsaturated[i] +
                         self.travel_time_shallow_aquifer[i],
                         self.radial_distance[i],
                         -1*(self.depth_bottom_shallow_aquifer),
                         self.redox_shallow_aquifer, self.temperature,
                         self.thickness_shallow_aquifer,
                         self.porosity_shallow_aquifer,
                         self.dissolved_organic_carbon_shallow_aquifer,
                        self.pH_shallow_aquifer,
                         self.fraction_organic_carbon_shallow_aquifer
                         ]

            df.loc[3] = [flowline_id, self.total_travel_time[i],
                         self.well_screen_diameter/2,
                         -1*(self.depth_bottom_target_aquifer),
                         self.redox_target_aquifer, self.temperature,
                         self.thickness_target_aquifer,
                         self.porosity_target_aquifer,
                         self.dissolved_organic_carbon_target_aquifer,
                         self.pH_target_aquifer,
                         self.fraction_organic_carbon_target_aquifer
                         ]

            self.df_particle = self.df_particle.append(df, ignore_index=True)
            self.df_particle['density_aquifer'] = self.density_aquifer

        #------------------------------
        # Make df_flowline
        #------------------------------
        self.df_flowline = pd.DataFrame(columns=['flowline_id', 'discharge',
                                                 'particle_release_date',
                                                 'input_concentration',
                                                 'endpoint_id'])

        self.df_flowline['discharge'] = self.flowline_discharge
        self.df_flowline['flowline_id'] =  self.df_flowline.index + 1
        self.df_flowline['particle_release_date'] = self.particle_release_date
        self.df_flowline['input_concentration'] = self.input_concentration
        # #AH what is this? is this id zero?
        # self.df_flowline['endpoint_id'] = self.endpoint_id
        
        # AH which parameters for the 'microbial_parameters' option?

        if what_to_export == 'all' or what_to_export== 'omp_parameters':

            self.df_flowline['well_discharge_m3hour'] = self.well_discharge_m3hour
            self.df_flowline['recharge_rate'] = self.recharge_rate
            self.df_flowline['vertical_resistance_aquitard'] = self.vertical_resistance_aquitard
            self.df_flowline['KD'] = self.KD
            self.df_flowline['thickness_full_capillary_fringe'] = self.thickness_full_capillary_fringe
            self.df_flowline['recharge_concentration'] = self.recharge_concentration
            self.df_flowline['substance'] = self.substance
            self.df_flowline['input_concentration'] = self.input_concentration
            self.df_flowline['particle_release_date'] = self.particle_release_date
            self.df_flowline['soil_moisture_content_vadose_zone'] = self.soil_moisture_content_vadose_zone
            self.df_flowline['well_screen_diameter'] = self.well_screen_diameter
            self.df_flowline['removal_function'] = self.removal_function
            self.df_flowline['density_aquifer'] = self.density_aquifer

        if what_to_export == 'all':
            self.df_flowline['borehole_diameter'] = self.borehole_diameter
            self.df_flowline['k_hor_aquifer'] = self.k_hor_aquifer
            self.df_flowline['vani_aquifer'] = self.vani_aquifer
            self.df_flowline['k_hor_confining'] = self.k_hor_confining
            self.df_flowline['vani_confining'] = self.vani_confining
            self.df_flowline['k_hor_gravelpack'] = self.k_hor_gravelpack
            self.df_flowline['vani_gravelpack'] = self.vani_gravelpack
            self.df_flowline['k_hor_clayseal'] = self.k_hor_clayseal
            self.df_flowline['vani_clayseal'] = self.vani_clayseal
            self.df_flowline['dz_well'] = self.dz_well

        # #delete the unwanted columns depending on what the user asks for here
        return self.df_flowline, self.df_particle

    def plot_travel_time_versus_radial_distance(self,
                                                xlim=[0, 4000],
                                                ylim=[1, 5000]):
        ''' Plot the travel time versus the radial distance '''

        fig = plt.figure(figsize=[10, 5])
        plt.plot(self.radial_distance, self.total_travel_time, 'r', label=self.schematisation)
        # plt.plot(radial_distance, total_travel_time, 'b', label = 'Phreatic')
        plt.xlim(xlim)
        plt.ylim(ylim) 
        plt.yscale('log')
        plt.xlabel('Radial distance to well(s) (m)')
        plt.ylabel('Total travel time (years)')
        plt.title('Aquifer type: ' + self.schematisation)
        plt.grid()
        plt.savefig('travel_time_versus_radial_distance_'+self.schematisation+'.png', dpi=300, bbox_inches='tight')  # save_results_to + '/


    def plot_travel_time_versus_cumulative_abstracted_water(self,
                                                            xlim,
                                                            ylim=[1, 5000]):
        
        ''' Plot the travel time versus the cumulative abstracted water '''

        fig = plt.figure(figsize=[10, 5])
        plt.plot(self.cumulative_percent_abstracted_water, self.total_travel_time, 'r', label=self.schematisation)
        # plt.plot(radial_distance, total_travel_time, 'b', label = 'Phreatic')
        plt.xlim(xlim)  # [0.01,10000])
        plt.ylim(ylim) 
        plt.yscale('log')
        plt.xlabel('Cumulative percentage of abstracted water [%]')
        plt.ylabel('Total travel time (years)')
        plt.title('Aquifer type: ' + self.schematisation)
        plt.grid()
        plt.savefig('travel_time_versus_cumulative_percent_abstracted_water_'+self.schematisation+'.png', dpi=300, bbox_inches='tight')  # save_results_to + '/


class Concentration():
    """ Returns concentration in a groundwater well for a given Organic Micro Pollutant or microbial species.

    Parameters
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
        Column 'Kow'  # only necessary for OMP
        Column 'Labda'  # only necessary for microbiology 

    Returns
    -------    
    """
    def __init__(self, analytical_instance): #, substance, df_particle, df_flowline):
        # def __init__(self, substance: Substance, df_particle, df_flowline, removel_function?):
        self.omp_inialized = False
        self.df_particle = analytical_instance.df_particle
        self.df_flowline = analytical_instance.df_flowline


    def _init_omp(self):
        # AH this to be filled in from the substance class, which comes from existing AquaPriori tool?
        if self.omp_inialized:
            pass
        else:
            # self.df_part['Kow'] = self.df_part['redox_zone'].apply(lambda x: substance.get_Kow(x)
            #PLACEHOLDER to get some values for the Kow to test the rest of the function, FAKE VALUES
            self.df_particle['log_Koc'] = 1.92  # Koc for benzene PLACEHOLDER VALUE
            self.df_particle['pKa'] = 99        # pKa for benzene PLACEHOLDER VALUE

            # self.df_particle['log_Koc'] = -0.36 # Koc for AMPA PLACEHOLDER VALUE
            # self.df_particle['pKa'] = 0.4       # pKa for AMPA PLACEHOLDER VALUE

            # self.df_particle['log_Koc'] = 6.43  # Koc for benzo(a)pyrene PLACEHOLDER VALUE
            # self.df_particle['pKa'] = 99        # pKa for benzo(a)pyrene PLACEHOLDER VALUE

            # self.df_particle.loc[self.df_particle['red0x_zone'] == 0, 'Kow'] = 1.7
            # self.df_particle.loc[self.df_particle['red0x_zone'] == 1, 'Kow'] = 2.3
            # self.df_particle.loc[self.df_particle['red0x_zone'] == 2, 'Kow'] = 3.1
        self.omp_inialized = True


    def _init_microbiology():
        pass


    def compute_omp_removal(self):
        """ Returns the concentrain at each particle point.

        Paramneters
        -----------
        df_flowline, df_particle

        Returns
        -------
        df_particle: pandas.DataFrame
        extra column: 
        extra column: concentration
        extra column: Retardation
        extra column: break_through_time
        """
        self._init_omp()

        # Retardation equation based on Karickhoff (1981) and Schwarzenbach et al. (1993)
        # (section 10.3 in Appelo & Postma 2005), however with addition of 
        # the effects of (i) DOC-binding according to Kan & Tomson (1990),
        # and (ii) OMP ionization (dissociation) according to Schellenberg et al. (1984)
        # Equation 4.8-4.10 in in report

        # AH how are we incorporating the OMP database sheet into the model?
        # temperature correction for Koc
        '''Assuming the relation to be similar to the Van ‘t Hoff equation and equally performing for other OMPs yields:'''

        self.df_particle['Koc_temperature_correction'] = 10 **self.df_particle.log_Koc * 10 ** (1913 * (1 / (self.df_particle.temperature + 273.15) - 1 / (20 + 273.15)))

        self.df_particle['retardation'] = (1 + (1 / (1 + 10 ** (self.df_particle.pH - self.df_particle.pKa)) * self.df_particle.density_aquifer / 1000
                            * (1 - self.df_particle.porosity_layer)
                            * self.df_particle.fraction_organic_carbon * self.df_particle.Koc_temperature_correction)
                       / (self.df_particle.porosity_layer * (1 + (self.df_particle.Koc_temperature_correction * 1 / (1 + 10 ** (self.df_particle.pH - self.df_particle.pKa))
                                                            * 0.2 * self.df_particle.dissolved_organic_carbon * 0.000001))))

        # self.df_part...
            # LEFT OFF HERE, retardation works, need to add the rest to calculate concentration

    def compute_microbiology_removal(self):
        pass

    # def compute_well_concentration(self, evaluation_time = None)
    #     """ Returns the concentration in the raw water of each well."""
    #     if evaluation_time is None:
    #         select all flowlines
    #     else:
    #         select flowline with break_through_time < evaluation_time                                                  
    #         conc_flowline = concentration at end of selected flowline
    #         concentration_well = sum (conc_selected_flowline_i * discharge_flowline_i) / sum discharge_all_flowline                                                             


    def plot_concentration(self):
        pass

    def plot_age_distribution(self):
        pass

    def plot_logremoval(self):
        pass

# ------------------------------------------------------------------------------
# Basin Aquifer Recharge (BAR) Functions
# ------------------------------------------------------------------------------


def calculate_travel_time_recharge_basin(length,
                                width, 
                                depth,
                                annual_basin_discharge):
    '''
    Calculated the median travel time of infiltration water in a
    recharge basin, according to Huisman & Olsthoorn (1983)

    Parameters
    ----------
    width = width of recharge basin [m];
    length = total length of recharge basin and drainage gallery [m];
    depth = mean water-filled depth of recharge basin [m];
    annual_basin_discharge = total volume of infiltration water discharged
        into basin [m3/yr]

    x = longitudinal distance from inlet, situated on one of the ends of an
        elongate basin [m];
    calculated for median distance for x = 0.632 * length #AH why 0.632?

    Returns
    -------
    travel_time: [years]
    '''

    travel_time_recharge_basin = (length * width * depth
                                  / (annual_basin_discharge * 10 ** 6
                                     / 365.24) * np.log(1 / (1 - 0.632)))
    # AH why 0.632 = median?

    return travel_time_recharge_basin


def calculate_travel_time_upper_aquifer( 
                    horizontal_distance_basin_gallery,   
                    hydraulic_conductivity_reference_value,
                    temperature,
                    porosity_recharge_basin,
                    groundwater_level_above_saturated_zone,
                     ):
    '''
    The travel time in the upper aquifer is based on design
    characteristics of the BAR system, which is rigorously
    schematized as indicated in Fig.A3.1. The schematisation
    and analytical solutions are based on Huisman & Olsthoorn
    (1983) and Stuyfzand & van der Schans (2018), with a
    MAR system composed of parallel spreading basins and
    drainage galleries in between, with phreatic aquifer
    characteristics and typical recharge system requirements.

    The design follows from the equations given below, assuming that:
    (i) flow is predominantly horizontal (Dupuit assumption;
    also excluding unsaturated zone beneath basin), requiring
    that both the spreading basin and drainage gallery fully
    penetrates the saturated thickness of the aquifer;
    (ii) the aquifer is homogeneous and isotropic,
    (iii) rainfall and evaporation can be ignored;
    (iv) all infiltration water is recovered (QIN = QOUT).

    Parameters
    ----------
    horizontal_distance_basin_gallery = horizontal distance between basin bank and drainage gallery [m];
    hydraulic_conductivity_reference_value = m/d, at specific deg. C
    temperature = deg C
    porosity_recharge_basin = porosity [-]
    groundwater_level_above_saturated_zone = normal maximum rise of watertable above H0 [m];
    
    Outputs
    -------
    median_retention_time = minimum required or median retention time in aquifer [d];
    travel_time_distribution_upper_phreatic = travel time in the upper aquifer [d]

    Other paramters
    --------------
    hydraulic_conductivity = hydraulic conductivity of aquifer [m/d];


    Not used but calculated/input in spreadsheet AH check what to do with these?
    ---------------------------
    tIP = maximum bridging period from the dynamic reservoir,
        which allows abstraction during an intake (and spreading) interruption [d];
    VDYN,A, VDYN,B = amount of water in dynamic storage, within the aquifer
        and basin respectively, which can be used to cover a (short) intake stop [Mm3];
    vIN = entry rate [m/d]
    QOUT = pumping rate, approximately equal to QIN prior to intake stop,
        however on a daily basis [m3/d];
    SY = phreatic specific yield [-];
    H0 = saturated thickness of aquifer before spreading [m];

    '''

    # horizontal_distance_basin_gallery = np.sqrt(
    #     hydraulic_conductivity * groundwater_level_above_saturated_zone * median_retention_time / porosity_recharge_basin)
    hydraulic_conductivity = (hydraulic_conductivity_reference_value
                              * ((temperature + 42.5) / (11 + 42.5)) ** 1.5)
    # AH where do these numbers come from? 11, 42.5, 1.5?

    median_retention_time = (porosity_recharge_basin * horizontal_distance_basin_gallery
                             ** 2 / (hydraulic_conductivity * groundwater_level_above_saturated_zone))

    travel_time_array = np.arange(0.50, 4.6, 0.1)

    travel_time_distribution_upper_phreatic = median_retention_time * travel_time_array

    # EPM = exponential piston model, AH percent travel time?
    percent_flux_EPM = 1 - np.exp((travel_time_distribution_upper_phreatic[0] - travel_time_distribution_upper_phreatic) / (
        median_retention_time - travel_time_distribution_upper_phreatic[0]))

    # LMP = linear piston model, AH percent travel time?
    percent_flux_LPM = np.arange(0, 1.1, 0.1)

    # AH same as the percent flux? keep variables consistent...
    # percent_travel_time_distribution = np.arange()
    percent_flux = np.append(percent_flux_LPM[0:9], percent_flux_EPM[9:])

    aquifer_dictionary = {'travel_time_basin': median_retention_time,
                          'travel_time_upper_aquifer': travel_time_distribution_upper_phreatic,
                          'percent_flux': percent_flux,
                          'percent_flux_EPM':  percent_flux_EPM,
                          'percent_flux_LPM': percent_flux_LPM,
                          'aquifer_type': 'BAR_upper_phreatic',
                          # 'travel_time_array': travel_time_array,
                          }

    return aquifer_dictionary


def calculate_travel_time_upper_aquifer_and_basins(length,
                                                   width,
                                                   depth,
                                                   annual_basin_discharge,
                                                   horizontal_distance_basin_gallery,
                                                   hydraulic_conductivity_reference_value,
                                                   temperature,
                                                   porosity_recharge_basin,
                                                   groundwater_level_above_saturated_zone,
                                                   ):

    '''Adds the median travel time in the recharge basins to the travel 
    time calculated in the upper, phreatic aquifer
    
    Inputs
    ------
    width = width of recharge basin [m];
    length = total length of recharge basin and drainage gallery [m];
    depth = mean water-filled depth of recharge basin [m];


    Returns
    -------
    Travel time distribution in the upper aquifer and basin '''

    travel_time_recharge_basin = calculate_travel_time_recharge_basin(
        length=length,
        width=width,
        depth=depth,
        annual_basin_discharge=annual_basin_discharge)

    upper_phreatic_aquifer_dict = calculate_travel_time_upper_aquifer(
        horizontal_distance_basin_gallery=horizontal_distance_basin_gallery,
        hydraulic_conductivity_reference_value=hydraulic_conductivity_reference_value,
        temperature=temperature,
        porosity_recharge_basin=porosity_recharge_basin,
        groundwater_level_above_saturated_zone=groundwater_level_above_saturated_zone,
        )

    travel_time_upper_aquifer_and_basins = travel_time_recharge_basin + upper_phreatic_aquifer_dict['travel_time_upper_aquifer']

    return travel_time_upper_aquifer_and_basins

    #AH LEFT OFF HERE, NEED TO CONFIRM WITH MARTIN vdS WHAT WE WANT OUTPUT FROM THE BAR SHEET


def calculate_travel_time_deeper_aquifer(horizontal_distance_basin_gallery,
                                         hydraulic_conductivity_reference_value,
                                         temperature,
                                         porosity_recharge_basin,
                                         groundwater_level_above_saturated_zone,
                                         vertical_resistance_aquitard,
                                         thickness_second_aquifer,
                                         porosity_second_aquifer,
                                         vertical_hydraulic_head_difference,
                                         ):
    '''
    For the deeper part we have 2 options: (i) it represents
    the second aquifer, which is semiconfined and (deeply)
    anoxic situated below an aquitard that separates it from
    the upper aquifer, or (ii) it represents the deeper, anoxic
    parts of the upper aquifer, often with lower hydraulic
    conductivity and longer transit times. In the latter
    case, we neglect the contribution from the second aquifer
    if existing at all
    
    Paramters
    ---------
    vertical_resistance_aquitard   # [d], c_V
    thickness_second_aquifer: mean thickness [m]
    porosity_second_aquifer (n), [-]
    vertical_hydraulic_head_difference: mean vertical jump 
        in hydraulic head between both aquifers (Δh).
    
    Outputs
    -------
    '''

    travel_time_second_aquifer = (0.3 *  thickness_second_aquifer
        * vertical_resistance_aquitard / vertical_hydraulic_head_difference)
    
    percent_flux_deep_aquifer = 2000 / travel_time_second_aquifer

    multiplication_factor = np.sqrt(travel_time_second_aquifer)

    percent_flux = np.arange(0,110,10)

    # EPM = exponential piston model
    travel_time_distribution_EPM = percent_flux / 100

    travel_time_upper_aquifer_dict = calculate_travel_time_upper_aquifer(
        horizontal_distance_basin_gallery=horizontal_distance_basin_gallery,
        hydraulic_conductivity_reference_value=hydraulic_conductivity_reference_value,
        temperature=temperature,
        porosity_recharge_basin=porosity_recharge_basin,
        groundwater_level_above_saturated_zone=groundwater_level_above_saturated_zone,
        )

    travel_time_distribution_deeper_aquifer = (travel_time_upper_aquifer_dict['travel_time_upper_aquifer']
                                                * multiplication_factor)

    return travel_time_distribution_deeper_aquifer


def travel_time_distribution_BAR(length,
                                 width,
                                 depth,
                                 annual_basin_discharge,
                                 horizontal_distance_basin_gallery,
                                 hydraulic_conductivity_reference_value,
                                 temperature,
                                 porosity_recharge_basin,
                                 groundwater_level_above_saturated_zone,
                                 vertical_resistance_aquitard,
                                 thickness_second_aquifer,
                                 porosity_second_aquifer,
                                 vertical_hydraulic_head_difference,
                                 ):
    
    ''' Calculate travel time distribution for Basin Aquifer Recharge (BAR)
    systems

    Calculates the median travel time distribution (TTD) of the basin, 
    the TTD of the upper aquifer and TTD of the deeper aquifer.

    Assumptions: from \cite{Stuyfzand2020}
    The following assumptions are made: 
    (i) the water inlet is situated at one of the ends of an elongated 
    recharge basin 
    (ii) the infiltration rate is uniformly distributed
    (iii) rainfall and evaporation do not have a significant impact.

    Parameters
    -----------
    width = width of recharge basin [m];
    length = total length of recharge basin and drainage gallery [m];
    depth = mean water-filled depth of recharge basin [m];

    horizontal_distance_basin_gallery = horizontal distance between basin bank and drainage gallery [m];
    hydraulic_conductivity_reference_value = m/d, at specific deg. C
    temperature = deg C
    porosity_recharge_basin = porosity [-]
    groundwater_level_above_saturated_zone = normal maximum rise of watertable above H0 [m];

    '''

    median_retention_time = calculate_travel_time_recharge_basin(
        length=length,
        width=width,
        depth=depth,
        annual_basin_discharge=annual_basin_discharge,
        )

    upper_phreatic_aquifer_dict = calculate_travel_time_upper_aquifer(
        horizontal_distance_basin_gallery=horizontal_distance_basin_gallery,
        hydraulic_conductivity_reference_value=hydraulic_conductivity_reference_value,
        temperature=temperature,
        porosity_recharge_basin=porosity_recharge_basin,
        groundwater_level_above_saturated_zone=groundwater_level_above_saturated_zone,
        )

    travel_time_upper_aquifer_and_basins = calculate_travel_time_upper_aquifer_and_basins(
        length=length,
        width=width,
        depth=depth,
        annual_basin_discharge=annual_basin_discharge,
        horizontal_distance_basin_gallery=horizontal_distance_basin_gallery,
        hydraulic_conductivity_reference_value=hydraulic_conductivity_reference_value,
        temperature=temperature,
        porosity_recharge_basin=porosity_recharge_basin,
        groundwater_level_above_saturated_zone=groundwater_level_above_saturated_zone,
        )

    travel_time_deeper_aquifer = calculate_travel_time_deeper_aquifer(
        horizontal_distance_basin_gallery=horizontal_distance_basin_gallery,
        hydraulic_conductivity_reference_value=hydraulic_conductivity_reference_value,
        temperature=temperature,
        porosity_recharge_basin=porosity_recharge_basin,
        groundwater_level_above_saturated_zone=groundwater_level_above_saturated_zone,
        vertical_resistance_aquitard=vertical_resistance_aquitard,
        thickness_second_aquifer=thickness_second_aquifer,
        porosity_second_aquifer=porosity_second_aquifer,
        vertical_hydraulic_head_difference=vertical_hydraulic_head_difference,
        )

    column_names = ["travel_time_upper_aquifer", "travel_time_upper_aquifer_and_basins",
                        "travel_time_deeper_aquifer", "percent_flux", "well_discharge_m3hour"]
    
    # AH check the flux calculations for the streamline... Pr*Q? 
    data = [upper_phreatic_aquifer_dict['travel_time_upper_aquifer'],
            travel_time_upper_aquifer_and_basins,
            travel_time_deeper_aquifer,
            upper_phreatic_aquifer_dict['percent_flux'],
            upper_phreatic_aquifer_dict['percent_flux'] * annual_basin_discharge]

    df = pd.DataFrame (data = np.transpose(data), columns=column_names)

    aquifer_dictionary = {'travel_time_basin': median_retention_time,
                          'travel_time_upper_aquifer': upper_phreatic_aquifer_dict['travel_time_upper_aquifer'],
                          'travel_time_upper_aquifer_and_basins': travel_time_upper_aquifer_and_basins,
                          'travel_time_deeper_aquifer': travel_time_deeper_aquifer,
                          'percent_flux': upper_phreatic_aquifer_dict['percent_flux'],
                          'df_output': df,
                          'aquifer_type': 'BAR',
                         }

    return aquifer_dictionary




#%%
# LEFT OFF HERE
# ask Martin K how to format/structure all these functions/classes

#%%

#   QUESITONS
# 2. input permeability, thickness seperately?
# 4. naming conventions: zone 1, zone 2, zone 3 or unsaturated zone, zone 1, zone 2 .. 
#   gets more complicated in the BAR (basin, phreatic/first/upper aquifer, deep/aquitard/second aquifer)
# 5. error in spreadsheet in "Peter's t2 fit", does not seem to take into account changes in the thickness of the aquifer and porosity (fixed n2 = 0.32 and D2 = 95)...?
# 6. Calculating the well_discharge_m3hour (q) for ech streamline, Phreatic %Qrx*Q, Semi-confinded: Qr/Qx*Qx?? see notes and df_output for each
# %%
