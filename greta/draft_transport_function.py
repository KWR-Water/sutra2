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
# specific questions flagged for;
# @MartinvdS // @steven //@martinK
####

# LEFT OFF : APRIL 29 
# -> finished adding the breakthrough concentration calculations (time and concentrations at each layer) 
# for the phreatic case. Confim with Martin vd S about the inidividual questions marked in code and next steps

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
schematisation_type: string
    'phreatic', 'semi-confined', 'riverbankfiltration', 'basinfiltration'
removal_function: string
    'omp' -> get parameters for OMP
    'microbiology' -> get parameters for microbiology
"""


#%%
class HydroChemicalSchematisation:

    """ Converts input parameters of AquaPriori GUI to a complete parameterisation."""

    #everything initialized to None to start...
    def __init__(self,
                 schematisation_type=None,
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
                 KD=None,  # AH premeability*thicnkeness or input seperately? @MartinvdS
                 groundwater_level_ASL = None,
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
        # AH @MartinvdS can/should we (user) be able to alter these?
        # if so include these as default arguements for the class?
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
        self.schematisation_type = schematisation_type 
        self.groundwater_level_ASL = groundwater_level_ASL
        self.thickness_vadose_zone = thickness_vadose_zone
        self.thickness_shallow_aquifer = thickness_shallow_aquifer
        self.thickness_target_aquifer = thickness_target_aquifer

        self.porosity_vadose_zone = porosity_vadose_zone
        self.porosity_shallow_aquifer = porosity_shallow_aquifer
        self.porosity_target_aquifer = porosity_target_aquifer

        self.dissolved_organic_carbon_vadose_zone = dissolved_organic_carbon_vadose_zone
        self.dissolved_organic_carbon_shallow_aquifer = dissolved_organic_carbon_shallow_aquifer
        self.dissolved_organic_carbon_target_aquifer = dissolved_organic_carbon_target_aquifer

        # 'suboxic'; 'anoxic'; 'deeply_anoxic'
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

        # AH @MartinvdS input these seperately? # [m2/d], permeability * thickness of vertical layer
        self.KD = KD
        self.groundwater_level = self.groundwater_level_ASL # 0-self.thickness_vadose_zone
        self.thickness_full_capillary_fringe = thickness_full_capillary_fringe

        self.recharge_concentration = recharge_concentration

        self.substance = substance

        # AH @MartinvdS are the input concentration and release dates single values? or do they vary over time/distance?
        self.input_concentration=input_concentration
        self.particle_release_date=particle_release_date

    def make_dictionary(self,):
        ''' Dicitonaries of the different paramters '''
    # make dictionaries of each grouping of parameters
    
        recharge = {
            'recharge1': {
                'recharge_rate': self.recharge_rate
                # 'compound': ....
                # 'concentration': ....
                #'x_start':,
                #x_end':,
            }}

        geo_parameters = {
            'geolayer1': {
                'vadose': True,  # not part of modflow
                'top': 0,  # flow parameters
                'bot': 0 - self.thickness_vadose_zone,
                # 'K_hor': self.k_hor_aquifer, # AH what goes here? @MartinvdS
                # 'K_vert': self.K_vert, # AH what goes here? @MartinvdS
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
                'K_hor': self.k_hor_aquifer, #these paramters need to be updated with the name of hte layers #AH
                # 'K_vert': self.K_vert, # AH what goes here? @MartinvdS
                'porosity': self.porosity_shallow_aquifer,
                'redoxzone': self.redox_shallow_aquifer,  # concentration parameters
                'dissolved_organic_carbon': self.dissolved_organic_carbon_shallow_aquifer,
                'vertical_resistance_aquitard': self.vertical_resistance_aquitard, # AH with the 2nd layer or 3rd layer? @MartinvdS
            },
            'geolayer3': {
                'vadose': False,  # not part of modflow
                'top': 0 - self.thickness_vadose_zone - self.thickness_shallow_aquifer,  # flow parameters
                'bot': (0 - self.thickness_vadose_zone
                        - self.thickness_shallow_aquifer
                        - self.thickness_target_aquifer),
                'K_hor': self.k_hor_aquifer,
                # 'K_vert': self.K_vert, # AH what goes here? @MartinvdS
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
                'discharge': - self.well_discharge_m3hour / 24, # AH  @MartinvdS should this be '* 24' to make  m3/day?
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

        # AH @MartinvdS what to put in these dicts (not used in the analytical functions)
        # backfill_parameters = {
        #     'backfill1': {
        #         'boreholeradius': self.borehole_diameter / 2,
        #         # 'K_hor': self.k_hor_clayseal # AH which clayseal? or gravelpack? @MartinvdS
        #     }
        # }

        # drain_parameters = {
        #     'depth':
        #     'discharge':
        #     }

        # AH what goes here? @MartinvdS
        # surfacewater_parameters = {
        # 'water_level': ....
        #     }

        # AH how do we want the dictionary made? @MartinvdS
        # it is not very convenient for the analytical well models to have in one dictionary, 
        # but outputting many dictionaries is also inconvenient
        hydrochemical_dictionary = {'well_parameters': well_parameters,
                                    'geo_parameters': geo_parameters,
                                    'recharge': recharge, }
        return hydrochemical_dictionary


#%%
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

    def _check_required_variables(self,required_variables):
        for req_var in required_variables:
            value = getattr(self.schematisation, req_var)
            if value is None:
                raise KeyError(f'Error, required variable {req_var} is not defined.')

    def _check_init_phreatic(self):
        '''check the variables that we need for the individual aquifer types are not NONE aka set by the user'''

        required_variables = ["schematisation_type", #repeat for all
                              "thickness_vadose_zone",
                              "thickness_shallow_aquifer",
                              "thickness_target_aquifer",
                              "porosity_vadose_zone",
                              "porosity_shallow_aquifer",
                              "porosity_target_aquifer",
                              "well_discharge_m3hour",
                              "recharge_rate",
                              "vertical_resistance_aquitard",
                              "soil_moisture_content_vadose_zone",
                              "KD",
                              "groundwater_level",
                              "thickness_full_capillary_fringe",
                            #   self.test_variable,
                            ]
        self._check_required_variables(required_variables)

    def _check_init_confined(self):
        '''check the variables that we need for the individual aquifer types are not NONE aka set by the user'''
       
        required_variables = ["schematisation_type", #repeat for all
                              "thickness_vadose_zone",
                              "thickness_shallow_aquifer",
                              "thickness_target_aquifer",
                              "porosity_vadose_zone",
                              "porosity_shallow_aquifer",
                              "porosity_target_aquifer",
                              "well_discharge_m3hour",
                              "recharge_rate",
                              "vertical_resistance_aquitard",
                              "soil_moisture_content_vadose_zone",
                              "KD",
                              "groundwater_level",
                              "thickness_full_capillary_fringe",
                            #   self.test_variable,
                            ]
        self._check_required_variables(required_variables)

    def _create_radial_distance_array(self):

        ''' Calculate the radial distance from a well,
        specifying the aquifer type (semi-confined, phreatic)
        '''
        percent_flux = np.array([0.001, 0.01, 0.1, 0.5])
        percent_flux = np.append(percent_flux, np.arange(1, 100, 1))
        self.percent_flux = np.append(percent_flux, [99.5, 99.99])

        radial_distance = self.radial_distance_recharge * \
            np.sqrt(self.percent_flux / 100)

        if self.schematisation.schematisation_type == 'semi-confined':

            radial_distance = np.append(radial_distance,
                                    [(radial_distance[-1] + ((self.spreading_length * 3) - radial_distance[-1]) / 3),
                                    (radial_distance[-1] + 2 *
                                    ((self.spreading_length * 3) - radial_distance[-1]) / 3),
                                        (self.spreading_length * 3)])
        
        self.radial_distance = radial_distance
        return self.radial_distance, self.percent_flux

    def _calculate_travel_time_unsaturated_zone(self):
        '''Equation A.11 in report '''

        # AH
        # (groundwater_level +thickness_vadose_zone) -
        # > before was the top of the unsaturated zone, 
        # confirm if this is the same?
        # aka is the top of the unsaturated zone always
        # the thickness_vadose_zone + groundwater_level
        self.travel_time_unsaturated = ((((self.schematisation.groundwater_level + self.thickness_vadose_zone_drawdown)
                                    - self.schematisation.groundwater_level
                                    - self.schematisation.thickness_full_capillary_fringe)
                                    * self.schematisation.soil_moisture_content_vadose_zone
                                    + self.schematisation.thickness_full_capillary_fringe
                                    * self.schematisation.porosity_vadose_zone)
                                / self.schematisation.recharge_rate)

        return self.travel_time_unsaturated


    def _calculate_travel_time_aquitard_zone1(self):
        ''' Calculate the travel time in zone 1 (aquitard in semi-confined case)
        using the the Peters (1985) solution( eq. 8.8 in Peters \cite{Peters1985})
        Equation A.12 in report
        '''                                   

        self.travel_time_shallow_aquifer = (2 * math.pi * self.schematisation.KD * self.schematisation.vertical_resistance_aquitard
                            / (self.schematisation.well_discharge_m3hour * 24)
                            * (self.schematisation.thickness_shallow_aquifer
                                / besselk(0, self.radial_distance
                                        / math.sqrt(self.schematisation.KD * self.schematisation.vertical_resistance_aquitard)))
                            / 365.25)

        return self.travel_time_shallow_aquifer


    def _calculate_travel_time_target_aquifer(self):
        '''Calculate the travel time in zone 2 (aquifer with the production well)
        using the the Peters (1985) solution
        Equation A.13/A.14 in report

        '''
        # porosity_target_aquifer, AH #we keep this number for camparing to excel, but will change later to be the user defined porosities
        porosity_target_aquifer=0.32  
        thickness_target_aquifer=95

        self.travel_time_target_aquifer = (2 * math.pi * self.spreading_length ** 2 / (self.schematisation.well_discharge_m3hour * 24)
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

        self.head = (-self.schematisation.well_discharge_m3hour * 24 / (2 * math.pi * self.schematisation.KD)
                * besselk(0, self.radial_distance / self.spreading_length)
                + self.schematisation.groundwater_level)

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
        # AH, confirm this is how to calculate the discharge of each flowline?
        # Percent flux ( array of 1,,23,....,99,99.5,99.9), row 7 in TTD phreatic, difference between cells divided by the well discharge
        # @MartinvdS
        self.flowline_discharge = (np.diff( np.insert(self.cumulative_percent_abstracted_water,0,0., axis=0))/100)*self.schematisation.well_discharge_m3hour

        column_names = ["total_travel_time", "travel_time_unsaturated", 
                        "travel_time_shallow_aquifer", "travel_time_target_aquifer",
                        "radial_distance", "head", "cumulative_percent_abstracted_water", 
                        "flowline_discharge"
                        ]

        # AH check the well_discharge_m3hour calculations for the streamline... Pr*Q? 
        data = [self.total_travel_time, self.travel_time_unsaturated, 
                    self.travel_time_shallow_aquifer, self.travel_time_target_aquifer,
                    self.radial_distance, self.head, self.cumulative_percent_abstracted_water, 
                    self.flowline_discharge
                    ]

        self.df_output = pd.DataFrame (data = np.transpose(data), columns=column_names)
        
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

        self.spreading_length = math.sqrt(self.schematisation.vertical_resistance_aquitard * self.schematisation.KD)

        self.radial_distance_recharge = (math.sqrt(self.schematisation.well_discharge_m3hour
                                                    / (math.pi * self.schematisation.recharge_rate / (365.25 * 24))))

        self.radial_distance, self.percent_flux = self._create_radial_distance_array()

        self.head = (self.schematisation.groundwater_level - self.schematisation.well_discharge_m3hour * 24 
                    / (2 * math.pi * self.schematisation.KD) 
                    * np.log(self.radial_distance_recharge / self.radial_distance))

        self.depth_bottom_vadose_aquifer = self.schematisation.thickness_vadose_zone #no drawdown
        self.depth_bottom_shallow_aquifer = self.schematisation.thickness_vadose_zone + self.schematisation.thickness_shallow_aquifer
        self.depth_bottom_target_aquifer = self.schematisation.thickness_vadose_zone + self.schematisation.thickness_shallow_aquifer + self.schematisation.thickness_target_aquifer

        self.thickness_vadose_zone_drawdown = (self.schematisation.groundwater_level 
                                              + self.schematisation.thickness_vadose_zone) - self.head

        self.travel_time_unsaturated = self._calculate_travel_time_unsaturated_zone()

        self.travel_time_shallow_aquifer = ((self.schematisation.thickness_shallow_aquifer - (self.schematisation.groundwater_level - self.head))
                            * self.schematisation.porosity_shallow_aquifer / self.schematisation.recharge_rate)

        self.travel_time_target_aquifer = (self.schematisation.porosity_target_aquifer * self.schematisation.thickness_target_aquifer / self.schematisation.recharge_rate
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

        self.spreading_length = math.sqrt(self.schematisation.vertical_resistance_aquitard * self.schematisation.KD)

        self.radial_distance_recharge_semiconfined = self.spreading_length * 3

        self.radial_distance_recharge = (math.sqrt(self.schematisation.well_discharge_m3hour
                                                    / (math.pi * self.schematisation.recharge_rate / (365.25 * 24))))

        self.radial_distance, self.percent_flux = self._create_radial_distance_array()
        
        self.depth_bottom_shallow_aquifer = self.schematisation.thickness_vadose_zone + self.schematisation.thickness_shallow_aquifer
        self.depth_bottom_target_aquifer = self.schematisation.thickness_vadose_zone + self.schematisation.thickness_shallow_aquifer + self.schematisation.thickness_target_aquifer


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
        self.head_minus_max = 100*(self.head - self.schematisation.groundwater_level) / (self.head[0] - self.schematisation.groundwater_level)

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
        df_particle = pd.DataFrame(columns=['flowline_id', 'zone', 'travel_time_zone', 'total_travel_time', 
                                                #  'xcoord','ycoord', # AH convert from radial distance? @MartinvdS
                                                 "radial_distance",
                                                 'zcoord',
                                                 'redox_zone', 'temperature', 
                                                 'thickness_layer', 'porosity_layer', 
                                                 'dissolved_organic_carbon', 
                                                 'pH', 'fraction_organic_carbon', 
                                                 'input_concentration', 'steady_state_concentration'])

        df = df_particle.copy()

        for i in range(len(self.df_output)):
            flowline_id = i+1
            input_concentration = self.schematisation.input_concentration #df_flowline.loc[self.df_flowline['flowline_id'] == 1, 'input_concentration'].values[0]

            df.loc[0] = [flowline_id, None, 0, 0, self.radial_distance[i],
                         0, None, self.schematisation.temperature, None, None, None, 
                         None, None,input_concentration, input_concentration]

            if self.schematisation == 'phreatic':
                vadose_zone = self.schematisation.thickness_vadose_zone[i]
                bottom_vadose = self.depth_bottom_vadose_aquifer
            else: 
                vadose_zone = self.schematisation.thickness_vadose_zone
                bottom_vadose = self.schematisation.thickness_vadose_zone

            df.loc[1] = [flowline_id, 
                         "vadose_zone",
                         self.travel_time_unsaturated[i],
                         self.travel_time_unsaturated[i],
                         self.radial_distance[i],
                         -1*vadose_zone,
                         self.schematisation.redox_vadose_zone, 
                         self.schematisation.temperature,
                         bottom_vadose,
                         self.schematisation.porosity_vadose_zone,
                         self.schematisation.dissolved_organic_carbon_vadose_zone,
                         self.schematisation.pH_vadose_zone,
                         self.schematisation.fraction_organic_carbon_vadose_zone,
                         input_concentration,
                         None,
                         ]

            df.loc[2] = [flowline_id, "shallow_aquifer", self.travel_time_shallow_aquifer[i],
                         self.travel_time_unsaturated[i] + self.travel_time_shallow_aquifer[i],
                         self.radial_distance[i],
                         -1*(self.depth_bottom_shallow_aquifer),
                         self.schematisation.redox_shallow_aquifer, 
                         self.schematisation.temperature,
                         self.schematisation.thickness_shallow_aquifer,
                         self.schematisation.porosity_shallow_aquifer,
                         self.schematisation.dissolved_organic_carbon_shallow_aquifer,
                         self.schematisation.pH_shallow_aquifer,
                         self.schematisation.fraction_organic_carbon_shallow_aquifer,
                         input_concentration, 
                         None,
                         ]

            df.loc[3] = [flowline_id, "target_aquifer",  self.travel_time_target_aquifer[i],
                         self.total_travel_time[i],
                         self.schematisation.well_screen_diameter/2,
                         -1*(self.depth_bottom_target_aquifer),
                         self.schematisation.redox_target_aquifer, 
                         self.schematisation.temperature,
                         self.schematisation.thickness_target_aquifer,
                         self.schematisation.porosity_target_aquifer,
                         self.schematisation.dissolved_organic_carbon_target_aquifer,
                         self.schematisation.pH_target_aquifer,
                         self.schematisation.fraction_organic_carbon_target_aquifer,
                         input_concentration,
                         None
                         
                         ]

            df_particle = df_particle.append(df, ignore_index=True)
            df_particle['density_aquifer'] = self.schematisation.density_aquifer
            df_particle['redox_zone'] = df_particle['redox_zone'].fillna('').astype(str)
            df_particle['flowline_id'] = df_particle['flowline_id'].astype(int)
            # df_particle = df_particle.replace({np.nan: None})
            self.df_particle = df_particle

        #------------------------------
        # Make df_flowline
        #------------------------------
        df_flowline = pd.DataFrame(columns=['flowline_id', 'discharge',
                                                 'particle_release_date',
                                                 'input_concentration',
                                                 'endpoint_id'])

        df_flowline['discharge'] = self.flowline_discharge
        df_flowline['flowline_id'] =  df_flowline.index + 1
        df_flowline['particle_release_date'] = self.schematisation.particle_release_date
        df_flowline['input_concentration'] = self.schematisation.input_concentration
        # #AH what is this? is this id zero? @MartinvdS
        # df_flowline['endpoint_id'] = self.endpoint_id
        
        # AH which parameters for the 'microbial_parameters' option? @MartinvdS or @steven

        if what_to_export == 'all' or what_to_export== 'omp_parameters':

            df_flowline['well_discharge_m3hour'] = self.schematisation.well_discharge_m3hour
            df_flowline['recharge_rate'] = self.schematisation.recharge_rate
            df_flowline['vertical_resistance_aquitard'] = self.schematisation.vertical_resistance_aquitard
            df_flowline['KD'] = self.schematisation.KD
            df_flowline['thickness_full_capillary_fringe'] = self.schematisation.thickness_full_capillary_fringe
            df_flowline['recharge_concentration'] = self.schematisation.recharge_concentration
            df_flowline['substance'] = self.schematisation.substance
            df_flowline['input_concentration'] = self.schematisation.input_concentration
            df_flowline['particle_release_date'] = self.schematisation.particle_release_date
            df_flowline['soil_moisture_content_vadose_zone'] = self.schematisation.soil_moisture_content_vadose_zone
            df_flowline['well_screen_diameter'] = self.schematisation.well_screen_diameter
            df_flowline['removal_function'] = self.schematisation.removal_function
            df_flowline['density_aquifer'] = self.schematisation.density_aquifer

        if what_to_export == 'all':
            df_flowline['borehole_diameter'] = self.schematisation.borehole_diameter
            df_flowline['k_hor_aquifer'] = self.schematisation.k_hor_aquifer
            df_flowline['vani_aquifer'] = self.schematisation.vani_aquifer
            df_flowline['k_hor_confining'] = self.schematisation.k_hor_confining
            df_flowline['vani_confining'] = self.schematisation.vani_confining
            df_flowline['k_hor_gravelpack'] = self.schematisation.k_hor_gravelpack
            df_flowline['vani_gravelpack'] = self.schematisation.vani_gravelpack
            df_flowline['k_hor_clayseal'] = self.schematisation.k_hor_clayseal
            df_flowline['vani_clayseal'] = self.schematisation.vani_clayseal
            df_flowline['dz_well'] = self.schematisation.dz_well

        self.df_flowline = df_flowline
        # #delete the unwanted columns depending on what the user asks for here
        return self.df_flowline, self.df_particle

    def plot_travel_time_versus_radial_distance(self,
                                                xlim=[0, 4000],
                                                ylim=[1, 5000]):
        ''' Plot the travel time versus the radial distance '''

        fig = plt.figure(figsize=[10, 5])
        plt.plot(self.radial_distance, self.total_travel_time, 'r', label=self.schematisation.schematisation_type)
        # plt.plot(radial_distance, total_travel_time, 'b', label = 'Phreatic')
        plt.xlim(xlim)
        plt.ylim(ylim) 
        plt.yscale('log')
        plt.xlabel('Radial distance to well(s) (m)')
        plt.ylabel('Total travel time (years)')
        plt.title('Aquifer type: ' + self.schematisation.schematisation_type)
        plt.grid()
        plt.savefig('travel_time_versus_radial_distance_'+self.schematisation.schematisation_type+'.png', dpi=300, bbox_inches='tight')  # save_results_to + '/


    def plot_travel_time_versus_cumulative_abstracted_water(self,
                                                            xlim,
                                                            ylim=[1, 5000]):
        
        ''' Plot the travel time versus the cumulative abstracted water '''

        fig = plt.figure(figsize=[10, 5])
        plt.plot(self.cumulative_percent_abstracted_water, self.total_travel_time, 'r', label=self.schematisation.schematisation_type)
        # plt.plot(radial_distance, total_travel_time, 'b', label = 'Phreatic')
        plt.xlim(xlim)  # [0.01,10000])
        plt.ylim(ylim) 
        plt.yscale('log')
        plt.xlabel('Cumulative percentage of abstracted water [%]')
        plt.ylabel('Total travel time (years)')
        plt.title('Aquifer type: ' + self.schematisation.schematisation_type)
        plt.grid()
        plt.savefig('travel_time_versus_cumulative_percent_abstracted_water_'+self.schematisation.schematisation_type+'.png', dpi=300, bbox_inches='tight')  # save_results_to + '/

class Substance:
    def __init__(self, name, ):
        """
        name: String, 
            name of the substance (for now limited dictionary to 'benzene', 'AMPA', 'benzo(a)pyrene'
        substance_dict: dictionary
            log Koc: distribution coefficient of organic carbon and water ([-]
            molar_mass: molar mass of substance [g/mol]
            pKa: disassociation constant for acic H-OMP [-]
            omp_half_life: per redox zone, [days]) 
        """
        self.name = name
    
        # Substance dict here as placeholder for the actual database
        substances_dict = { 
            'benzene': {
                'log_Koc': 1.92,
                'molar_mass': 78.1, 
                'pKa': 99,
                'omp_half_life': {
                    'suboxic': 10.5,
                    'anoxic': 420,
                    'deeply_anoxic': 1e99,
                    },
                },
            'AMPA': {
                'log_Koc': -0.36,
                'molar_mass': 111.04 , 
                'pKa': 0.4,
                'omp_half_life': {
                    'suboxic': 46,
                    'anoxic': 46,
                    'deeply_anoxic': 1e99,
                    },
                },
            'benzo(a)pyrene': {
                'log_Koc': 6.43,
                'molar_mass': 252.3, 
                'pKa': 99,
                'omp_half_life': {
                    'suboxic': 530,
                    'anoxic': 2120,
                    'deeply_anoxic': 2120,
                    },
                },
            }

        self.substance_dict = substances_dict[name]
        self.log_Koc = self.substance_dict['log_Koc']
        self.pKa = self.substance_dict['pKa']


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
    def __init__(self, schematisation, substance: Substance): #, substance, df_particle, df_flowline):
        # def __init__(self, substance: Substance, df_particle, df_flowline, removel_function?):
        self.omp_inialized = False
        self.df_particle = schematisation.df_particle
        self.df_flowline = schematisation.df_flowline
        self.substance = Substance(substance)

    def _init_omp(self):
        if self.omp_inialized:
            pass
        else:
            self.df_particle['omp_half_life'] = self.df_particle['redox_zone'].map(self.substance.substance_dict['omp_half_life'])
            self.df_particle['log_Koc'] = self.substance.log_Koc
            self.df_particle['pKa'] = self.substance.pKa
        self.omp_inialized = True


    def _init_microbiology():
        pass

    def _calculate_retardation(self):
        ''' Equation 4.8-4.10 in report
        Retardation equation based on Karickhoff (1981) and Schwarzenbach et al. (1993)
        (section 10.3 in Appelo & Postma 2005), however with addition of 
        the effects of (i) DOC-binding according to Kan & Tomson (1990),
        and (ii) OMP ionization (dissociation) according to Schellenberg et al. (1984)

        0.2 -> fraction of binding sites supplied by DOC which bind the OMP 
        and prevent sortion to aquifer'''




        self.df_particle['retardation'] = (1 + (1 / (1 + 10 ** (self.df_particle.pH - self.df_particle.pKa)) * self.df_particle.density_aquifer / 1000
                            * (1 - self.df_particle.porosity_layer)
                            * self.df_particle.fraction_organic_carbon * self.df_particle.Koc_temperature_correction)
                    / (self.df_particle.porosity_layer * (1 + (self.df_particle.Koc_temperature_correction * 1 / (1 + 10 ** (self.df_particle.pH - self.df_particle.pKa))
                                                            * 0.2 * self.df_particle.dissolved_organic_carbon * 0.000001))))

    def _calculate_omp_half_life_temperature_correction(self):
        '''Equation 3.2 in report
        R = 8.314 J/K/mol
        Ea = activation energy = 63*10^3 J/mol'''

        self.df_particle['omp_half_life_temperature_corrected'] = self.df_particle['omp_half_life'] * 10 ** (-63000 / (2.303 * 8.314) * (1 / (20 + 273.15) - 1 / (self.df_particle.temperature + 273.15)))
        self.df_particle.loc[ self.df_particle.omp_half_life == 1e99, 'omp_half_life_temperature_corrected'] = 1e99

    def _calculate_Koc_temperature_correction(self):
        ''' Equation 3.1 in report, 
        from Luers and Ten Hulscher (1996): Assuming the relation to be similar 
        to the Van ‘t Hoff equation and equally performing for other OMPs yields'''

        self.df_particle['Koc_temperature_correction'] = 10 ** self.df_particle.log_Koc * 10 ** (1913 * (1 / (self.df_particle.temperature + 273.15) - 1 / (20 + 273.15)))
    
    def _calculate_state_concentration_in_zone(self):
        '''Equation 4.11 in report '''
    
        for i in range(len(self.df_particle)-1):
            if self.df_particle.steady_state_concentration.loc[i+1] is None:
                
                # if omp is persistent, value at end of zone equal to value incoming to zone
                if self.df_particle.omp_half_life.loc[i+1] == 1e99:
                    self.df_particle.at[i+1, 'steady_state_concentration'] = self.df_particle.steady_state_concentration.loc[i]
                
                # AH @MartinvdS why 300 the limit here? just to avoid very small numnbers?
                # Column O in Phreatic excel sheet
                elif (self.df_particle.travel_time_zone.loc[i+1] * 365.25 * self.df_particle.retardation.loc[i+1]
                                                                            / self.df_particle.omp_half_life_temperature_corrected.loc[i+1]) >300:
                    self.df_particle.at[i+1, 'steady_state_concentration'] = 0

                # otherwise, calculate the outcoming concentration from the zone, given the input concentration to the zone. 
                # in the case of the vadose zone, the incoming concentration is the initial concentration
                else:
                    self.df_particle.at[i+1, 'steady_state_concentration'] = (self.df_particle.steady_state_concentration.loc[i]
                                                                            / (2 ** (self.df_particle.travel_time_zone.loc[i+1] * 365.25 * self.df_particle.retardation.loc[i+1]
                                                                            / self.df_particle.omp_half_life_temperature_corrected.loc[i+1])))

    def _calculcate_total_breakthrough_travel_time(self):
        ''' Calculate the total time (years) for breakthrough at the well'''
        self.df_flowline['total_breakthrough_travel_time']  = ""
        for i in range(len(self.df_flowline)):
            flowline_id = i + 1

            df = self.df_particle.loc[self.df_particle['flowline_id'] == flowline_id]
            df.fillna(0)['breakthrough_travel_time']
            self.df_flowline.at[i, 'total_breakthrough_travel_time'] = sum(df.fillna(0)['breakthrough_travel_time'])

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

        self._calculate_Koc_temperature_correction()
       
        self._calculate_omp_half_life_temperature_correction()

        self._calculate_retardation()

        self._calculate_state_concentration_in_zone()

        self.df_particle['breakthrough_travel_time'] = self.df_particle.retardation * self.df_particle.travel_time_zone

        self._calculcate_total_breakthrough_travel_time()
       
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

#%%

#   QUESITONS
# 2. input permeability, thickness seperately?
# 4. naming conventions: zone 1, zone 2, zone 3 or unsaturated zone, zone 1, zone 2 .. 
#   gets more complicated in the BAR (basin, phreatic/first/upper aquifer, deep/aquitard/second aquifer)
# 5. error in spreadsheet in "Peter's t2 fit", does not seem to take into account changes in the thickness of the aquifer and porosity (fixed n2 = 0.32 and D2 = 95)...?
# 6. Calculating the well_discharge_m3hour (q) for ech streamline, Phreatic %Qrx*Q, Semi-confinded: Qr/Qx*Qx?? see notes and df_output for each
# %%
