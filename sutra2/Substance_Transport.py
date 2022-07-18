#%% -----------------------------------------------------------------------------
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

#### Notes ####

# things which must be checked indicated in comments with AH
# specific questions flagged for;
# @MartinvdS // @steven //@martinK

####

#%% ----------------------------------------------------------------------------
# INITIALISATION OF PYTHON e.g. packages, etc.
# ------------------------------------------------------------------------------

# Plotting modules
import matplotlib.pyplot as plt
import matplotlib 
import matplotlib.colors as colors

import numpy as np
import pandas as pd
import os
import sys
from pandas import read_csv
from pandas import read_excel
import math
from scipy.special import kn as besselk
import datetime
from datetime import timedelta

# module_path = os.path.abspath(os.path.join("..","..","sutra2
#"))
# if module_path not in sys.path:
#     sys.path.insert(0,module_path)

from sutra2.Analytical_Well import AnalyticalWell 
from sutra2.ModPath_functions import ModPathWell

# from Analytical_Well import AnalyticalWell
# from ModPath_functions import ModPathWell
# # from sutra2
#.ModPath_functions import ModPathWell
# # from Analytical_Well import AnalyticalWell

path = os.getcwd()  # path of working directory

class Organism:
    ''' 
    Placeholder class which will later be replaced by the QSAR functionality of AquaPriori ?!

    For removal of microbial organisms ('mbo') / pathogen 
    
    microbial_species_dict: dict
    Attributes
    ---------
    organism_name: str
        species_name of the substance
        
    'alpha0': float
        reference_collision_efficiency [-]
        per redox zone ('suboxic', 'anoxic', deeply_anoxic')
    'pH0': float
        reference pH for calculating collision efficiency [-]
        per redox zone ('suboxic', 'anoxic', deeply_anoxic')
    'organism_diam': float
        diameter of pathogen/species [m]
    'mu1': float
        inactivation coefficient [1/day]
        per redox zone ('suboxic', 'anoxic', deeply_anoxic')
    '''  
    def __init__(self, organism_name, 
                    removal_function = 'mbo'):
        """
        Parameters
        ----------
        organism: str
            name of the organism (for now limited dictionary to 
            'solani','carotovorum', 'solanacearum')

        Returns
        -------
        organism_dict: dictionary
            'alpha0': float
                reference_collision_efficiency [-]
                per redox zone ('suboxic', 'anoxic', deeply_anoxic')
            'pH0': float
                reference pH for calculating collision efficiency [-]
                per redox zone ('suboxic', 'anoxic', deeply_anoxic')
            'organism_diam': float
                diameter of pathogen/species [m]
            'mu1': float
                inactivation coefficient [1/day]
                per redox zone ('suboxic', 'anoxic', deeply_anoxic')

        @SR (21-1-2022): Adjust documentation for microbial organisms (mbo)        
        """
        self.organism_name = organism_name

        # Dict    # Acties Steven 25-1-22

        # Naming convention organism: Uppercamelcase species
        micro_organism_dict = {
            "solani": 
                {"organism_name": "solani",
                    "alpha0": {
                        "suboxic": 0.037, 
                        "anoxic": 0.037e-2,     # NOT reported: factor 100 smaller than suboxic
                        "deeply_anoxic": 0.037e-2
                    },
                    "pH0": {
                        "suboxic": 7.5, 
                        "anoxic": 7.5,          # NOT reported: assumed equal to suboxic
                        "deeply_anoxic": 7.5    # NOT reported: assumed equal to suboxic
                    },
                    "organism_diam": 2.731e-6,
                    "mu1": {
                        "suboxic": 1.2472, 
                        "anoxic": 0.1151, 
                        "deeply_anoxic": 0.1151
                    }
                },
            "carotovorum": 
                {"organism_name": "carotovorum",
                    "alpha0": {
                        "suboxic": 0.300, 
                        "anoxic": 0.577, 
                        "deeply_anoxic": 0.577
                    },
                    "pH0": {
                        "suboxic": 7.5, 
                        "anoxic": 7.5, 
                        "deeply_anoxic": 7.5
                    },
                    "organism_diam": 1.803e-6,
                    "mu1": {
                        "suboxic": 1.2664, 
                        "anoxic": 0.1279, 
                        "deeply_anoxic": 0.1279
                    }
                },
            "solanacearum": 
                {"organism_name": "solanacearum",
                    "alpha0": {
                        "suboxic": 0.011, 
                        "anoxic": 0.456, 
                        "deeply_anoxic": 0.456
                    },
                    "pH0": {
                        "suboxic": 7.5, 
                        "anoxic": 7.5, 
                        "deeply_anoxic": 7.5
                    },
                    "organism_diam": 1.945e-6,
                    "mu1": {
                        "suboxic": 0.3519, 
                        "anoxic": 0.1637, 
                        "deeply_anoxic": 0.1637
                    }
                },
            }

        if self.organism_name in micro_organism_dict.keys():
            self.organism_dict = micro_organism_dict[self.organism_name]
        else: # return empty dict
            self.organism_dict = \
                {"organism_name": self.organism_name,
                 "alpha0": {
                    "suboxic": None, 
                    "anoxic": None, 
                    "deeply_anoxic": None
                    },
                 "pH0": {
                    "suboxic": None, 
                    "anoxic": None, 
                    "deeply_anoxic": None
                    },
                 "organism_diam": None,
                 "mu1": {
                    "suboxic": None, 
                    "anoxic": None, 
                    "deeply_anoxic": None
                     }
                }
        
class Substance:
    ''' 
    Placeholder class which will later be replaced by the QSAR functionality of AquaPriori.

    For removal of organic micro pollutants (omp)
    substances_dict: dict

    Attributes
    ---------
    substance_name: String,
        substance_name of the substance (for now limited dictionary to 'benzene', 'AMPA', 'benzo(a)pyrene'
    substance_dict: dictionary
        Nested dictionary with the following keys per substance.

    
    substance_name: String,
        substance_name of the substance (for now limited dictionary to 'benzene', 'AMPA', 'benzo(a)pyrene')
    log Koc: float
        distribution coefficient of organic carbon and water, [-]
    molar_mass: float
        molar mass of substance, [g/mol]
    pKa: float
        disassociation constant for acid H-OMP, [-]
    omp_half_life: float
        per redox zone ('suboxic', 'anoxic', deeply_anoxic'), [days]

    '''
    def __init__(self, substance_name, 
                    removal_function = 'omp'):
        """
        Parameters
        ----------
        substance_name: String,
            substance_name of the substance (for now limited dictionary to 'benzene', 'AMPA', 'benzo(a)pyrene'

        Returns
        -------
        substance_dict: dictionary
            log Koc: float
                distribution coefficient of organic carbon and water ([-]
            molar_mass: float
                molar mass of substance [g/mol]
            pKa: float
                disassociation constant for acic H-OMP [-]
            omp_half_life: float
                per redox zone, [days])
     
        """
        self.substance_name = substance_name

        # Substance dict here as placeholder for the actual database
        substances_dict = {
            'benzene': {
                'substance_name': 'benzene',
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
                'substance_name': 'AMPA',
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
                'substance_name': 'benzo(a)pyrene',
                'log_Koc': 6.43,
                'molar_mass': 252.3,
                'pKa': 99,
                'omp_half_life': {
                    'suboxic': 530,
                    'anoxic': 2120,
                    'deeply_anoxic': 2120,
                    },
                },
            'OMP-X': {
                'substance_name': 'OMP-X',
                'log_Koc': 0,
                'molar_mass': 100,
                'pKa': 99,
                'omp_half_life': {
                    'suboxic': 1e99,
                    'anoxic': 1e99,
                    'deeply_anoxic': 1e99,
                    },
                },
            }
        if self.substance_name in substances_dict.keys():
            self.substance_dict = substances_dict[self.substance_name]
        else: # return empty dict
            self.substance_dict = \
                {"substance_name": self.substance_name,
                 "log_Koc": 0,
                 "molar_mass": 0,
                 'pKa': None,
                 "omp_half_life": {
                    "suboxic": None, 
                    "anoxic": None, 
                    "deeply_anoxic": None
                     }
                }

#ah_todo @MartinK, MartinvdS -> let the user specify the chemical in the Substance transport file instead of schematisation?
# also let them feed it a dictionary with their own substance?


class SubstanceTransport():
    """ 
    Returns concentration in a groundwater well for a given Organic Micro Pollutant or microbial species.

    Attributes
    ----------
    well: object
        The AnalyticalWell object for the schematisation of the aquifer type.
    omp_inialized: bool
        Boolian indicating whether the Substance object has been initialized
    df_flowline: pandas.DataFrame
        Column 'flowline_id': int
        Column 'flowline_type': string
        Column 'flowline_discharge': float
        Column 'particle_release_day': float
        Column 'input_concentration': float
        Column 'endpoint_id': int
        Column 'well_discharge': float
        Column 'name': string
        Column 'removal_function': string
        Column 'total_breakthrough_travel_time': float
        Column 'breakthrough_concentration': float

    df_particle: pandas.DataFrame
        Column 'flowline_id': int
        Column 'zone': string
        Column 'travel_time': float
        Column 'xcoord': float
        Column 'ycoord': float
        Column 'zcoord': float
        Column 'redox': float
        Column 'temp_water': float
        Column 'travel_distance': float
        Column 'porosity': float
        Column 'dissolved_organic_carbon': float
        Column 'pH': float
        Column 'fraction_organic_carbon': float
        Column 'solid_density': float
        Column 'input_concentration': float
        Column 'steady_state_concentration': float
        Column 'omp_half_life': float
        Column 'log_Koc': float
        Column 'pKa': float
        Column 'Koc_temperature_correction': float
        Column 'omp_half_life_temperature_corrected': float
        Column 'retardation': float
        Column 'breakthrough_travel_time': float


    substance: object
        The Substance object with the OMP of interest.

    substance_dict: dictionary
        Nested dictionary with the following per substance.

        substance_name: String,
            substance_name of the substance (for now limited dictionary to 'benzene', 'AMPA', 'benzo(a)pyrene'
        log Koc: float
            distribution coefficient of organic carbon and water [-]
        molar_mass: float
            molar mass of substance [g/mol]
        pKa: float
            disassociation constant for acic H-OMP [-]
        omp_half_life: float
            per redox zone ('suboxic', 'anoxic', deeply_anoxic'), [days]

    organism: object
        The microbial organism object with the microbial organism (mbo) of interest

        organism_dict: dictionary
            'alpha0': float
                reference_collision_efficiency [-]
                per redox zone ('suboxic', 'anoxic', deeply_anoxic')
            'pH0': float
                reference pH for calculating collision efficiency [-]
                per redox zone ('suboxic', 'anoxic', deeply_anoxic')
            'organism_diam': float
                diameter of pathogen/species [m]
            'mu1': float
                inactivation coefficient [1/day]
                per redox zone ('suboxic', 'anoxic', deeply_anoxic')
    """


    def __init__(self,
                well: AnalyticalWell or ModPathWell,
                substance: Substance = 'benzo(a)pyrene',
                organism: Organism = 'solani',
                partition_coefficient_water_organic_carbon=None,
                dissociation_constant=None,
                halflife_suboxic=None,
                halflife_anoxic=None,
                halflife_deeply_anoxic=None,
                alpha0_suboxic=None,
                alpha0_anoxic=None,
                alpha0_deeply_anoxic=None,
                pH0_suboxic=None,
                pH0_anoxic=None,
                pH0_deeply_anoxic=None,
                mu1_suboxic=None,
                mu1_anoxic=None,
                mu1_deeply_anoxic=None,
                organism_diam=None,
                # Add 'removal_function'=None?
                removal_function=None
                ):


        '''
        Initialization of the SubstanceTransport class. Checks for either user-defined OMP substance parameters (removal_function == 'omp')
        or user-defined microbial organism removal parameters (removal_function == 'mbo') and overrides the database values, if any.

        Parameters
        ----------
        well: object
            The AnalyticalWell object for the schematisation of the aquifer type.
        substance: object
            The Substance object with the OMP of interest.
        organism: object
            The Organism object with microbial organism (MBO) of interest
        partition_coefficient_water_organic_carbon: float
            Distribution coefficient of OMP between organic carbon and water, dimensionless.
        dissociation_constant: float
            Dissociation equilibirum constant of the OMP, dimensionless.
        halflife_suboxic, halflife_anoxic, halflife_deeply_anoxic: float
            Time required to reduce the concentration of the OMP by half, from any concentration point in time [days].
        alpha0_suboxic, alpha0_anoxic, alpha0_deeply_anoxic: float
            reference_collision_efficiency [-]
            per redox zone ('suboxic', 'anoxic', deeply_anoxic')
        pH0_suboxic, pH0_anoxic, pH0_deeply_anoxic: float
            reference pH for calculating collision efficiency [-]
            per redox zone ('suboxic', 'anoxic', deeply_anoxic')
        mu1_suboxic, mu1_anoxic, mu1_deeply_anoxic: float
            inactivation coefficient [1/day]
            per redox zone ('suboxic', 'anoxic', deeply_anoxic')
        organism_diam: float
            diameter of pathogen/species [m]
        removal_function: str 
            removal_function: ['omp' or 'mbo']
            Calculate removal of either organic micro pollutants ('omp') or microbial organisms ('mbo')
        '''
        self.well = well
        self.df_particle = well.df_particle
        self.df_flowline = well.df_flowline

        # Substance
        self.substance_name = substance
        self.partition_coefficient_water_organic_carbon = partition_coefficient_water_organic_carbon
        self.dissociation_constant = dissociation_constant
        self.halflife_suboxic = halflife_suboxic
        self.halflife_anoxic = halflife_anoxic
        self.halflife_deeply_anoxic = halflife_deeply_anoxic
        
        # Organism
        self.organism_name = organism
        self.alpha0_suboxic=alpha0_suboxic
        self.alpha0_anoxic=alpha0_anoxic
        self.alpha0_deeply_anoxic=alpha0_deeply_anoxic
        self.pH0_suboxic=pH0_suboxic
        self.pH0_anoxic=pH0_anoxic
        self.pH0_deeply_anoxic=pH0_deeply_anoxic
        self.mu1_suboxic=mu1_suboxic
        self.mu1_anoxic=mu1_anoxic
        self.mu1_deeply_anoxic=mu1_deeply_anoxic
        self.organism_diam=organism_diam        

        # Removal function
        if removal_function is not None:
            self.removal_function = removal_function
        else:
            self.removal_function = well.schematisation.removal_function    

        # Load substance data
        self.substance = Substance(substance_name = substance, 
                                removal_function = self.removal_function)
        # Load microbial organism data
        self.organism = Organism(organism_name = organism, 
                                removal_function = self.removal_function)

        # Run init of 'omp'
        if self.removal_function == 'omp':
            self.omp_initialized = False
            self.micro_organism_initialized = True
      
        # Run init of 'pathogen'
        elif self.removal_function == 'mbo':
            self.omp_initialized = True
            self.micro_organism_initialized = False

        # Create user dict with 'removal_parameters' from input
        # if self.removal_function == 'omp':
        user_removal_parameters = {
            self.substance_name:
                {'substance_name': self.substance_name,
                'log_Koc': self.partition_coefficient_water_organic_carbon,
                'pKa': self.dissociation_constant,
                'omp_half_life': {
                    'suboxic': self.halflife_suboxic,
                    'anoxic': self.halflife_anoxic,
                    'deeply_anoxic': self.halflife_deeply_anoxic,
                    }
                },
            self.organism_name:
                {"organism_name": self.organism_name,
                    "alpha0": {
                        "suboxic": self.alpha0_suboxic, 
                        "anoxic": self.alpha0_anoxic, 
                        "deeply_anoxic": self.alpha0_deeply_anoxic
                        },
                    "pH0": {
                        "suboxic": self.pH0_suboxic, 
                        "anoxic": self.pH0_anoxic, 
                        "deeply_anoxic": self.pH0_deeply_anoxic
                        },
                    "organism_diam": self.organism_diam,
                    "mu1": {
                        "suboxic": self.mu1_suboxic, 
                        "anoxic": self.mu1_anoxic, 
                        "deeply_anoxic": self.mu1_deeply_anoxic
                        }
                    },
                }

        # Compare the removal_parameters dictionaries and override the default values if the user inputs a value
        
        if self.removal_function == 'omp':
            # User defined removal parameters [omp]
            user_removal_parameters = user_removal_parameters[self.substance_name]
            # assumes that default dict contains the substance input by the user (we only have three right now though!)
            default_removal_parameters = self.substance.substance_dict
        elif self.removal_function == 'mbo':
            # User defined removal parameters [omp]
            user_removal_parameters = user_removal_parameters[self.organism_name]
            # assumes that default dict contains microbial organism input 
            # (for now limited dictionary to 'solani','carotovorum', 'solanacearum')
            default_removal_parameters = self.organism.organism_dict

        # iterate through the dictionary keys
        # @Steven: check of over user_removal_parameters moet worden geïtereerd?! --> ipv over default_parameters
        for key, value in user_removal_parameters.items():
            if type(value) is dict:
                for tkey, cvalue in value.items():
                    if cvalue is None: #reassign the value from the default dict if not input by the user
                        user_removal_parameters[key][tkey] = default_removal_parameters[key][tkey]
                    elif type(cvalue) is dict:
                        for subkey, subval in cvalue.items():
                            if subval is None:
                                user_removal_parameters[key][tkey][subkey] = default_removal_parameters[key][tkey][subkey]
                    # else: no assignment from default dict required...
            else:
                if value is None:
                    user_removal_parameters[key] = default_removal_parameters[key]
            
        #assign updated dict as attribute of the class to be able to access later
        self.removal_parameters = user_removal_parameters 

        # microbial organism init
        # @SR : init of microbial organisms (pathogen)

        # Check if ModPathWell (Class) can be used instead of 'AnalyticalWell'
        
        # # AH need to make sure here that the substance passed is the same, e.g. comapre the dictionaries BUT ALSO
        # # make sure that user doesn't call one substance in the hydrochemicalschematisation class and another in the concentration class
        # # probably only a problem for ourselves, this should be written into a larger "run" class for the model which could avoid this

    def _init_omp(self):
        ''' 
        Initialisation if the Substance is an OMP
        '''
        if self.omp_initialized:
            pass
        else:
            self.df_particle['omp_half_life'] = self.df_particle['redox'].map(self.removal_parameters['omp_half_life'])
            self.df_particle['log_Koc'] = self.removal_parameters['log_Koc']
            self.df_particle['pKa'] = self.removal_parameters['pKa']

        self.omp_initialized = True


    def _init_micro_organism(self):
        ''' Initialisation if the Species is a microbial organism'''
        
        if self.micro_organism_initialized:
            pass
        else:
            self.df_particle['mu1'] = self.df_particle['redox'].map(self.removal_parameters['mu1'])
            self.df_particle['alpha0'] = self.df_particle['redox'].map(self.removal_parameters['alpha0'])
            self.df_particle['pH0'] = self.df_particle['redox'].map(self.removal_parameters['pH0'])
            self.df_particle['organism_diam'] = self.removal_parameters['organism_diam']

        self.micro_organism_initialized = True

    #@Steven_todo, @MartinK: include in a function (used for both omp_removal as well as mbo removal)     
    def _contamination_date_vs_well_abstraction(self):
        # reduce the amount of text per line by extracting the following parameters
        self.compute_contamination_for_date = self.well.schematisation.compute_contamination_for_date
        start_date_well = self.well.schematisation.start_date_well
        start_date_contamination = self.well.schematisation.start_date_contamination
        self.end_date_contamination = self.well.schematisation.end_date_contamination

        if start_date_well > start_date_contamination:
            self.start_date = start_date_well
            self.back_date_start = start_date_contamination

        elif start_date_well <= start_date_contamination:
            self.start_date = start_date_contamination
            self.back_date_start = start_date_well

        self.compute_date = self.compute_contamination_for_date - self.start_date
        self.back_compute_date = self.start_date - self.back_date_start
    
        # Particle release day
        self.particle_release_day = (self.start_date - start_date_contamination).days


    def _calculate_retardation(self):
        ''' Calculates the retardation of the OMP due to sorption and biodegradation.
        Adds a column to the 'df_particle' with the retardation value.

        Equation 4.8-4.10 in TRANSATOMIC report
        Retardation equation based on Karickhoff (1981) and Schwarzenbach et al. (1993)
        (section 10.3 in Appelo & Postma 2005), however with addition of
        the effects of (i) DOC-binding according to Kan & Tomson (1990),
        and (ii) OMP ionization (dissociation) according to Schellenberg et al. (1984)

        Returns
        -------
        df_particle: pandas.dataframe
            Column 'retardation': float

        '''
        #0.2 -> fraction of binding sites supplied by DOC which bind the OMP
        #and prevent sortion to aquifer
        if self.well.schematisation.biodegradation_sorbed_phase:
            self.df_particle['retardation'] = (1 + (1 / (1 + 10 ** (self.df_particle.pH - self.df_particle.pKa)) * self.df_particle.solid_density
                            * (1 - self.df_particle.porosity)
                            * self.df_particle.fraction_organic_carbon * self.df_particle.Koc_temperature_correction)
                    / (self.df_particle.porosity * (1 + (self.df_particle.Koc_temperature_correction * 1 / (1 + 10 ** (self.df_particle.pH - self.df_particle.pKa))
                                                            * 0.2 * self.df_particle.dissolved_organic_carbon * 0.000001))))
        else:
            self.df_particle['retardation'] = 1

    def _calculate_omp_half_life_temperature_correction(self):
        '''
        Corrects the OMP half-life for temperature if 'temp_correction_halflife' is 'True' in the HydroChemicalSchematisation.
        Adds column to 'df_particle' with corrected value.

        Equation 3.2 in TRANSATOMIC report
        R = 8.314 J/K/mol
        Ea = activation energy = 63*10^3 J/mol

        Returns
        -------
        df_particle: pandas.dataframe
            Column 'omp_half_life_temperature_corrected': float'''

        if self.well.schematisation.temp_correction_halflife:
            self.df_particle['omp_half_life_temperature_corrected'] = self.df_particle['omp_half_life'] * 10 ** (-63000 / (2.303 * 8.314) * (1 / (20 + 273.15) - 1 / (self.df_particle.temp_water + 273.15)))
        else:
            self.df_particle['omp_half_life_temperature_corrected'] = self.df_particle['omp_half_life']

        self.df_particle.loc[self.df_particle.omp_half_life == 1e99, 'omp_half_life_temperature_corrected'] = 1e99

    def _calculate_Koc_temperature_correction(self):
        ''' Corrects the OMP Koc for temperature if 'temp_correction_Koc' is 'True' in the HydroChemicalSchematisation.
        Adds column to 'df_particle' with corrected value.

        Equation 3.1 in TRANSATOMIC report,
        from Luers and Ten Hulscher (1996): Assuming the relation to be similar
        to the Van ‘t Hoff equation and equally performing for other OMPs yields

        Returns
        -------
        df_particle: pandas.dataframe
            Column 'Koc_temperature_correction': float
        '''


        # empty Koc_temperature_correction series to start with
        self.df_particle['Koc_temperature_correction'] = 0
        if self.well.schematisation.temp_correction_Koc:
            self.df_particle.loc[self.df_particle.log_Koc != 0,'Koc_temperature_correction'] = 10 ** self.df_particle.loc[self.df_particle.log_Koc != 0,"log_Koc"] * \
                                                                        10 ** (1913 * (1 / (self.df_particle.loc[self.df_particle.log_Koc != 0,"temp_water"] + 273.15) - 1 / (20 + 273.15)))
        else:
            self.df_particle.loc[:,'Koc_temperature_correction'] = self.df_particle.loc[:,"log_Koc"]
        
        # if log_Koc is zero, assign value of zero
        self.df_particle.loc[self.df_particle.log_Koc == 0,'Koc_temperature_correction'] = 0

    def _calculate_steady_state_concentration_in_zone_omp(self):
        '''
        Calculates the steady state concentration in the well for each flowline.
        Add column to 'df_particle' with the steady state concentration

        Equation 4.11 in TRANSATOMIC report

        Returns
        -------
        df_particle: pandas.dataframe
            Column 'steady_state_concentration': float

        '''

        #check if there is degradation prior to infiltration
        DOC_inf = self.well.schematisation.dissolved_organic_carbon_infiltration_water
        TOC_inf = self.well.schematisation.total_organic_carbon_infiltration_water

        if DOC_inf and TOC_inf > 0:
            DOC_TOC_ratio = DOC_inf / TOC_inf
            K_oc = self.df_particle['Koc_temperature_correction'].iloc[0]
            c_in = 100 - 100 * (1 - (DOC_TOC_ratio + (1 - DOC_TOC_ratio) / (1 + K_oc * TOC_inf  * 0.000001)))
            self.df_particle.loc[self.df_particle.zone=='surface', 'input_concentration']=c_in

        for i in range(len(self.df_particle)):
            if self.df_particle.steady_state_concentration.iloc[i] is None:

                # if omp is persistent, value at end of zone equal to value incoming to zone
                if (self.df_particle.omp_half_life.iloc[i] == 1e99) | np.isnan(self.df_particle.omp_half_life_temperature_corrected.iloc[i]):
                    self.df_particle.steady_state_concentration.iloc[i] = self.df_particle.steady_state_concentration.iloc[i-1]

                # Column O in Phreatic excel sheet
                # # AH 300 limit only to avoid very small numnbers, makes no difference for other calculations therefore left in
                # Put back in, otherwise there is an error there are too many numbers in the output
                # can't reproduce the error, so take out again
                
                elif (self.df_particle.travel_time.iloc[i] * self.df_particle.retardation.iloc[i]
                                                                            / self.df_particle.omp_half_life_temperature_corrected.iloc[i]) >300:
                    self.df_particle.steady_state_concentration.iloc[i] = 0

                # otherwise, calculate the outcoming concentration from the zone, given the input concentration to the zone.
                # in the case of the vadose zone, the incoming concentration is the initial concentration
                else:
                    self.df_particle.steady_state_concentration.iloc[i] = (self.df_particle.steady_state_concentration.iloc[i-1]
                                                                            / (2 ** (self.df_particle.travel_time.iloc[i] * self.df_particle.retardation.iloc[i]
                                                                            / self.df_particle.omp_half_life_temperature_corrected.iloc[i])))

    def _calculate_total_breakthrough_travel_time(self):
        ''' Calculate the total time for breakthrough for each flowline at the well

        Returns
        -------
        df_flowline: pandas.dataframe
            Column 'total_breakthrough_travel_time': float

        '''
        self.df_flowline['total_breakthrough_travel_time']  = 0
        self.df_flowline['breakthrough_concentration']  = 0.

        for fid in self.df_flowline.index: # range(len(self.df_flowline)):
            # flowline_id = i + 1

            # df = self.df_particle.loc[fid,:]
            # df.fillna(0)['breakthrough_travel_time']
            self.df_flowline.at[fid, 'total_breakthrough_travel_time'] = sum(self.df_particle.loc[fid,:].fillna(0)['breakthrough_travel_time'])
            self.df_flowline.at[fid, 'breakthrough_concentration'] = self.df_particle.loc[fid,'steady_state_concentration'].iloc[-1]

    def compute_omp_removal(self):
        """ 
        Calculates the concentration in the well of each flowline. Returns
        the values in 'df_flowline' and 'df_particle' as attributes of the object.

        Returns
        -------
        df_flowline: pandas.DataFrame
            Column 'flowline_id': Integer
            Column 'flowline_type': string
            Column 'flowline_discharge': Float
            Column 'particle_release_day': Float
            Column 'input_concentration': float
            Column 'endpoint_id': Integer
            Column 'well_discharge': float
            Column 'name': string
            Column 'removal_function': string
            Column 'total_breakthrough_travel_time': float
            The breakthrough concentration in the well for the OMP taking into account retardation.
            Column 'breakthrough_concentration': float
            The breakthrough concentration in the well for the OMP taking into account sorption
            and biodegradation.

        df_particle: pandas.DataFrame
            Column 'flowline_id': int
            Column 'zone': string
            Column 'travel_time': float
            Column 'xcoord': float
            Column 'ycoord': float
            Column 'zcoord': float
            Column 'redox': float
            Column 'temp_water': float
            Column 'travel_distance': float
            Column 'porosity': float
            Column 'dissolved_organic_carbon': float
            Column 'pH': float
            Column 'fraction_organic_carbon': float
            Column 'solid_density': float
            Column 'input_concentration': float
            Column 'steady_state_concentration': float
            The steady state concentration at the well of the OMP for the flowline, [mass/L]
            Column 'omp_half_life': float
            Column 'log_Koc': float
            Column 'pKa': float
            Column 'Koc_temperature_correction': float
            The temperature corrected Koc value, only if 'temp_correction_Koc' is 'True' in the HydroChemicalSchematisation.
            Column 'omp_half_life_temperature_corrected': float
            The temperature corrected OMP half-life value, if 'temp_correction_halflife' is 'True' in the HydroChemicalSchematisation.
            Column 'retardation': float
            Column 'breakthrough_travel_time': float

            """

        # load (diffuse) input concentration    
        self.df_flowline.loc[:,'input_concentration'] = self.well.schematisation.diffuse_input_concentration
        self.df_particle.loc[:,'input_concentration'] = None
        if 'steady_state_concentration' not in self.df_particle.columns:
            # Steady state concentration [-]
            self.df_particle.loc[:,'steady_state_concentration'] = None

        # First value of input_concentration/steady_state_concentration equals diffuse input concentration of pathline
        self.df_particle.loc[self.df_particle.total_travel_time == 0.,"input_concentration"] = self.df_flowline.loc[:,'input_concentration'].values
        self.df_particle.loc[self.df_particle.total_travel_time == 0.,"steady_state_concentration"] = self.df_flowline.loc[:,'input_concentration'].values


        if self.well.schematisation.point_input_concentration:
            ''' point contamination '''
            # need to take into account the depth of the point contamination here....
            # need to change the df_particle and df_flowline to only be the flowlines for the point contamination flowline(s)
            # use a single point contamination for now
            # FIRST recalculate the travel times for the contamination, then initialize the class

            #only for a SINGLE point contamination
            distance = self.well.schematisation.distance_point_contamination_from_well
            depth = self.well.schematisation.depth_point_contamination
            cumulative_fraction_abstracted_water = (math.pi * self.well.schematisation.recharge_rate
                                                        * distance ** 2)/abs(self.well.schematisation.well_discharge)
            ind = self.df_particle.flowline_id.iat[-1]

            if self.well.schematisation.schematisation_type == 'phreatic':
                head = self.well.schematisation._calculate_hydraulic_head_phreatic(distance=distance)
                df_flowline_points, df_particle_points = self.well._add_phreatic_point_sources(distance=distance,
                                            depth_point_contamination=depth,
                                            cumulative_fraction_abstracted_water=cumulative_fraction_abstracted_water)

            elif self.well.schematisation.schematisation_type == 'semiconfined':
                bottom_vadose_zone = self.well.schematisation.bottom_vadose_zone_at_boundary

                df_flowline_points, df_particle_points = self.well._add_semiconfined_point_sources(distance=distance,
                                        depth_point_contamination=depth,  )
            # Assign flowline id to point
            df_particle_points.loc[:,'flowline_id'] = df_particle_points.loc[:,'flowline_id'].values + ind
            df_flowline_points.loc[:,'flowline_id'] = df_flowline_points.loc[:,'flowline_id'].values + ind
            # reindex df_particle and df_flowline
            df_particle_points.index = df_particle_points.flowline_id
            df_flowline_points.index = df_flowline_points.flowline_id

            # Add (input) concentration to point source df's
            df_flowline_points.loc[:,'input_concentration'] = self.well.schematisation.point_input_concentration
            df_particle_points.loc[:,'input_concentration'] = None
            df_particle_points.loc[:,'steady_state_concentration'] = None

            # First value of input_concentration/steady_state_concentration equals point input concentration of pathline
            # for fid in df_particle_points.flowline_id.unique():
                
            df_particle_points.loc[df_particle_points.total_travel_time == 0., 'input_concentration'] = df_flowline_points.loc[:,'input_concentration'].values
            df_particle_points.loc[df_particle_points.total_travel_time == 0., 'steady_state_concentration'] = df_flowline_points.loc[:,'input_concentration'].values

            # Add flow/discharge data to point particle df's
            df_flowline_points.loc[:,'flowline_type'] = "point_source"
            df_flowline_points.loc[:,'flowline_discharge'] = abs(self.well.schematisation.discharge_point_contamination)

            #AH_todo, something here to loop through different point sources if more than one.

            # Add point particles to self.df_particle and self.df_flowline & reindex based on flowline_id
            self.df_particle = self.df_particle.append(df_particle_points)
            self.df_particle.index = self.df_particle.loc[:,"flowline_id"].values
            # self.df_particle.reset_index(drop=True, inplace=True)

            self.df_flowline = self.df_flowline.append(df_flowline_points)
            self.df_flowline.index = self.df_flowline.loc[:,"flowline_id"].values
            # self.df_flowline.reset_index(drop=True, inplace=True)

            # Add substance name
            self.df_flowline.loc[:,'name'] = self.removal_parameters['substance_name']

        self._init_omp()

        self._calculate_Koc_temperature_correction()

        self._calculate_omp_half_life_temperature_correction()

        self._calculate_retardation()

        self._calculate_steady_state_concentration_in_zone_omp()

        self.df_particle.loc[:,'breakthrough_travel_time'] = self.df_particle.loc[:,"retardation"].values * self.df_particle.loc[:,"total_travel_time"].values

        self._calculate_total_breakthrough_travel_time()

        # Calculate contamination relationship to well abstraction date
        self._contamination_date_vs_well_abstraction()

        # add the particle release date
        self.df_flowline.loc[:,'particle_release_day'] = self.particle_release_day


    def compute_concentration_in_well_at_date(self):
        #@Martink, this function is quite slow. I'm not sure how to make it go faster?
        ''' 
        Calculates the concentration in the well up to a specific date,
        taking into account the start and end date of the contamiantion and
        start date of the well.

        Returns
        -------
        df_well_concentration: pandas.dataframe
            Column 'time': float
                Array of time starting at minimum of the 'start_date_well' or 'start_date_contamination' as time = 0
                 and the other value is set as negative.
            Column 'date': datetime.date
                Array of dates starting at the minimum of the 'start_date_well' or 'start_date_contamination'
            Column total_concentration_in_well: float
                Summed concentration of the OMP in the well.
        '''

        # calculate the time after which no more contamination
        if self.end_date_contamination is None:
            pass
        else:
            end_time = self.end_date_contamination- self.start_date
            self.df_flowline['end_time_contamination_breakthrough'] = self.df_flowline['total_breakthrough_travel_time'] + end_time.days

        # AH_todo, Solution to make this faster is to erduce this down to xx timesteps, e.g. once a month
        time_array = np.arange(0, self.compute_date.days+1, 1)
        back_date_array = np.arange(-self.back_compute_date.days,0, 1)
        time_array = np.append(back_date_array,time_array)
        time_array_dates = pd.date_range(start=self.back_date_start,end=self.compute_contamination_for_date)

        #Calculate the concentration in the well,
        self.df_flowline['concentration_in_well'] = (self.df_flowline['breakthrough_concentration']
                            * self.df_flowline['flowline_discharge']/ self.df_flowline['well_discharge'])
        df_flowline = self.df_flowline

        # AH_todo, preset the length of list to improve some time
        well_concentration = []

        #sum the concentration in the well for each timestep
        for i in range(len(time_array)):
            t = time_array[i]
            if self.end_date_contamination is None:
                well_conc = sum(df_flowline['concentration_in_well'].loc[(df_flowline['total_breakthrough_travel_time'] <= t)])
            else:
                well_conc = sum(df_flowline['concentration_in_well'].loc[(df_flowline['total_breakthrough_travel_time'] <= t) & (df_flowline['end_time_contamination_breakthrough'] >= t)])
            well_concentration.append(well_conc)
        df_well_concentration = pd.DataFrame({'time':time_array, 'date':time_array_dates, 'total_concentration_in_well': well_concentration})

        return df_well_concentration

    def plot_concentration(self,
                            xlim=None,
                            ylim=None,
                            as_fraction_input = None,
                            x_axis = 'Date'):
        ## @Alex: adjust name of function to 'plot_concentration_timeseries'
        ''' Plot the concentration of the given OMP as a function of time since the start of the contamination

        Parameters
        ----------
        x_axis: string
            Choice of the x-axis as 'Time' in years starting at 0, or as the 'Date' since the
            minimum of 'start_date_well' or 'start_date_contamination'
        as_fraction_input: bool
            If 'True' plots concentration on y-axis as a fraction of the sum of the
            input concentration (diffuse and point source), [C/C0]
        xlim: array
            The x-axis limits
        ylim: array
            The y-axis limits

        ReturnVage @MartinK what to put here?
        '''



        # reduce the amount of text per line by extracting the following parameters
        point_input_concentration = self.well.schematisation.point_input_concentration
        diffuse_input_concentration = self.well.schematisation.diffuse_input_concentration
        schematisation_type = self.well.schematisation.schematisation_type
        compute_contamination_for_date = self.well.schematisation.compute_contamination_for_date
        start_date_well = self.well.schematisation.start_date_well
        start_date_contamination = self.well.schematisation.start_date_contamination
        end_date_contamination = self.well.schematisation.end_date_contamination

        start_date = max(start_date_well,start_date_contamination)
        back_date_start = min(start_date_well,start_date_contamination)
        compute_date = compute_contamination_for_date - start_date

        if point_input_concentration is None:
            input_concentration = diffuse_input_concentration
        else:
            input_concentration = diffuse_input_concentration + point_input_concentration

        df_well_concentration = self.compute_concentration_in_well_at_date()

        # as fraction of the input concentration
        if as_fraction_input:
            df_well_concentration[:] = [x / input_concentration for x in df_well_concentration]
            ylabel = 'Fraction of input concentration'
        else:
            ylabel = 'Concentration (ug/L)'

        fig = plt.figure(figsize=[10, 5])
        if x_axis == 'Date':
            plt.plot(df_well_concentration.date, df_well_concentration.total_concentration_in_well, 'b', label =str(self.substance.substance_name))
            plt.axvline(x=start_date_well, color= 'k', label = 'Start date well')
            plt.axvline(x=start_date_contamination, color= 'r', label = 'Start date contamination')
            if end_date_contamination is None:
                pass
            else:
                plt.axvline(x=end_date_contamination, color= 'g', label = 'End date contamination')

            plt.xlabel('Date')
            if xlim == None:
                plt.xlim([datetime.date((back_date_start.year-5), 1, 1), compute_contamination_for_date])
            else: plt.xlim(xlim)

        elif x_axis == 'Time':
            plt.plot(df_well_concentration.time/365.24, df_well_concentration.total_concentration_in_well,  'b', label =str(self.substance.substance_name))
            plt.axvline(x=(start_date_well-start_date).days/365.24, color= 'k', label = 'Start date well')
            plt.axvline(x=(start_date_contamination-start_date).days/365.24, color= 'r', label = 'Start date contamination')
            if end_date_contamination is None:
                pass
            else:
                plt.axvline(x=(end_date_contamination-start_date).days/365.24, color= 'g', label = 'End date contamination')

            plt.xlabel('Time (years)')
            if xlim == None:
                plt.xlim([(back_date_start-start_date).days/365-5,compute_date.days/365.24])
            else: plt.xlim(xlim)

        if ylim == None:
            plt.ylim([0, input_concentration])
        else: plt.ylim(ylim)
        plt.legend(loc=2)
        plt.ylabel(ylabel)
        plt.title('Aquifer type: ' + schematisation_type)
        plt.grid()
        plt.legend()
        # plt.savefig('well_concentration_over_time_'+str(self.substance.substance_name)+'_'+schematisation_type+'.png', dpi=300, bbox_inches='tight')
        return fig

    # Check for parameters in df_flowline #
    def _df_fillna(self, df,df_column: str, value = 0., dtype_ = 'float'):
        ''' Check dataframe for missing values for
            calculation of removal.
            df: the pandas.DataFrame to check
            df_column: the required dataframe column
            value: the default value in case other alues are missing
            dtype_: dtype of df column
            Return adjusted dataframe 'df'
        '''

        if not df_column in df.columns:
            # Add dataframe column with default value
            df[df_column] = value
        else:
            # Fill dataframe series (if needed with default value)
            if df[df_column].dropna().empty:
                df[df_column] = df[df_column].fillna(value)
            else: 
                # fill empty rows (with mean value of other records)
                if not dtype_ == 'object':
                    value_mean = df[df_column].values.mean()
                    df[df_column] = df[df_column].fillna(value_mean)
                else: 
                    df[df_column] = df[df_column].fillna(value)
        # Adjust dtype
        df = df.astype({df_column: dtype_})

        return df

    def calc_lambda(self, df_particle, df_flowline, 
                redox = 'anoxic',
                mu1 = 0.149,
                por_eff = 0.33, 
                grainsize = 0.00025,
                alpha0 = 0.001,
                pH = 7.5,
                pH0 = 7.5,
                temp_water = 10., rho_water = 999.703,
                organism_diam = 2.33e-8, v_por = 0.01):

        ''' For more information about the advective microbial removal calculation: 
            BTO2012.015: Ch 6.7 (page 71-74)

            Calculate removal coefficient lambda [/day].
            
            Parameters
            -----------
            redox: str
                redox condition ['suboxic','anoxic','deeply_anoxic']
            
            mu1: float
                inactivation coefficient [day-1]
            
            por_eff: float
                effective porosity [-]
            
            grainsize: float
                grain diameter of sediment [m]

            alpha0: float
                'reference sticky coefficient', for a reference pH [pH0]

            pH: float
                pH of the water [-]

            pH0: float
                reference pH for which alpha0 was determined 
            
            temp_water: float
                Water temperature [degrees celcius]
            
            rho_water: float
                Water density [kg m-3]
            
            alpha: float
                'sticky coefficient' [-], pH corrected
            
            organism_diam: float
                Organism/species diameter [m]

            v_por: float
                porewater velocity [m/d]

            const_BM: float
                Boltzmann constant [1,38 × 10-23 J K-1] 
        
            Return
            ---------
            df_flowline: pandas.DataFrame
                Column "organism_diam" for each node

            df_particle: pandas.DataFrame
                Column "relative_distance" for each node
                Column "porewater_velocity" for each node
                Column "k_att" for each node
                Column "lamda" for each node
                Column "alpha" for each node
        '''

        # Boltzmann coefficient [J K-1]
        const_BM = 1.38e-23    

        # Create empty column 'relative_distance' in df_particle
        df_particle['relative_distance'] = None  
        # Create empty column 'porewater_velocity' in df_particle
        df_particle['porewater_velocity'] = None     


        # Pathogen diameter [m]
        df_flowline = self._df_fillna(df_flowline,
                                df_column = 'organism_diam',
                                value = organism_diam)

        # Water temperature [degrees Celsius]
        df_particle = self._df_fillna(df_particle,
                                df_column = 'temp_water',
                                value = temp_water)

        # Water density [kg m-3]
        df_particle = self._df_fillna(df_particle,
                                df_column = 'rho_water',
                                value = rho_water) 
        
        # Effective porosity [-]
        df_particle = self._df_fillna(df_particle,
                                df_column = 'porosity',
                                value = por_eff)   

        # grain size [m]
        df_particle = self._df_fillna(df_particle,
                                df_column = 'grainsize',
                                value = grainsize)

        # pH of groundwater [-]
        df_particle = self._df_fillna(df_particle,
                                df_column = 'pH',
                                value = pH)

        # redox condition ['suboxic','anoxic','deeply_anoxic']
        df_particle = self._df_fillna(df_particle,
                                df_column = 'redox',
                                value = redox, dtype_ = 'object')

        # alpha0 [-]
        df_particle = self._df_fillna(df_particle,
                                df_column = 'alpha0',
                                value = alpha0)
        # reference pH [-]
        df_particle = self._df_fillna(df_particle,
                                df_column = 'pH0',
                                value = pH0)

        # mu1 [day -1]
        df_particle = self._df_fillna(df_particle,
                                df_column = 'mu1',
                                value = mu1)

        # Create empty column 'alpha' in df_particle
        df_particle['alpha'] = None  

        # Calculate 'sticky coefficient' alpha [np.array]
        alpha = {}
        for pid in df_flowline.index:
            alpha[pid] = df_particle.loc[pid,"alpha0"].values * 0.9**((df_particle.loc[pid,"pH"].values - df_particle.loc[pid,"pH0"].values )/0.1)

            # Fill df_particle 'alpha'
            df_particle.loc[pid,"alpha"] = alpha[pid]

        # Create empty column 'k_att' in df_particle
        df_particle['k_att'] = None
        # Create empty column 'lambda' in df_particle
        df_particle['lamda'] = None        
        
        # Calculate porewater velocity
        v_por, dist, tdiff = {}, {}, {}   
        for pid in df_flowline.index:

            # Distance squared argument
            dist_ = (df_particle.loc[pid,"xcoord"].diff().values**2 + 
                                df_particle.loc[pid,"ycoord"].diff().values**2 + 
                                df_particle.loc[pid,"zcoord"].diff().values**2)
            # Replace nan values
            dist_ = np.nan_to_num(dist_.astype('float'),nan=0.)
            #Distance array between nodes
            dist[pid] = np.sqrt(dist_)
            # Time difference array
            tdiff[pid] = np.nan_to_num(df_particle.loc[pid,"total_travel_time"].diff().values.astype('float'),0.001)

            # Calculate porewater velocity [m/day] (do not include effective porosity)
            v_por[pid] = abs(dist[pid][1:]/tdiff[pid][1:])
            # correct "0" velocity values
            v_por[pid][v_por[pid] == 0.] = 0.001

            # Fill column relative_distance in 'df_particle' 
            df_particle.loc[pid,'relative_distance'] = np.array([0] + list(dist[pid][1:]))
            # df_particle.loc[pid,'relative_distance'][-1] = 0

            # Fill porewater velocity in 'df_particle'
            df_particle.loc[pid,"porewater_velocity"] = np.array([v_por[pid][0]] + list(v_por[pid]))
            # In the last pathline row the porewater_velocity is equal to previous velocity
            # df_particle.loc[pid,"porewater_velocity"][-1] = v_por[pid][-1]

            # number of nodes
            n_nodes = len(df_particle.loc[pid,:])            

        # Create empty dicts and keys to keep track of arrays per particle id
        k_coll, gamma, As_happ, mu = {}, {}, {}, {}
        D_BM, k_diff, k_att, lamda = {}, {}, {}, {}
        # particle ids
        for pid in df_flowline.index:

            # Collision term 'k_coll'
            k_coll[pid] = (3/2.)*((1-df_particle.loc[pid,"porosity"].values) / \
                        df_particle.loc[pid,"grainsize"].values) * df_particle.loc[pid,"alpha"].values

            # Porosity dependent variable 'gamma'
            gamma[pid] = (1-df_particle.loc[pid,"porosity"].values)**(1/3)

            # Calculate Happel’s porosity dependent parameter 'A_s' (Eq. 5: BTO2012.015)
            ''' !!! Use correct formula:-> As =  2 * (1-gamma**5) /  (2 - 3 * gamma + 3 * gamma**5 - 2 * gamma**6)
                instead of... 2 * (1-gamma)**5 / (.......) 
            '''
            As_happ[pid] = 2 * (1-gamma[pid]**5) / \
                    (2 - 3 * gamma[pid] + 3 * gamma[pid]**5 - 2 * gamma[pid]**6)

            # Dynamic viscosity (mu) [kg m-1 s-1]
            mu[pid] = (df_particle.loc[pid,"rho_water"].values * 497.e-6) / \
                        (df_particle.loc[pid,"temp_water"].values + 42.5)**(3/2)
 
            # Diffusion constant 'D_BM' (Eq.6: BTO2012.015) --> unit: m2 s-1
            D_BM[pid] = (const_BM * (df_particle.loc[pid,"temp_water"].values + 273.)) / \
                        (3 * np.pi * df_flowline.loc[pid,"organism_diam"] * mu[pid])
            # Diffusion constant 'D_BM' (Eq.6: BTO2012.015) --> unit: m2 d-1
            D_BM[pid] *= 86400.

            # Diffusion related attachment term 'k_diff'
            k_diff[pid] = ((D_BM[pid] /
                        (df_particle.loc[pid,"grainsize"].values * df_particle.loc[pid,"porosity"].values * df_particle.loc[pid,"porewater_velocity"].values))**(2/3) * 
                            df_particle.loc[pid,"porewater_velocity"].values)

            # 'attachment coefficient' k_att [/dag]
            k_att[pid] = k_coll[pid] * 4 * As_happ[pid]**(1/3) * k_diff[pid]
            # removal coefficient 'lamda' [/day], using the 'mu1' mean.
            lamda[pid] = k_att[pid] + mu1

            # Fill df_particle 'k_att'
            df_particle.loc[pid,"k_att"] = k_att[pid]
            # Fill df_particle 'lambda'
            df_particle.loc[pid,"lamda"] = lamda[pid]

        # return (adjusted) df_particle and df_flowline
        return df_particle, df_flowline
        
    def calc_advective_microbial_removal(self,df_particle,df_flowline, 
                                        endpoint_id, trackingdirection = "forward",
                                        grainsize = 0.00025,
                                        temp_water = 11., rho_water = 999.703,
                                        redox = 'anoxic', por_eff = 0.33,
                                        pH = 7.5,
                                        organism_name = None,
                                        organism_diam = 2.33e-8,
                                        mu1 = 0.023, alpha0 = 1.E-5, pH0 = 6.8,
                                        conc_start = 1., conc_gw = 0.):
                                        
        ''' 
            Calculate the steady state concentration along traveled distance per 
            node for each pathline from startpoint to endpoint_id' after advective microbial removal.

            For more information about the advective microbial removal calculation: 
                BTO2012.015: Ch 6.7 (page 71-74)

            #SR_todo: remove df_flowline 'substance', 'removal_function', 'input_concentration' [species independent vars]
            #SR todo: remove df_particle 'lamda' 
            # --> not input, but output (replace 'substance' by 'name' of OMP or MBO)
            
            Attributes
            -----------
            df_particle: pandas.DataFrame
                Column 'flowline_id': int
                Column 'xcoord': float
                Column 'ycoord': float
                Column 'zcoord': float
                Column 'relative_distance': float
                Column 'travel_time': float
                Column 'total_travel_time': float
                Column 'porosity': float
                Column 'porewater_velocity': float

            df_flowline: pandas.DataFrame
                Column 'flowline_id': int or str
                    Identifier of pathline
                Column 'flowline_type': str
                    Source of contamination ['diffuse_source','point_source',]
                Column 'flowline_discharge': float
                    Flux per flowline (='pathline') [m^3 d-1]
                Column 'particle_release_day': float or int
                    Delay in release of contaminants [days]
                Column 'endpoint_id': int or str
                    Endpoint id of the flowline
                Column 'well_discharge': float
                    Cumulative flux of flowlines to well [m^3 d-1]

            endpoint_id: str or int
                ID for which to calculate the final concentration.

            tracking_direction: str
                tracking direction of particles ['forward','backward']

            grainsize: float
                grain diameter of sediment [m]

            temp_water: float
                Water temperature [degrees celcius]

            rho_water: float
                Water density [kg m-3]

            redox: str
                redox condition ['suboxic','anoxic','deeply_anoxic']

            por_eff: float
                effective porosity [-]

            pH: float
                pH of the water [-]
            
            mu1: float
                inactivation coefficient [day-1]
            
            pH0: float
                reference pH for which alpha0 was determined 
            
            alpha: float
                'sticky coefficient' [-], pH corrected
            
            alpha0: float
                'reference sticky coefficient', for a reference pH [pH0]
            
            organism_diam: float
                organism/species diameter [m]
            
            v_por: float
                porewater velocity [m/d]
            
            conc_start: float
                starting concentration
            
            conc_gw: float
                initial groundwater concentration
            
            distance_traveled: float
                distance between points [m]
            
            traveltime: float
                time between start and endpoint [days]

            Returns
            --------
            df_flowline: pandas.DataFrame
                Column "name": str
                    organism name of organism to calculate removal for
                Column "organism_diam": str
                    organism diam of organism to calculate removal for
                Column "removal_function": str
                    filled with 'mbo' -> removal function for microbial organisms
                    other option -> 'omp' -> for removal of organic micropollutants
                Column "input_concentration": float
                    Starting concentration of OMP or mbo of interest [N/L]
                Column "input_concentration_gw": float
                    Starting concentration in ambient groundwater [N/L]
                Column "breakthrough_concentration": float
                    Breakthrough_concentration over flowline (after reaching endpoint_id) [N/L]
                Column "breakthrough_travel_time": float
                    Cumulative traveltime per flowline from start to endpoint [days]
                Column "concentration_in_well": float
                    Steady-state concentration in well for endpoint_id
                    (assuming all flowlines were contaminated a traveltime duration ago)

            df_particle: pandas.DataFrame
                Column "relative_distance": float
                    distance between previous and current node
                Column "porewater_velocity": float
                    groundwater velocity between previous and current node
                Column "k_att": float
                    'attachment coefficient' [day-1]
                Column "lamda": float
                    'removal rate' [day-1] (redox dependent) --> calculated
                Column "alpha": float
                    'sticky coefficient' [-], pH corrected
                Column "steady_state_concentration": float
                    Steady-state concentration at each node, calculated [N/L]

        '''
        if organism_name is None:
            organism_name = self.organism_name
        
        # Organism/Species name
        self.df_flowline['name'] = organism_name

        # Initialize mbo-parameters
        self._init_micro_organism()



        # Calculate removal coefficient 'lamda' [/day].
        df_particle, df_flowline = self.calc_lambda(df_particle, df_flowline,
                                                    redox = redox, mu1 = mu1, por_eff = por_eff,
                                                    grainsize = grainsize, alpha0 = alpha0, 
                                                    pH = pH, pH0 = pH0, temp_water = temp_water, 
                                                    rho_water = rho_water, organism_diam = organism_diam)

        # Starting concentration [-]
        df_flowline = self._df_fillna(df_flowline,
                                df_column = 'input_concentration',
                                value = conc_start,
                                dtype_ = 'float')

        # Original concentration in groundwater [-]
        df_flowline = self._df_fillna(df_flowline,
                                df_column = 'input_concentration_gw',
                                value = conc_gw,
                                dtype_ = 'float')

        if 'breakthrough_concentration' not in df_flowline.columns:
            # breakthrough concentration for particles along a certain flowline
            df_flowline['breakthrough_concentration'] = None

        if 'breakthrough_traveltime' not in df_flowline.columns:
            # breakthrough traveltime for particles along a certain flowline
            df_flowline['breakthrough_travel_time'] = None

        if 'concentration_in_well' not in df_particle.columns:
            # Final (averaged) steady state concentration in well at endpoint_id [N/L]
            df_particle['concentration_in_well'] = None

        if 'steady_state_concentration' not in df_particle.columns:
            # steady-state concentration for particles along a certain flowline
            df_particle['steady_state_concentration'] = None
        
        for pid in df_flowline.index:
            df_particle.loc[pid,"steady_state_concentration"].iloc[0] = df_flowline.loc[pid,"input_concentration"]

        # Start dictionaries for relative and steady-state concentrations
        conc_rel = {}
        conc_steady = {}  # concentration corrected for initial concentration
        for pid in df_flowline.index:
            
            # Calculate the relative removal along pathlines
            exp_arg = -((df_particle.loc[pid,"lamda"].values/df_particle.loc[pid,"porewater_velocity"].values) *
                                df_particle.loc[pid,'relative_distance'].values).astype('float')

            conc_rel[pid] = df_flowline.loc[pid,"input_concentration_gw"] + \
                            (df_flowline.loc[pid,"input_concentration"] - df_flowline.loc[pid,"input_concentration_gw"]) * \
                                np.exp(exp_arg)

            if trackingdirection == 'forward':
                conc_steady[pid] = (conc_rel[pid]/df_flowline.loc[pid,"input_concentration"]).cumprod() * \
                                    df_flowline.loc[pid,"input_concentration"]
                # Replace 'nan'-values by 0.
                conc_steady[pid]  = np.nan_to_num(conc_steady[pid],0.)
                
                # Index of average final concentration "C_"  idx=-1[C0,...,..,...,Ct-1,C_final]
                idx_final = -1
                                
            elif trackingdirection == 'backward':
                conc_steady[pid] = (conc_rel[pid][::-1]/df_flowline.loc[pid,"input_concentration"]).cumprod() * \
                                    df_flowline.loc[pid,"input_concentration"]
                # Replace 'nan'-values by 0.
                conc_steady[pid]  = np.nan_to_num(conc_steady[pid],0.)
                # First value represents concentration of endpoint and v.v.
                # Reverse array--> [C_final,Ct-1,..,...,C0]
                conc_steady[pid] = conc_steady[pid][::-1]
                # Index of average final concentration "C_" idx=0 [C_final,Ct-1,..,...,C0]
                idx_final = 0

            # Add steady_state_concentration to df_particle
            df_particle.loc[pid,'steady_state_concentration'] = conc_steady[pid]

        # Calculate average final concentration
        C_final = 0.
        for pid in df_flowline.index:
            if df_flowline.loc[pid,"endpoint_id"] == endpoint_id:
                # Steady-state concentration after breakthrough
                C_breakthrough = conc_steady[pid][idx_final]
                # Traveltime from contamination to endpoint location [days]
                time_breakthrough = df_particle.loc[pid,"total_travel_time"].iloc[-1]

                # Steady-state final concentration
                C_final += (C_breakthrough * df_flowline.loc[pid,"flowline_discharge"]) / \
                            df_flowline.loc[pid,"well_discharge"]

                # Add breakthrough_concentration to df_flowline
                df_flowline.loc[pid,'breakthrough_concentration'] = C_breakthrough
                # Add breakthrough_traveltime to df_flowline
                df_flowline.loc[pid,'breakthrough_travel_time'] = time_breakthrough
        
        # Add final concentration in well (at endpoint_id)
        df_flowline.loc[df_flowline.endpoint_id == endpoint_id,"concentration_in_well"] = C_final

        # SR_todo: @MartinvdS required to calc total_breakthrough_travel_time mbo?
        # self._calculate_total_breakthrough_travel_time()

        self._contamination_date_vs_well_abstraction()
        
        # add the particle release date
        self.df_flowline['particle_release_day'] = self.particle_release_day # (self.start_date - start_date_contamination).days

        # return (adjusted) df_particle and df_flowline
        return df_particle, df_flowline, C_final

    # Plot the distribution of traveltimes from the well to model radius
    def plot_travel_time_distribution(self, df_particle, times_col = "total_travel_time",
                                      distance_col = "xcoord", index_col = "flowline_id",
                                      fpath_fig = None):
        ''' Plot cumululative travel times using column 'times_col' relative to distance 'distance_col' 
            per pathline with index column 'index_col'.'''

        # Pathline indices
        flowline_ids = list(df_particle.loc[:,index_col].unique())
        
        ## modpath solution ##
        # Distance points
        distance_points = np.empty((len(flowline_ids)), dtype = 'float')
        # Time points
        time_points = np.empty((len(flowline_ids)), dtype = 'float')

        for idx, fid in enumerate(flowline_ids):

            # Fill distance_points array
            distance_points[idx] = df_particle.loc[df_particle[index_col] == fid, [distance_col,times_col]].sort_values(
                                        by = times_col, ascending = True).iloc[0,0]
            
            # Fill time_points array
            time_points[idx] = df_particle.loc[df_particle[index_col] == fid, [distance_col,times_col]].sort_values(
                                        by = times_col, ascending = True).iloc[-1,1]
        
        
        # Create travel time distribution plot
        fig, ax = plt.subplots(figsize= (10,10),dpi=300)

        # Plot time-radius plot
        plt.plot(distance_points, time_points, lw = 0., marker = '.')  # analytical
        plt.plot(distance_points, time_points, lw = 0.2)  # modpath
        plt.xlabel("Radial distance (m)")
        plt.ylabel("Travel time (days)")

        # travel distribution plot (file name)
        if fpath_fig is None:
            fpath_fig = os.path.join(self.well.dstroot,self.well.schematisation_type + "_traveltime_distribution_modpath.png")
        # Save figure
        fig.savefig(fpath_fig)

    def plot_age_distribution(self, df_particle: pd.DataFrame,
                                vmin = 0.,vmax = 1.,orientation = {'row': 0},
                                fpathfig = None, figtext = None,x_text = 0,
                                y_text = 0, lognorm = True, xmin = 0., xmax = None,
                                line_dist = 1, dpi = 192, trackingdirection = "forward",
                                cmap = 'viridis_r'):
        ''' Create pathline plots with residence times 
            using colours as indicator.
            with: 
            - df_particle: dataframe containing xyzt-points of the particle paths.
            - fpathfig = output location of plots
            figtext: figure text to show within plot starting at
            x-location 'x_text' and y-location 'y_text'.
            lognorm: True/False (if vmin <= 0. --> vmin = 1.e-5)
            line_dist: min. distance between well pathlines at source (m)
            line_freq: show every 'X' pathlines
            dpi: pixel density
            trackingdirection: direction of calculating flow along pathlines"
            cmap: Uses colormap 'viridis_r' (viridis reversed as default)
            '''
   
        if lognorm:
            if vmin <= 0.:
                vmin = 1.e-2

        # Keep track of minimum line distance between pathlines ('line_dist')
        xmin_plot = 0.
        # Flowline IDs
        flowline_ID = list(df_particle.index.unique())
        for fid in flowline_ID:

            x_points = df_particle.loc[df_particle.index == fid, "xcoord"].values
#           y_points = df_particle.loc[df_particle.index == fid, "ycoord"].values
            z_points = df_particle.loc[df_particle.index == fid, "zcoord"].values
            # Cumulative travel times
            time_points = df_particle.loc[df_particle.index == fid, "total_travel_time"].values

            # Plot every 'line_dist' number of meters one line
            if trackingdirection == "forward":
                # x_origin: starting position of pathline
                x_origin = x_points[0]  
            else:
                # x_origin: starting position of pathline
                x_origin = x_points[-1]

            # Check if line should be plotted
            if x_origin > xmin_plot:
                xmin_plot += line_dist
                pass
            else:
                continue

            # Combine to xyz scatter array with x,y,values
            xyz_scatter = np.stack((x_points,z_points,time_points), axis = 1)

            # Plot function (lines))
            marker_size=0.5  # 's' in functie scatter
            if lognorm:
                # formatting of values (log or linear?)
                norm_vals = colors.LogNorm()
                # norm_labels = [str(iLog) for iLog in [8,7,6,5,4,3,2,1,0]]
            else:
                # formatting of values (log or linear: None?)
                norm_vals = None
                
            # Mask values outside of vmin & vmax
            time_vals = np.ma.masked_where((xyz_scatter[:,2] > vmax), xyz_scatter[:,2])
            plt.scatter(xyz_scatter[:,0],
                        xyz_scatter[:,1],
                        s = marker_size,
                        cmap = cmap,
                        c=time_vals,
                        marker = 'o',
                        norm= norm_vals)
            plt.plot(xyz_scatter[:,0],
                        xyz_scatter[:,1], c = 'k', lw = 0.1)  
                        # 'o', markersize = marker_size,
                        # markerfacecolor="None", markeredgecolor='black') #, lw = 0.1)
            plt.xlim(xmin,xmax)
            plt.clim(vmin,vmax)
        # Voeg kleurenbalk toe
        cbar = plt.colorbar()
        ticklabs_old = [t.get_text() for t in cbar.ax.get_yticklabels()]
        # print(ticklabs_old)
        # ticklabs_new = [str(iLab).replace("-0.0","0.0") for iLab in \
        #                 np.linspace(-np.log10(vmin),0.,num = len(ticklabs_old), endpoint = True)]
        # cbar.ax.set_yticklabels(ticklabs_new)
        
        # cbar.set_yticks([mn,md,mx])
        # cbar.ax.set_yticklabels(norm_labels)
        cbar.set_label("Residence time [days]")
    #    cbar.set_label("Concentratie [N/L]")
        # Titel
        # plt.title("Pathogenenverwijdering in grondwater")
        # Label x-as
        plt.xlabel("Distance (m)")
        # Label y-as
        plt.ylabel("Depth (m)")
        # Tekst voor in de figuur
        if figtext is not None:
            plt.text(x = x_text,y = y_text,s = figtext,
                    bbox={'facecolor': 'gray', 'alpha': 0.5, 'pad': 10})
            plt.subplots_adjust(left=0.5)
        if fpathfig is None:
            plt.show()
        else:
            plt.savefig(fpathfig, dpi = dpi)
        # Sluit figuren af
        plt.close('all')    


    def plot_logremoval(self, df_particle: pd.DataFrame,
                        df_flowline: pd.DataFrame,
                        vmin = 0.,vmax = 1.,orientation = {'row': 0},
                        fpathfig = None, figtext = None,x_text = 0,
                        y_text = 0, lognorm = True, xmin = 0., xmax = None,
                        line_dist = 1, dpi = 192, trackingdirection = "forward",
                        cmap = 'viridis_r'):
        ''' Create pathline plots with log removal
            using colours as indicator.
            with: 
            - df_particle: dataframe containing xyzt-points of the particle paths.
            - fpathfig = output location of plots
            figtext: figure text to show within plot starting at
            x-location 'x_text' and y-location 'y_text'.
            lognorm: True/False (if vmin <= 0. --> vmin = 1.e-5)
            line_dist: min. distance between well pathlines at source (m)
            line_freq: show every 'X' pathlines
            dpi: pixel density
            trackingdirection: direction of calculating flow along pathlines"
            cmap: Uses colormap 'viridis_r' (viridis reversed as default)
            '''
            
        if lognorm:
            if vmin <= 0.:
                vmin = 1.e-2

        # Keep track of minimum line distance between pathlines ('line_dist')
        xmin_plot = 0.
        # Flowline IDs
        flowline_ID = list(df_particle.index.unique())
        for fid in flowline_ID:

            x_points = df_particle.loc[fid, "xcoord"].values
#           y_points = df_particle.loc[fid, "ycoord"].values
            z_points = df_particle.loc[fid, "zcoord"].values
            # Concentration along pathlines
            conc_points = df_particle.loc[fid, "steady_state_concentration"].values
            # Starting concentration of pathline
            input_concentration = df_flowline.loc[fid,"input_concentration"]
            # Concentration relative to input
            relative_conc = conc_points / input_concentration

            # Plot every 'line_dist' number of meters one line
            if trackingdirection == "forward":
                # x_origin: starting position of pathline
                x_origin = x_points[0]  
            else:
                # x_origin: starting position of pathline
                x_origin = x_points[-1]

            # Check if line should be plotted
            if x_origin > xmin_plot:
                xmin_plot += line_dist
                pass
            else:
                continue

            # Combine to xyz scatter array with x,y,values
            xyz_scatter = np.stack((x_points,z_points,relative_conc), axis = 1)

            # Plot function (lines))
            marker_size=0.5  # 's' in functie scatter
            if lognorm:
                # formatting of values (log or linear?)
                norm_vals = colors.LogNorm()
                norm_labels = [str(iLog) for iLog in [20,18,16,14,12,10,8,6,4,2,0]]
            else:
                # formatting of values (log or linear: None?)
                norm_vals = None
                
            # Mask values outside of vmin & vmax
            conc_vals = np.ma.masked_where((xyz_scatter[:,2] < vmin) | (xyz_scatter[:,2] > vmax), xyz_scatter[:,2])
            plt.scatter(xyz_scatter[:,0],
                        xyz_scatter[:,1],
                        s = marker_size,
                        cmap = cmap,
                        c = conc_vals,
                        marker = 'o',
                        norm= norm_vals)
            plt.plot(xyz_scatter[:,0],
                        xyz_scatter[:,1], c = 'k', lw = 0.1)  
                        # 'o', markersize = marker_size,
                        # markerfacecolor="None", markeredgecolor='black') #, lw = 0.1)
            plt.xlim(xmin,xmax)
            plt.clim(vmin,vmax)
        # Voeg kleurenbalk toe
        cbar = plt.colorbar()
        ticklabs_old = [t.get_text() for t in cbar.ax.get_yticklabels()]
        # print(ticklabs_old)
        # ticklabs_new = [str(iLab).replace("-0.0","0.0") for iLab in \
        #                 np.linspace(-np.log10(vmin),0.,num = len(ticklabs_old), endpoint = True)]
        # cbar.ax.set_yticklabels(ticklabs_new)
        
        # cbar.set_yticks([mn,md,mx])
        cbar.ax.set_yticklabels(norm_labels)
        cbar.set_label("Log removal [-]")
    #    cbar.set_label("Concentratie [N/L]")
        # Titel
        # plt.title("Pathogenenverwijdering in grondwater")
        # Label x-as
        plt.xlabel("Distance (m)")
        # Label y-as
        plt.ylabel("Depth (m)")
        # Tekst voor in de figuur
        if figtext is not None:
            plt.text(x = x_text,y = y_text,s = figtext,
                    bbox={'facecolor': 'gray', 'alpha': 0.5, 'pad': 10})
            plt.subplots_adjust(left=0.5)
        if fpathfig is None:
            plt.show()
        else:
            plt.savefig(fpathfig, dpi = dpi)
        # Sluit figuren af
        plt.close('all')    
        
    def compute_pathogen_removal(self):
        #AH_todo
        pass

    # def plot_age_distribution(self):
    #     #AH_todo
    #     pass

    # def plot_logremoval(self):
    #     #AH_todo
    #     pass
