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
import datetime as dt
from datetime import timedelta

# module_path = os.path.abspath(os.path.join("..","..","sutra2
#"))
# if module_path not in sys.path:
#     sys.path.insert(0,module_path)

from sutra2.Analytical_Well import *
from sutra2.ModPath_Well import ModPathWell

# from Analytical_Well import AnalyticalWell
# from ModPath_functions import ModPathWell
# # from sutra2
#.ModPath_functions import ModPathWell
# # from Analytical_Well import AnalyticalWell

path = os.getcwd()  # path of working directory

class MicrobialOrganism:
    '''
    Placeholder class which will later be replaced by the QSAR functionality of AquaPriori ?!

    For removal of microbial organisms ('mbo') / pathogen

    Return
    -------
    microbial_species_dict: dict
    removal_function = 'mbo'
    '''

    def __init__(self, organism_name,
                alpha0_suboxic: int or float = None,
                alpha0_anoxic: int or float = None,
                alpha0_deeply_anoxic: int or float = None,
                pH0_suboxic: int or float = None,
                pH0_anoxic: int or float = None,
                pH0_deeply_anoxic: int or float = None,
                mu1_suboxic: int or float = None,
                mu1_anoxic: int or float = None,
                mu1_deeply_anoxic: int or float = None,
                organism_diam: int or float = None
                ):
        """
        Parameters
        ----------
        organism_name: str
            name of the organism (for now limited dictionary to
            'solani','carotovorum', 'solanacearum')
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

        """
        self.organism_name = organism_name
        # (default removal function)
        self.removal_function = 'mbo'

        # Dict    # Acties Steven 25-1-22

        # Naming convention organism: Uppercamelcase species
        self.micro_organism_database = {
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

        # Microbial Organism in db (default values)
        if self.organism_name in self.micro_organism_database.keys():
            self.organism_dict = self.micro_organism_database[self.organism_name]
        else: # return default values if microbial organism not in database
            self.organism_dict = \
                {"organism_name": self.organism_name,
                 "alpha0": {
                    "suboxic": 1.e-2,
                    "anoxic": 1.e-3,
                    "deeply_anoxic": 1.e-4
                    },
                 "pH0": {
                    "suboxic": 7.5,
                    "anoxic": 7.5,
                    "deeply_anoxic": 7.5
                    },
                 "organism_diam": 1.e-6,
                 "mu1": {
                    "suboxic": 1.,
                    "anoxic": 0.1,
                    "deeply_anoxic": 0.1
                     }
                }

        # add/overwrite parameters 'if values are not None'
        if alpha0_suboxic is not None:
            self.organism_dict["alpha0"]["suboxic"] = alpha0_suboxic
        if alpha0_anoxic is not None:
            self.organism_dict["alpha0"]["anoxic"] = alpha0_anoxic
        if alpha0_deeply_anoxic is not None:
            self.organism_dict["alpha0"]["deeply_anoxic"] = alpha0_deeply_anoxic
        if pH0_suboxic is not None:
            self.organism_dict["pH0"]["suboxic"] = pH0_suboxic
        if pH0_anoxic is not None:
            self.organism_dict["pH0"]["anoxic"] = pH0_anoxic
        if pH0_deeply_anoxic is not None:
            self.organism_dict["pH0"]["deeply_anoxic"] = pH0_deeply_anoxic
        if organism_diam is not None:
            self.organism_dict["organism_diam"] = organism_diam
        if mu1_suboxic is not None:
            self.organism_dict["mu1"]["suboxic"] = mu1_suboxic
        if mu1_anoxic is not None:
            self.organism_dict["mu1"]["anoxic"] = mu1_anoxic
        if mu1_deeply_anoxic is not None:
            self.organism_dict["mu1"]["deeply_anoxic"] = mu1_deeply_anoxic


class Substance:
    '''
    Placeholder class which will later be replaced by the QSAR functionality of AquaPriori.

    For removal of organic micro pollutants (omp)

    Return
    -------
    microbial_species_dict: dict
    removal_function = 'omp'
    '''

    def __init__(self, substance_name,
                molar_mass: int or float = None,
                partition_coefficient_water_organic_carbon: int or float = None,
                dissociation_constant: int or float = None,
                halflife_suboxic: int or float = None,
                halflife_anoxic: int or float = None,
                halflife_deeply_anoxic: int or float = None):
        '''
        Parameters
        -----------
        substance_name: String,
            substance_name of the substance (for now limited dictionary to 'benzene', 'AMPA', 'benzo(a)pyrene'
        molar_mass: float or int
            Molar mass of the substance [g/mol]
        partition_coefficient_water_organic_carbon: float
            Distribution coefficient of OMP between organic carbon and water, dimensionless.
        dissociation_constant: float
            Dissociation equilibirum constant of the OMP, dimensionless.
        halflife_suboxic, halflife_anoxic, halflife_deeply_anoxic: float
            Time required to reduce the concentration of the OMP by half, from any concentration point in time [days].

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

        '''

        # substance name
        self.substance_name = substance_name
        # (default removal function)
        self.removal_function = 'omp'

        # Substance dict here as placeholder for the actual database
        self.substance_database = {
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

        # Substance in db (default values)
        if self.substance_name in self.substance_database.keys():
            self.substance_dict = self.substance_database[self.substance_name]

        else: # return default values if substance not in database
            self.substance_dict = \
                {"substance_name": self.substance_name,
                "log_Koc": 0,
                "molar_mass": 10,
                'pKa': 99,
                "omp_half_life": {
                    "suboxic": 1e99,
                    "anoxic": 1e99,
                    "deeply_anoxic": 1e99
                    }
                }


        # add/overwrite parameters 'if values are not None'
        if partition_coefficient_water_organic_carbon is not None:
            self.substance_dict["log_Koc"] = partition_coefficient_water_organic_carbon
        if molar_mass is not None:
            self.substance_dict["molar_mass"] = molar_mass
        if dissociation_constant is not None:
            self.substance_dict["pKa"] = dissociation_constant
        if halflife_suboxic is not None:
            self.substance_dict["omp_half_life"]["suboxic"] = halflife_suboxic
        if halflife_anoxic is not None:
            self.substance_dict["omp_half_life"]["anoxic"] = halflife_anoxic
        if halflife_deeply_anoxic is not None:
            self.substance_dict["omp_half_life"]["deeply_anoxic"] = halflife_deeply_anoxic


#ah_todo @MartinK, MartinvdS -> let the user specify the chemical in the Substance transport file instead of schematisation?
# also let them feed it a dictionary with their own substance?


class Transport():
    """

    # @Steven --> verander in 'Transport'
    Substance + Organism --> Pollutant
    Returns concentration in a groundwater well for a given Organic Micro Pollutant or microbial species.

    Attributes
    ----------
    well: object
        The AnalyticalWell object for the schematisation of the aquifer type.
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

    pollutant: Class Substance or class MicrobialOrganism 

    Substance: object
        The Substance object with the OMP of interest.

    substance_dict: dictionary
        Nested dictionary with the following per substance.

        substance_name: String,
            substance_name of the substance (for now limited dictionary to 'benzene', 'AMPA', 'benzo(a)pyrene'
        log_Koc: float
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
                pollutant: Substance or MicrobialOrganism,
                start_date_well = dt.datetime.strptime('1950-01-01', "%Y-%m-%d"),
                start_date_contamination=None,
                end_date_contamination= None,
                compute_contamination_for_date=None,

                ):


        '''
        removal_parameters:

        Initialization of the Transport class. Checks for either user-defined OMP substance parameters (removal_function == 'omp')
        or user-defined microbial organism removal parameters (removal_function == 'mbo') and overrides the database values, if any.

        Parameters
        ----------
        well: object
            The AnalyticalWell object for the schematisation of the aquifer type.
        pollutant: object
            Substance: object
                The Substance object with the OMP of interest.
            MicrobialOrganism: object
                The Organism object with microbial organism (MBO) of interest
        start_date_well: dt.datetime.strptime('YYYY-MM-DD', "%Y-%m-%d"),
            Start date of the well.
        start_date_contamination, end_date_contamination: dt.datetime.strptime('YYYY-MM-DD', "%Y-%m-%d"),
            Start and end date for both the diffuse and point sources.
            @MartinvdS, is this a problem to use the start/end date for both?
        compute_contamination_for_date: dt.datetime.strptime('YYYY-MM-DD', "%Y-%m-%d"),
            Date for which to compute the contamination in the well.

        removal_function (inherited from Substance ('omp') or MicrobialOrganism ('mbo')): str
            removal_function: ['omp' or 'mbo']
            Calculate removal of either organic micro pollutants ('omp') or microbial organisms ('mbo')
        '''

        # well: AnalyticalWell or Modpath_Well object
        self.well = well
        # Load input dataframes
        self.df_particle = well.df_particle
        self.df_flowline = well.df_flowline

        # pollutant properties (Substance or MicrobialOrganism)
        self.pollutant = pollutant
        # removal_function ('omp' or 'mbo')
        self.removal_function = self.pollutant.removal_function

        # Run init of 'omp'
        if self.removal_function == 'omp':
            self.omp_initialized = False
            self.micro_organism_initialized = True

        # Run init of 'pathogen'
        elif self.removal_function == 'mbo':
            self.omp_initialized = True
            self.micro_organism_initialized = False

        if type(self.pollutant) == Substance:
            # Substance
            self.pollutant_name = self.pollutant.substance_name
            self.removal_parameters = self.pollutant.substance_dict


        elif type(self.pollutant) == MicrobialOrganism:
            # Microbial Organism
            self.pollutant_name = self.pollutant.organism_name
            self.removal_parameters = self.pollutant.organism_dict

        # Inilitialize date information
        def check_date_format(check_date):
            for var in check_date:
                value = getattr(self, var)
                if not type(value) is dt.datetime:
                    raise TypeError(f"Error invalid date input, please enter a new {var} using the format dt.datetime.strptime('YYYY-MM-DD', '%Y-%m-%d')")

        self.start_date_well = start_date_well

        check_date_format(check_date=['start_date_well'])

        # Contamination
        if start_date_contamination is None:
            self.start_date_contamination = self.start_date_well
        else:
            self.start_date_contamination = start_date_contamination
            check_date_format(check_date=['start_date_contamination'])

        if compute_contamination_for_date is None:
            self.compute_contamination_for_date = self.start_date_well + timedelta(days=365.24*50)
        else:
            self.compute_contamination_for_date = compute_contamination_for_date
            check_date_format(check_date=['compute_contamination_for_date'])

        if end_date_contamination is None:
            self.end_date_contamination = end_date_contamination
        else:
            self.end_date_contamination = end_date_contamination
            check_date_format(check_date=['end_date_contamination'])

            if self.end_date_contamination < self.start_date_contamination:
                raise ValueError('Error, "end_date_contamination" is before "start_date_contamination". Please enter an new "end_date_contamination" or "start_date_contamination" ')


        ''' Check logical things here for contamination'''
        if self.compute_contamination_for_date < self.start_date_contamination:
            raise ValueError('Error, "compute_contamination_for_date" is before "start_date_contamination". Please enter an new "compute_contamination_for_date" or "start_date_contamination" ')
        if self.compute_contamination_for_date < self.start_date_well:
            raise ValueError('Error, "compute_contamination_for_date" is before "start_date_well". Please enter an new "compute_contamination_for_date" or "start_date_well" ')
        #AH_todo @MartinvdS -> if end_date_contamination < start_date_well what to do?


        # Relative dates of well abstraction and pollutant contamination date
        self._contamination_date_vs_well_abstraction()

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
            # Add "organism_diam" to df_flowline (is independent of node)
            self.df_flowline['organism_diam'] = self.removal_parameters['organism_diam']

        self.micro_organism_initialized = True

    #@Steven_todo, @MartinK: include in a function (used for both omp_removal as well as mbo removal)
    # @MartinK / MartinvdS: add as attributes to 'Transport' class:
    #       compute_contamination_for_date,start_date_well,start_date_contamination,end_date_contamination

    def _contamination_date_vs_well_abstraction(self):
        # reduce the amount of text per line by extracting the following parameters
        self.compute_contamination_for_date = self.compute_contamination_for_date
        start_date_well = self.start_date_well
        start_date_contamination = self.start_date_contamination
        self.end_date_contamination = self.end_date_contamination

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
        if self.well.biodegradation_sorbed_phase:
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

        if self.well.temp_correction_halflife:
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
        if self.well.temp_correction_Koc:
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
        DOC_inf = self.well.dissolved_organic_carbon_infiltration_water
        TOC_inf = self.well.total_organic_carbon_infiltration_water

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
            self.df_flowline.at[fid, 'total_breakthrough_travel_time'] = sum(self.df_particle.loc[fid,:].fillna(0)['travel_time'])
            self.df_flowline.at[fid, 'breakthrough_concentration'] = self.df_particle.loc[fid,'steady_state_concentration'].iloc[-1]

    # this will introduce an error if user calls it while he has a microb class. so I would make a get_removal method that chekcs whether it is a
    # omp or mbo simulation and call their resepecive removal method; these methods should be private (prefix _)
    def compute_omp_removal(self):
        #MWK: the description of df_flowline and df_particle is not the best place here. Better in __init__. only say here
        # what you will add
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
        self.df_flowline.loc[:,'input_concentration'] = self.well.diffuse_input_concentration
        self.df_particle.loc[:,'input_concentration'] = None
        if 'steady_state_concentration' not in self.df_particle.columns:
            # Steady state concentration [-]
            self.df_particle.loc[:,'steady_state_concentration'] = None

        # First value of input_concentration/steady_state_concentration equals diffuse input concentration of pathline
        self.df_particle.loc[self.df_particle.total_travel_time == 0.,"input_concentration"] = self.df_flowline.loc[:,'input_concentration'].values
        self.df_particle.loc[self.df_particle.total_travel_time == 0.,"steady_state_concentration"] = self.df_flowline.loc[:,'input_concentration'].values


        if self.well.point_input_concentration:
            ''' point contamination '''
            # need to take into account the depth of the point contamination here....
            # need to change the df_particle and df_flowline to only be the flowlines for the point contamination flowline(s)
            # use a single point contamination for now
            # FIRST recalculate the travel times for the contamination, then initialize the class

            #only for a SINGLE point contamination
            distance = self.well.distance_point_contamination_from_well
            depth = self.well.depth_point_contamination
            cumulative_fraction_abstracted_water = (math.pi * self.well.recharge_rate
                                                        * distance ** 2)/abs(self.well.well_discharge)
            ind = self.df_particle.flowline_id.iat[-1]

            if self.well.schematisation_type == 'phreatic':
                head = calculate_hydraulic_head_phreatic(distance = distance,
                                             groundwater_level = self.well.groundwater_level,
                                             well_discharge = self.well.well_discharge,
                                             KD = self.well.KD,
                                             radial_distance_recharge = self.well.radial_distance_recharge)

                df_flowline_points, df_particle_points = self.well._add_phreatic_point_sources(distance=distance,
                                            depth_point_contamination=depth,
                                            cumulative_fraction_abstracted_water=cumulative_fraction_abstracted_water)

            elif self.well.schematisation_type == 'semiconfined':
                bottom_vadose_zone = self.well.bottom_vadose_zone_at_boundary

                df_flowline_points, df_particle_points = self.well._add_semiconfined_point_sources(distance=distance,
                                        depth_point_contamination=depth,  )
            # Assign flowline id to point
            df_particle_points.loc[:,'flowline_id'] = df_particle_points.loc[:,'flowline_id'].values + ind
            df_flowline_points.loc[:,'flowline_id'] = df_flowline_points.loc[:,'flowline_id'].values + ind
            # reindex df_particle and df_flowline
            df_particle_points.index = df_particle_points.flowline_id
            df_flowline_points.index = df_flowline_points.flowline_id

            # Add (input) concentration to point source df's
            df_flowline_points.loc[:,'input_concentration'] = self.well.point_input_concentration
            df_particle_points.loc[:,'input_concentration'] = None
            df_particle_points.loc[:,'steady_state_concentration'] = None

            # First value of input_concentration/steady_state_concentration equals point input concentration of pathline
            # for fid in df_particle_points.flowline_id.unique():

            df_particle_points.loc[df_particle_points.total_travel_time == 0., 'input_concentration'] = df_flowline_points.loc[:,'input_concentration'].values
            df_particle_points.loc[df_particle_points.total_travel_time == 0., 'steady_state_concentration'] = df_flowline_points.loc[:,'input_concentration'].values

            # Add flow/discharge data to point particle df's
            df_flowline_points.loc[:,'flowline_type'] = "point_source"
            df_flowline_points.loc[:,'flowline_discharge'] = abs(self.well.discharge_point_contamination)

            #AH_todo, something here to loop through different point sources if more than one.

            # Add point particles to self.df_particle and self.df_flowline & reindex based on flowline_id
            self.df_particle = self.df_particle.append(df_particle_points)
            self.df_particle.index = self.df_particle.loc[:,"flowline_id"].values
            # self.df_particle.reset_index(drop=True, inplace=True)

            self.df_flowline = self.df_flowline.append(df_flowline_points)
            self.df_flowline.index = self.df_flowline.loc[:,"flowline_id"].values
            # self.df_flowline.reset_index(drop=True, inplace=True)

            # Add substance name
            self.df_flowline.loc[:,'name'] = self.pollutant_name

        # MWK you call it twice now... also in __inti__
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
        # --> MWK: Probably in the for loop over time_array.few things to try:
        #          - initialize well_concentration to the length it will get: well_concentration=[0]*len(time_array)
        #          - do the summation only once before the loop by adding a column cumulative_well_concentration = df_flowline['concentration_in_well'].cumsum()
        #            and get the right well_conc from that column instead. I think you dont even need to loop over t anymore
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

    # I would add a show_fig = True (as default) to show the fig directly with this method
    def plot_concentration(self,
                            xlim=None,
                            ylim=None,
                            as_fraction_input = None,
                            x_axis = 'Date'):
        ## @Alex: adjust name of function to 'plot_concentration_timeseries'
        # MWK why is this not done yet?
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
        point_input_concentration = self.well.point_input_concentration
        diffuse_input_concentration = self.well.diffuse_input_concentration
        schematisation_type = self.well.schematisation_type
        compute_contamination_for_date = self.compute_contamination_for_date
        start_date_well = self.well.start_date_well
        start_date_contamination = self.well.start_date_contamination
        end_date_contamination = self.well.end_date_contamination

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
            df_well_concentration.loc[:,"total_concentration_in_well"] = df_well_concentration.loc[:,"total_concentration_in_well"].values / input_concentration
            ylabel = 'Fraction of input concentration'
        else:
            ylabel = 'Concentration (ug/L)'

        fig = plt.figure(figsize=[10, 5])
        if x_axis == 'Date':
            plt.plot(df_well_concentration.date, df_well_concentration.total_concentration_in_well, 'b', label =str(self.pollutant_name))
            plt.axvline(x=start_date_well, color= 'k', label = 'Start date well')
            plt.axvline(x=start_date_contamination, color= 'r', label = 'Start date contamination')
            if end_date_contamination is not None:
                plt.axvline(x=end_date_contamination, color= 'g', label = 'End date contamination')

            plt.xlabel('Date')
            if xlim is None:
                plt.xlim([dt.date((back_date_start.year-5), 1, 1), compute_contamination_for_date])
            else: plt.xlim(xlim)

        elif x_axis == 'Time':
            plt.plot(df_well_concentration.time/365.24, df_well_concentration.total_concentration_in_well,  'b', label =str(self.pollutant_name))
            plt.axvline(x=(start_date_well-start_date).days/365.24, color= 'k', label = 'Start date well')
            plt.axvline(x=(start_date_contamination-start_date).days/365.24, color= 'r', label = 'Start date contamination')
            if end_date_contamination is not None:
                plt.axvline(x=(end_date_contamination-start_date).days/365.24, color= 'g', label = 'End date contamination')

            plt.xlabel('Time (years)')
            if xlim is None:
                plt.xlim([(back_date_start-start_date).days/365-5,compute_date.days/365.24])
            else: plt.xlim(xlim)

        if ylim is None:
            plt.ylim([0, input_concentration])
        else: plt.ylim(ylim)
        plt.legend(loc=2)
        plt.ylabel(ylabel)
        plt.title('Aquifer type: ' + schematisation_type)
        plt.grid()
        plt.legend()
        # plt.savefig('well_concentration_over_time_'+str(self.pollutant_name)+'_'+schematisation_type+'.png', dpi=300, bbox_inches='tight')
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

        if df_column not in df.columns:
            # Add dataframe column with default value
            df[df_column] = value
        elif df[df_column].dropna().empty:
            df[df_column] = df[df_column].fillna(value)
        elif dtype_ == 'object':
            df[df_column] = df[df_column].fillna(value)
        else:
            value_mean = df[df_column].values.mean()
            df[df_column] = df[df_column].fillna(value_mean)
        # Adjust dtype
        df = df.astype({df_column: dtype_})

        return df

    # MWK: where is OMP and where is mbo removal determined?

    # this can be made a function, just likei others. is there only one df_particle df_flowline  associated with this class?
    # if so store them as ttribure and remove the as argumetns
    # this method is too long split it up
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
            tdiff[pid][tdiff[pid] == 0.] = 0.001
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

    # docstring too long, method to long, refactor in multiple methods
    def calc_advective_microbial_removal(self,df_particle: pd.DataFrame =None,
                                        df_flowline: pd.DataFrame =None,
                                        endpoint_id: str = "well1", trackingdirection = "forward",
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
            organism_name = self.pollutant_name

        if df_particle is not None:
            self.df_particle = df_particle
        if df_flowline is not None:
            self.df_flowline = df_flowline

        # Organism/Species name
        self.df_flowline['name'] = organism_name

        # Initialize mbo-parameters
        self._init_micro_organism()

        # Calculate removal coefficient 'lamda' [/day].
        self.df_particle, self.df_flowline = self.calc_lambda(self.df_particle, self.df_flowline,
                                                    redox = redox, mu1 = mu1, por_eff = por_eff,
                                                    grainsize = grainsize, alpha0 = alpha0,
                                                    pH = pH, pH0 = pH0, temp_water = temp_water,
                                                    rho_water = rho_water, organism_diam = organism_diam)

        # Starting concentration [-]
        self.df_flowline = self._df_fillna(self.df_flowline,
                                df_column = 'input_concentration',
                                value = conc_start,
                                dtype_ = 'float')

        # Original concentration in groundwater [-]
        self.df_flowline = self._df_fillna(self.df_flowline,
                                df_column = 'input_concentration_gw',
                                value = conc_gw,
                                dtype_ = 'float')

        if 'breakthrough_concentration' not in self.df_flowline.columns:
            # breakthrough concentration for particles along a certain flowline
            self.df_flowline['breakthrough_concentration'] = None

        if 'breakthrough_traveltime' not in self.df_flowline.columns:
            # breakthrough traveltime for particles along a certain flowline
            self.df_flowline['breakthrough_travel_time'] = None

        if 'concentration_in_well' not in self.df_particle.columns:
            # Final (averaged) steady state concentration in well at endpoint_id [N/L]
            self.df_particle['concentration_in_well'] = None

        if 'steady_state_concentration' not in self.df_particle.columns:
            # steady-state concentration for particles along a certain flowline
            self.df_particle['steady_state_concentration'] = None

        for pid in self.df_flowline.index:
            self.df_particle.loc[pid,"steady_state_concentration"].iloc[0] = self.df_flowline.at[pid,"input_concentration"]

        # Start dictionaries for relative and steady-state concentrations
        conc_rel = {}
        conc_steady = {}  # concentration corrected for initial concentration
        for pid in self.df_flowline.index:

            # Calculate the relative removal along pathlines
            exp_arg = -((self.df_particle.loc[pid,"lamda"].values / self.df_particle.loc[pid,"porewater_velocity"].values) *
                                self.df_particle.loc[pid,'relative_distance'].values).astype('float')

            conc_rel[pid] = self.df_flowline.loc[pid,"input_concentration_gw"] + \
                            (self.df_flowline.loc[pid,"input_concentration"] - self.df_flowline.loc[pid,"input_concentration_gw"]) * \
                                np.exp(exp_arg)

            if trackingdirection == 'forward':
                conc_steady[pid] = (conc_rel[pid]/self.df_flowline.loc[pid,"input_concentration"]).cumprod() * \
                                    self.df_flowline.loc[pid,"input_concentration"]
                # Replace 'nan'-values by 0.
                conc_steady[pid]  = np.nan_to_num(conc_steady[pid],0.)

                # Index of average final concentration "C_"  idx=-1[C0,...,..,...,Ct-1,C_final]
                idx_final = -1

            elif trackingdirection == 'backward':
                conc_steady[pid] = (conc_rel[pid][::-1]/self.df_flowline.loc[pid,"input_concentration"]).cumprod() * \
                                    self.df_flowline.loc[pid,"input_concentration"]
                # Replace 'nan'-values by 0.
                conc_steady[pid]  = np.nan_to_num(conc_steady[pid],0.)
                # First value represents concentration of endpoint and v.v.
                # Reverse array--> [C_final,Ct-1,..,...,C0]
                conc_steady[pid] = conc_steady[pid][::-1]
                # Index of average final concentration "C_" idx=0 [C_final,Ct-1,..,...,C0]
                idx_final = 0

            # Add steady_state_concentration to df_particle
            self.df_particle.loc[pid,'steady_state_concentration'] = conc_steady[pid]

        # Calculate average final concentration
        C_final = 0.
        for pid in self.df_flowline.index:
            if self.df_flowline.loc[pid,"endpoint_id"] == endpoint_id:
                # Steady-state concentration after breakthrough
                C_breakthrough = conc_steady[pid][idx_final]
                # Traveltime from contamination to endpoint location [days]
                time_breakthrough = self.df_particle.loc[pid,"total_travel_time"].iloc[-1]

                # Steady-state final concentration
                C_final += (C_breakthrough * self.df_flowline.loc[pid,"flowline_discharge"]) / \
                            self.df_flowline.loc[pid,"well_discharge"]

                # Add breakthrough_concentration to df_flowline
                self.df_flowline.loc[pid,'breakthrough_concentration'] = C_breakthrough
                # Add breakthrough_traveltime to df_flowline
                self.df_flowline.loc[pid,'breakthrough_travel_time'] = time_breakthrough

        # Add final concentration in well (at endpoint_id)
        self.df_flowline.loc[self.df_flowline.endpoint_id == endpoint_id,"concentration_in_well"] = C_final

        # No retardation included with mbo removal
        self.df_particle['retardation'] = 1.
        self.df_particle.loc[:,'breakthrough_travel_time'] = self.df_particle.loc[:,"retardation"].values * self.df_particle.loc[:,"total_travel_time"].values

        # SR_todo: @MartinvdS required to calc total_breakthrough_travel_time mbo?
        self._calculate_total_breakthrough_travel_time()

        self._contamination_date_vs_well_abstraction()

        # add the particle release date
        self.df_flowline['particle_release_day'] = self.particle_release_day # (self.start_date - start_date_contamination).days

        # return (adjusted) df_particle and df_flowline
        return self.df_particle, self.df_flowline, C_final

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
                                ymin = None, ymax = None,
                                line_dist = 1, dpi = 192, trackingdirection = "forward",
                                cmap = 'viridis_r',
                                show_vadose = True):
        ''' Create pathline plots with residence times
            using colours as indicator.
            with:
            - df_particle: dataframe containing xyzt-points of the particle paths.
            - fpathfig = output location of plots
            figtext: figure text to show within plot starting at
            x-location 'x_text' and y-location 'y_text'.
            lognorm: True/False (if vmin <= 0. --> vmin = 1.e-5)
            xmin: left boundary figure
            xmax: right boundary figure
            ymin: lower boundary figure
            ymax: upper boundary figure
            line_dist: min. distance between well pathlines at source (m)
            line_freq: show every 'X' pathlines
            dpi: pixel density
            trackingdirection: direction of calculating flow along pathlines"
            cmap: Uses colormap 'viridis_r' (viridis reversed as default)
            show_vadose: T/F (if True, include vadose zone flow lines; default = True)
            '''

        if lognorm:
            if vmin <= 0.:
                vmin = 1.e-2

        # Keep track of minimum line distance between pathlines ('line_dist')
        xmin_plot = 0.
        # Flowline IDs
        try:
            flowline_ID = list(df_particle.flowline_id.unique())
        except:
            # Use index as flowline_id instead
            df_particle.loc[:,"flowline_id"] = df_particle.index.values
            flowline_ID = list(df_particle.flowline_id.unique())


        for fid in flowline_ID:

            if show_vadose:
                x_points = df_particle.loc[df_particle.flowline_id == fid, "xcoord"].values
    #           y_points = df_particle.loc[df_particle.flowline_id == fid, "ycoord"].values
                z_points = df_particle.loc[df_particle.flowline_id == fid, "zcoord"].values
                # Cumulative travel times
                time_points = df_particle.loc[df_particle.flowline_id == fid, "total_travel_time"].values
            else:
                x_points = df_particle.loc[(df_particle.flowline_id == fid) & (~df_particle.zone.str.contains("vadose_zone")), "xcoord"].values
    #           y_points = df_particle.loc[(df_particle.flowline_id == fid) & (~df_particle.zone.str.contains("vadose_zone")), "ycoord"].values
                z_points = df_particle.loc[(df_particle.flowline_id == fid) & (~df_particle.zone.str.contains("vadose_zone")), "zcoord"].values
                # Cumulative travel times
                time_points = df_particle.loc[(df_particle.flowline_id == fid) & (~df_particle.zone.str.contains("vadose_zone")), "total_travel_time"].values

            # Plot every 'line_dist' number of meters one line
            if trackingdirection == "forward":
                # x_origin: starting position of pathline
                x_origin = x_points[0]
                if not show_vadose:
                    # Reduce time by initial time value out of vadose_zone
                    time_points -= time_points[0]
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

            # ymin, ymax
            plt.ylim(ymin,ymax)
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
        cbar.ax.set_yticklabels(norm_labels)
        cbar.set_label("Log removal [-]")

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
        raise NotImplementedError()

