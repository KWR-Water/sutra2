#%% ----------------------------------------------------------------------------
# INITIALISATION OF PYTHON e.g. packages, etc.
# ------------------------------------------------------------------------------

import numpy as np
import pandas as pd
import os
import sys
from pandas import read_csv
from pandas import read_excel
import math
import datetime
from datetime import timedelta

path = os.getcwd()



class Organism:
    ''' 
    Placeholder class which will later be replaced by the QSAR functionality of AquaPriori ?!

    For removal of microbial organisms ('mbo') / pathogen 
    
    microbial_species_dict: dict
    Attributes
    ---------
    organism_name: String
        species_name of the substance (for now limited dictionary to 'MS2' (MS2-virus)
        
    'alpha0': float
        reference_collision_efficiency [-]
        per redox zone ('suboxic', 'anoxic', deeply_anoxic')
    'reference_pH': float
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
        organism: String,
            name of the substance (for now limited dictionary to 'benzene', 'AMPA', 'benzo(a)pyrene'

        Returns
        -------
        organism_dict: dictionary
            'alpha0': float
                reference_collision_efficiency [-]
                per redox zone ('suboxic', 'anoxic', deeply_anoxic')
            'reference_pH': float
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
            "MS2": 
                {"organism_name": "MS2",
                    "alpha0": {
                        "suboxic": 1.e-3, 
                        "anoxic": 1.e-5, 
                        "deeply_anoxic": 1.e-5
                    },
                    "reference_pH": {
                        "suboxic": 6.6, 
                        "anoxic": 6.8, 
                        "deeply_anoxic": 6.8
                    },
                    "organism_diam": 2.33e-8,
                    "mu1": {
                        "suboxic": 0.039, 
                        "anoxic": 0.023, 
                        "deeply_anoxic": 0.023
                    }
                },
                # ALLES X 10
            "Escherichia coli": 
                {"organism_name": "Escherichia coli",
                    "alpha0": {
                        "suboxic": 1.e-2, 
                        "anoxic": 1.e-4, 
                        "deeply_anoxic": 1.e-4
                    },
                    "reference_pH": {
                        "suboxic": 6.6, 
                        "anoxic": 6.8, 
                        "deeply_anoxic": 6.8
                    },
                    "organism_diam": 1.5e-6,
                    "mu1": {
                        "suboxic": 0.39, 
                        "anoxic": 0.23, 
                        "deeply_anoxic": 0.23
                    }
                }
            }
        #@ Steven voeg toe: micro_organism_dict
        # self.micro_organism_dict = micro_organism_dict[substance_name]
        self.organism_dict = micro_organism_dict[self.organism_name]


class MicrobialRemoval():
    '''
    Returns concentration for a given microbial species.


    organism: object
        The microbial organism object with the microbial organism (mbo) of interest

        organism_dict: dictionary
            'alpha0': float
                reference_collision_efficiency [-]
                per redox zone ('suboxic', 'anoxic', deeply_anoxic')
            'reference_pH': float
                reference pH for calculating collision efficiency [-]
                per redox zone ('suboxic', 'anoxic', deeply_anoxic')
            'organism_diam': float
                diameter of pathogen/species [m]
            'mu1': float
                inactivation coefficient [1/day]
                per redox zone ('suboxic', 'anoxic', deeply_anoxic')
    
    '''

    def __init__(self,
                organism: Organism = 'MS2',
                partition_coefficient_water_organic_carbon=None,
                dissociation_constant=None,
                halflife_suboxic=None,
                halflife_anoxic=None,
                halflife_deeply_anoxic=None,
                alpha0_suboxic=None,
                alpha0_anoxic=None,
                alpha0_deeply_anoxic=None,
                reference_pH_suboxic=None,
                reference_pH_anoxic=None,
                reference_pH_deeply_anoxic=None,
                mu1_suboxic=None,
                mu1_anoxic=None,
                mu1_deeply_anoxic=None,
                organism_diam=None,
                ):

        '''
        Initialization of the MicrobialRemoval class, checks for user-defined 
        microbial organism removal parameters and overrides the database values.
        '''

        # Organism
        self.organism_name = organism
        self.alpha0_suboxic=alpha0_suboxic
        self.alpha0_anoxic=alpha0_anoxic
        self.alpha0_deeply_anoxic=alpha0_deeply_anoxic
        self.reference_pH_suboxic=reference_pH_suboxic
        self.reference_pH_anoxic=reference_pH_anoxic
        self.reference_pH_deeply_anoxic=reference_pH_deeply_anoxic
        self.mu1_suboxic=mu1_suboxic
        self.mu1_anoxic=mu1_anoxic
        self.mu1_deeply_anoxic=mu1_deeply_anoxic
        self.organism_diam=organism_diam        

        # Load (default) microbial organism data
        self.organism = Organism(organism_name = organism)

        # Create user dict with 'removal_parameters' from input
        # if self.removal_function == 'omp':
        user_removal_parameters = {
            self.organism_name:
                {"organism_name": self.organism_name,
                    "alpha0": {
                        "suboxic": self.alpha0_suboxic, 
                        "anoxic": self.alpha0_anoxic, 
                        "deeply_anoxic": self.alpha0_deeply_anoxic
                    },
                    "reference_pH": {
                        "suboxic": self.reference_pH_suboxic, 
                        "anoxic": self.reference_pH_anoxic, 
                        "deeply_anoxic": self.reference_pH_deeply_anoxic
                    },
                    "organism_diam": self.organism_diam,
                    "mu1": {
                        "suboxic": self.mu1_suboxic, 
                        "anoxic": self.mu1_anoxic, 
                        "deeply_anoxic": self.mu1_deeply_anoxic
                    }
                },
            }

        
        # User defined removal parameters [omp]
        user_removal_parameters = user_removal_parameters[self.organism_name]
        # assumes that default dict contains microbial organism input (only MS2 currently supported)
        default_removal_parameters = self.organism.organism_dict

        # iterate through the dictionary keys
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
        
    def calc_lambda(self, redox = 'anoxic',
                mu1 = 0.149, mu1_std = 0.0932,
                por_eff = 0.33,
                grainsize = 0.00025,
                pH = 7.5,
                temp_water = 10.,
                rho_water = 999.703,
                alpha0 = 0.001,
                reference_pH = 7.5,
                organism_diam = 2.33e-8, v_por = 0.01):
        
        ''' For more information about the advective microbial removal calculation: 
                BTO2012.015: Ch 6.7 (page 71-74)

        Calculate removal coefficient 'lambda_' [/day].

            lambda_ = k_att + mu_1
            with 'hechtingssnelheidscoëfficiënt k_att [/day] and
            inactivation constant 'mu1' [/day] - sub(oxic): 0.149 [/day]
            lognormal standarddev. 'mu1_std' - sub(oxic): 0.0932 [/day]

        First, calculate "hechtingssnelheidscoëfficiënt" 'k_att' [/day]
            for Particle Paths with id_vals 'part_idx'
        # Effective porosity ('por_eff') [-]
        # grain size 'grainsize' [m]
        # collision efficiency ('botsingsefficiëntie') 'bots_eff' [-]
        # Boltzmann constant (const_BM) [1,38 × 10-23 J K-1] 
        # Water temperature 'temp_water' [degrees celcius]
        # Water density 'rho_water' [kg m-3]
        # Organism/species diameter 'organism_diam' [m]
        # porewater velocity 'v_por' [m/d]
        
        # Check - (example 'E_coli'):
        >> k_att = calc_katt(part_idx = [0], por_eff = [0.33], korrelgrootte = [0.00025],
                bots_eff = 0.001, const_BM = 1.38e-23,
                temp_water = [10.], rho_water = [999.703],
                organism_diam = 2.33e-8, v_por = [0.01])
        >> print(k_att)
        >> {0: 0.7993188853572424} # [/day]
        '''

        # Boltzmann coefficient [J K-1]
        const_BM = 1.38e-23    

        # Collision efficiency
        coll_eff = alpha0 * 0.9**((pH - reference_pH)/0.1)

        # Collision term 'k_coll'
        k_coll = (3/2.)*((1-por_eff) / grainsize) * coll_eff

        # Porosity dependent variable 'gamma'
        gamma = (1-por_eff)**(1/3)

        # Calculate Happel’s porosity dependent parameter 'A_s' (Eq. 5: BTO2012.015)
        ''' !!! Use correct formula:-> As =  2 * (1-gamma**5) /  (2 - 3 * gamma + 3 * gamma**5 - 2 * gamma**6)
            instead of... 2 * (1-gamma)**5 / (.......) 
        '''
        As_happ = 2 * (1-gamma**5) / \
                (2 - 3 * gamma + 3 * gamma**5 - 2 * gamma**6)

        # Dynamic viscosity (mu) [kg m-1 s-1]
        mu = (rho_water * 497.e-6) / \
                    (temp_water + 42.5)**(3/2)


        # Diffusion constant 'D_BM' (Eq.6: BTO2012.015) --> unit: [m2 s-1]
        D_BM = (const_BM * (temp_water + 273.)) / \
                    (3 * np.pi * organism_diam * mu)
        # Diffusieconstante 'D_BM' (Eq.6: BTO2012.015) --> unit: [m2 d-1]
        D_BM *= 86400.

        # Diffusion related attachment term 'k_diff'
        k_diff = ((D_BM /
                    (grainsize * por_eff * v_por))**(2/3) * v_por)

        # hechtingssnelheidscoëfficiënt k_att [/dag] :'attachment coefficient related to porewater velocity'
        k_att = k_coll * 4 * As_happ**(1/3) * k_diff
        # removal coefficient 'lambda_' [/day], using the 'mu1' mean.
        lambda_ = k_att + mu1

        return lambda_, k_att
    
    def calc_advective_microbial_removal(self, grainsize = 0.00025,
                                        temp_water = 11., rho_water = 999.703,
                                        pH = 7.5, por_eff = 0.33,
                                        conc_start = 1., conc_gw = 0.,
                                        redox = 'anoxic',
                                        distance_traveled = 1., traveltime = 100.,
                                        mu1 = None, alpha0 = None, reference_pH = None,
                                        organism_diam = None):
        ''' Calculate the advective microbial removal of microbial organisms
            from source to end_point.

            For more information about the advective microbial removal calculation: 
                BTO2012.015: Ch 6.7 (page 71-74)

            Calculate the steady state concentration along traveled distance per 
            node for each pathline from startpoint to endpoint_id'.
            
            With removal coefficient 'lambda_' [/d] (redox dependent)
            effective porosity 'porosity' [-]
            Starting concentration 'conc_start' per pathline
            Initial groundwater concentration 'conc_gw'
            Distance between points 'distance_traveled' [m]
            Time between start and endpoint 'traveltime' [days]

            
        '''

        # mu1 [day -1]
        if mu1 is None:
            mu1 = self.removal_parameters['mu1'][redox]

        # alpha0 [-]
        if alpha0 is None:
            alpha0 = self.removal_parameters['alpha0'][redox]

        # reference pH [-]
        if reference_pH is None:
            reference_pH = self.removal_parameters['reference_pH'][redox]

        # organism diameter [m]
        if organism_diam is None:
            organism_diam = self.removal_parameters['organism_diam']

        # porewater_velocity
        v_por = distance_traveled / traveltime

        # Calculate removal coefficient lambda [day -1]
        self.lambda_, self.k_att = self.calc_lambda(redox = redox, mu1 = mu1,
                                    por_eff = por_eff, grainsize = grainsize, 
                                    pH = pH, 
                                    temp_water = temp_water, 
                                    rho_water = rho_water,
                                    alpha0 = alpha0, 
                                    reference_pH = reference_pH,
                                    organism_diam = organism_diam)

        # Calculate concentration after microbial removal in subsurface
        C_final = (conc_start - conc_gw) * np.exp(-(self.lambda_/v_por)*distance_traveled) + conc_gw


        # return final concentration 'C_final'
        return C_final

