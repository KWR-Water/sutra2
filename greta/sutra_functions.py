

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
