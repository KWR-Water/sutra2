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

class Substance:
    def __init__(self, substance_name, ):
        """
        substance_name: String, 
            substance_name of the substance (for now limited dictionary to 'benzene', 'AMPA', 'benzo(a)pyrene'
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

        self.substance_dict = substances_dict[substance_name]

class SubstanceTransport():
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
        self.schematisation = schematisation
        self.omp_inialized = False
        self.df_particle = schematisation.df_particle
        self.df_flowline = schematisation.df_flowline
        self.substance = Substance(substance) 


        # @MartinK - need to make sure here that the substance passed is the same, e.g. comapre the dictionaries BUT ALSO
        # make sure that user doesn't call one substance in the hydrochemicalschematisation class and another in the concentration class
        # probably only a problem for ourselves, this should be written into a larger "run" class for the model which could avoid this
        if self.substance.substance_name == self.schematisation.schematisation.substance:
            # Compare the dictionaries and override the default values if the user inputs a value
            # assumes that default dict contains the substance input by the user (we only have three right now though!)
            default_substance_dict = self.substance.substance_dict
            user_substance_dict = self.schematisation.schematisation.substance_parameters #user input dictionary of values

            # iterate through the dicitonary keys
            for key, value in user_substance_dict .items():
                if type(value) is dict:
                    for tkey, cvalue in value.items():
                        if cvalue is None: #reassign the value from the default dict if not input by the user
                            user_substance_dict[key][tkey] = default_substance_dict[key][tkey] 
                else:
                    if value is None:
                        user_substance_dict [key] = default_substance_dict[key]

            self.substance_dict = user_substance_dict #assign updated dict as attribute of the class to be able to access later
        else:
            self.substance_dict = self.substance.substance_dict

        # self.df_flowline['substance'] = self.substance_dict['substance_name']

    def _init_omp(self):
        if self.omp_inialized:
            pass
        else:
            self.df_particle['omp_half_life'] = self.df_particle['redox_zone'].map(self.substance_dict['omp_half_life'])
            self.df_particle['log_Koc'] = self.substance_dict['log_Koc']
            self.df_particle['pKa'] = self.substance_dict['pKa']

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

        if self.schematisation.schematisation.biodegradation_sorbed_phase:
            self.df_particle['retardation'] = (1 + (1 / (1 + 10 ** (self.df_particle.pH - self.df_particle.pKa)) * self.df_particle.solid_density_layer
                            * (1 - self.df_particle.porosity_layer)
                            * self.df_particle.fraction_organic_carbon * self.df_particle.Koc_temperature_correction)
                    / (self.df_particle.porosity_layer * (1 + (self.df_particle.Koc_temperature_correction * 1 / (1 + 10 ** (self.df_particle.pH - self.df_particle.pKa))
                                                            * 0.2 * self.df_particle.dissolved_organic_carbon * 0.000001))))
        else:
            self.df_particle['retardation'] = 1

    def _calculate_omp_half_life_temperature_correction(self):
        '''Equation 3.2 in report
        R = 8.314 J/K/mol
        Ea = activation energy = 63*10^3 J/mol'''
        
        #@MartinK this seems a bit roundabout way to access this?
        if self.schematisation.schematisation.temp_correction_halflife:
            self.df_particle['omp_half_life_temperature_corrected'] = self.df_particle['omp_half_life'] * 10 ** (-63000 / (2.303 * 8.314) * (1 / (20 + 273.15) - 1 / (self.df_particle.temperature + 273.15)))
        else: 
            self.df_particle['omp_half_life_temperature_corrected'] = self.df_particle['omp_half_life']

        self.df_particle.loc[ self.df_particle.omp_half_life == 1e99, 'omp_half_life_temperature_corrected'] = 1e99

    def _calculate_Koc_temperature_correction(self):
        ''' Equation 3.1 in report, 
        from Luers and Ten Hulscher (1996): Assuming the relation to be similar 
        to the Van ‘t Hoff equation and equally performing for other OMPs yields'''

        # if log_Koc is zero, assign value of zero
        if self.df_particle.log_Koc[0] == 0:
            self.df_particle['Koc_temperature_correction'] = 0
        elif self.schematisation.schematisation.temp_correction_Koc:
            self.df_particle['Koc_temperature_correction'] = 10 ** self.df_particle.log_Koc * 10 ** (1913 * (1 / (self.df_particle.temperature + 273.15) - 1 / (20 + 273.15)))
        else: 
            self.df_particle['Koc_temperature_correction'] = self.df_particle.log_Koc
           
    def _calculate_state_concentration_in_zone(self):
        '''Equation 4.11 in report '''

        #check if there is degradation prior to infiltration
        DOC_inf = self.schematisation.schematisation.dissolved_organic_carbon_infiltration_water 
        TOC_inf = self.schematisation.schematisation.total_organic_carbon_infiltration_water 

        if DOC_inf and TOC_inf > 0:
            DOC_TOC_ratio = DOC_inf / TOC_inf  
            K_oc = self.df_particle['Koc_temperature_correction'].iloc[0]
            c_in = 100 - 100 * (1 - (DOC_TOC_ratio + (1 - DOC_TOC_ratio) / (1 + K_oc * TOC_inf  * 0.000001)))
            self.df_particle.loc[self.df_particle.zone=='surface', 'input_concentration']=c_in

        # if self.schematisation.schematisation.concentration_point_contamination is None:
        #     ''' diffuse contamination'''
        for i in range(len(self.df_particle)-1):
            if self.df_particle.steady_state_concentration.loc[i+1] is None:
                
                # if omp is persistent, value at end of zone equal to value incoming to zone
                if self.df_particle.omp_half_life.loc[i+1] == 1e99:
                    self.df_particle.at[i+1, 'steady_state_concentration'] = self.df_particle.steady_state_concentration.loc[i]
                
                # AH 300 limit only to avoid very small numnbers, makes no difference for other calculations therefore removed
                # Column O in Phreatic excel sheet
                # Put back in, otherwise there is an error there are too many numbers in the output
                # @MartinK or @MartinvdS alternative to the error? or this method?
                elif (self.df_particle.travel_time_zone.loc[i+1] * self.df_particle.retardation.loc[i+1]
                                                                            / self.df_particle.omp_half_life_temperature_corrected.loc[i+1]) >300:
                    self.df_particle.at[i+1, 'steady_state_concentration'] = 0

                # otherwise, calculate the outcoming concentration from the zone, given the input concentration to the zone. 
                # in the case of the vadose zone, the incoming concentration is the initial concentration
                else:
                    self.df_particle.at[i+1, 'steady_state_concentration'] = (self.df_particle.steady_state_concentration.loc[i]
                                                                            / (2 ** (self.df_particle.travel_time_zone.loc[i+1] * self.df_particle.retardation.loc[i+1]
                                                                            / self.df_particle.omp_half_life_temperature_corrected.loc[i+1])))
        # else:
        #     '''point contamination'''


    def _calculcate_total_breakthrough_travel_time(self):
        ''' Calculate the total time (days) for breakthrough at the well'''
        self.df_flowline['total_breakthrough_travel_time']  = ""
        self.df_flowline['breakthrough_concentration']  = ""

        for i in range(len(self.df_flowline)):
            flowline_id = i + 1

            df = self.df_particle.loc[self.df_particle['flowline_id'] == flowline_id]
            df.fillna(0)['breakthrough_travel_time']
            self.df_flowline.at[i, 'total_breakthrough_travel_time'] = sum(df.fillna(0)['breakthrough_travel_time'])
            self.df_flowline.at[i, 'breakthrough_concentration'] = df['steady_state_concentration'].iloc[-1]

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
        
        self.df_flowline['input_concentration'] = self.schematisation.schematisation.diffuse_input_concentration
        self.df_particle['input_concentration'] = None
        self.df_particle['steady_state_concentration'] = None

        self.df_particle.loc[self.df_particle.zone=='surface', 'input_concentration'] = self.schematisation.schematisation.diffuse_input_concentration
        self.df_particle.loc[self.df_particle.zone=='surface', 'steady_state_concentration'] = self.schematisation.schematisation.diffuse_input_concentration

        if self.schematisation.schematisation.concentration_point_contamination:
            ''' point contamination '''
            # need to take into account the depth of the point contamination here....
            # need to change the df_particle and df_flowline to only be the flowlines for the point contamination flowline(s)
            # use a single point contamination for now
            # FIRST recalculate the travel times for the contamination, then initialize the class

            #only for a SINGLE point contamination
            distance = self.schematisation.schematisation.distance_point_contamination_from_well
            depth = self.schematisation.schematisation.depth_point_contamination
            cumulative_fraction_abstracted_water = (math.pi * self.schematisation.schematisation.recharge_rate 
                                                        * distance ** 2)/self.schematisation.schematisation.well_discharge
            ind = self.df_particle.flowline_id.iloc[-1]

            if self.schematisation.schematisation.schematisation_type == 'phreatic':
                self.head = self.schematisation.schematisation.calculate_hydraulic_head_phreatic(distance=distance)
                self.schematisation.phreatic(distance = distance, 
                                            depth_point_contamination=depth, 
                                            cumulative_fraction_abstracted_water=cumulative_fraction_abstracted_water )

            elif self.schematisation.schematisation.schematisation_type == 'semiconfined':
                bottom_vadose_zone = self.schematisation.schematisation.bottom_vadose_zone_at_boundary
                
                df_flowline, df_particle = self.schematisation.add_semiconfined_point_sources(distance=distance, 
                                        depth_point_contamination=depth,  )

            df_particle['flowline_id'] = df_particle['flowline_id'] + ind
            
            df_flowline['input_concentration'] = self.schematisation.schematisation.concentration_point_contamination
            df_particle['input_concentration'] = None
            df_particle['steady_state_concentration'] = None
            df_particle.loc[self.df_particle.zone=='surface', 'input_concentration'] = self.schematisation.schematisation.concentration_point_contamination
            df_particle.loc[self.df_particle.zone=='surface', 'steady_state_concentration'] = self.schematisation.schematisation.concentration_point_contamination

            df_flowline['flowline_id'] = df_flowline['flowline_id'] + ind
            df_flowline['flowline_type'] = "point_source"
            df_flowline['discharge'] = self.schematisation.schematisation.discharge_point_contamination

            #AH_todo, something here to loop through the different point sources?

            self.df_particle = self.df_particle.append(df_particle)
            self.df_particle.reset_index(drop=True, inplace=True)

            self.df_flowline = self.df_flowline.append(df_flowline) #
            self.df_flowline.reset_index(drop=True, inplace=True)

            self.df_flowline['substance'] = self.substance_dict['substance_name']

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


    def plot_concentration(self, xlim=[0, 500], ylim=[0,1 ]):
        ''' Plot the concentration of the given OMP as a function of time since the start of the contamination'''
        time_array = np.arange(0, 505, 1)*365.24
        
        # reduce the amount of text per line by extracting the following parameters
        concentration_point_contamination = self.schematisation.schematisation.concentration_point_contamination
        diffuse_input_concentration = self.schematisation.schematisation.diffuse_input_concentration
        schematisation_type = self.schematisation.schematisation.schematisation_type

        if concentration_point_contamination is None:
            input_concentration = diffuse_input_concentration 
        else:
            input_concentration = diffuse_input_concentration + concentration_point_contamination

        #Calculate the concentration in the well, 
        self.df_flowline['concentration_in_well'] = (self.df_flowline['breakthrough_concentration'] 
                            * self.df_flowline['discharge']/ self.df_flowline['well_discharge'])
        well_conc = []

        #sum the concentration in the well for each timestep
        for i in range(len(time_array)):
            t = time_array[i]
            well_conc.append(sum(self.df_flowline['concentration_in_well'].loc[self.df_flowline['total_breakthrough_travel_time'] <= t]))

        # well_conc = well_conc/input_concentration
        well_conc[:] = [x / input_concentration for x in well_conc]
        fig = plt.figure(figsize=[10, 5])
        plt.plot(time_array/365.24, well_conc, 'b', label =str(self.substance.substance_name))
        plt.xlim(xlim)
        plt.ylim(ylim) 
        plt.ylabel('Fraction of input concentration')
        plt.xlabel('Time since start of contamination (years)')
        plt.title('Aquifer type: ' + schematisation_type)
        plt.grid()
        plt.legend()
        plt.savefig('well_concentration_over_time_'+str(self.substance.substance_name)+'_'+schematisation_type+'.png', dpi=300, bbox_inches='tight')  # save_results_to + '/


    def plot_age_distribution(self):
        #AH_todo
        pass

    def plot_logremoval(self):
        #AH_todo
        pass
