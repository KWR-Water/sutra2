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

#### Notes ####

# things which must be checked indicated in comments with AH
# specific questions flagged for;
# @MartinvdS // @steven //@martinK

####

#%% ----------------------------------------------------------------------------
# INITIALISATION OF PYTHON e.g. packages, etc.
# ------------------------------------------------------------------------------

import matplotlib.pyplot as plt
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

module_path = os.path.abspath(os.path.join("..","..","greta"))
if module_path not in sys.path:
    sys.path.insert(0,module_path)

from greta.Analytical_Well import AnalyticalWell 
from greta.ModPath_functions import ModPathWell

# from Analytical_Well import AnalyticalWell
# from ModPath_functions import ModPathWell
# # from greta.ModPath_functions import ModPathWell
# # from Analytical_Well import AnalyticalWell

path = os.getcwd()  # path of working directory

class Substance:
    ''' 
    Placeholder class which will later be replaced by the QSAR functionality of AquaPriori.

    Attributes
    ---------
    substance_name: String,
        substance_name of the substance (for now limited dictionary to 'benzene', 'AMPA', 'benzo(a)pyrene'
    substance_dict: dictionary
        Nested dictionary with the following per substance.
        substance_name: String,
            substance_name of the substance (for now limited dictionary to 'benzene', 'AMPA', 'benzo(a)pyrene'
        log Koc: float
            distribution coefficient of organic carbon and water, [-]
        molar_mass: float
            molar mass of substance, [g/mol]
        pKa: float
            disassociation constant for acid H-OMP, [-]
        omp_half_life: float
            per redox zone ('suboxic', 'anoxic', deeply_anoxic'), [days]
    '''
    def __init__(self, substance_name, ):
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

        @SR (21-1-2022): ADD documentation microbial species       
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
            "norovirus": 
                {"substance_name": "norovirus",
                    "alpha0": 1.e-5,
                    "pH0": 6.8,
                    "pathogen_diam": 2.33e-8,
                    "mu1": {
                        "suboxic": 0.149, 
                        "anoxic": 0.023, 
                        "deeply_anoxic": 0.023
                        }
                }
            }
        self.substance_dict = substances_dict[substance_name]

        # Dict    # Acties Steven 25-1-22
        micro_organism_dict = {}

        #@ Steven voeg toe: micro_organism_dict
        self.micro_organism_dict = micro_organism_dict[substance_name]

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
        Column 'flowline_id': Integer
        Column 'flowline_type': string
        Column 'flowline_discharge': Float
        Column 'particle_release_day': Float
        Column 'input_concentration': float
        Column 'endpoint_id': Integer
        Column 'well_discharge': float
        Column 'substance': string
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
        Column 'temperature': float
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

    """


    def __init__(self,
                well: AnalyticalWell or ModPathWell,
                substance: Substance):

        '''
        Initialization of the Substanes class, checks for user-defined OMP substance paramters and overrides the database values.

        Parameters
        ----------
        well: object
            The AnalyticalWell object for the schematisation of the aquifer type.
        substance: object
            The Substance object with the OMP of interest.


        '''
        self.well = well
        self.omp_inialized = False
        self.df_particle = well.df_particle
        self.df_flowline = well.df_flowline
        self.substance = Substance(substance)

        # PATHOGEN init
        # @SR : init of microbial species (pathogen)
        # Check if ModPathWell (Class) can be used instead of 'AnalyticalWell'
        
        self.micro_organism_inialized = False

        # AH need to make sure here that the substance passed is the same, e.g. comapre the dictionaries BUT ALSO
        # make sure that user doesn't call one substance in the hydrochemicalschematisation class and another in the concentration class
        # probably only a problem for ourselves, this should be written into a larger "run" class for the model which could avoid this
        if self.substance.substance_name == self.well.schematisation.substance:
            # Compare the dictionaries and override the default values if the user inputs a value
            # assumes that default dict contains the substance input by the user (we only have three right now though!)
            default_substance_dict = self.substance.substance_dict
            user_substance_dict = self.well.schematisation.substance_parameters #user input dictionary of values

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
        ''' 
        Initialisation if the Substance is an OMP
        '''
        if self.omp_inialized:
            pass
        else:
            self.df_particle['omp_half_life'] = self.df_particle['redox'].map(self.substance_dict['omp_half_life'])
            self.df_particle['log_Koc'] = self.substance_dict['log_Koc']
            self.df_particle['pKa'] = self.substance_dict['pKa']

        self.omp_inialized = True


    def _init_micro_organism(self):
        ''' Initialisation if the Substance is a pathogen'''
        
        if self.micro_organism_inialized:
            pass
        else:
            self.df_particle['mu1'] = self.df_particle['redox'].map(self.substance_dict['mu1'])
            self.df_particle['alpha0'] = self.substance_dict['alpha0']
            self.df_particle['pH0'] = self.substance_dict['pH0']
            self.df_particle['pathogen_diam'] = self.substance_dict['pathogen_diam']

        self.pathogen_inialized = True


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
            self.df_particle['omp_half_life_temperature_corrected'] = self.df_particle['omp_half_life'] * 10 ** (-63000 / (2.303 * 8.314) * (1 / (20 + 273.15) - 1 / (self.df_particle.temperature + 273.15)))
        else:
            self.df_particle['omp_half_life_temperature_corrected'] = self.df_particle['omp_half_life']

        self.df_particle.loc[ self.df_particle.omp_half_life == 1e99, 'omp_half_life_temperature_corrected'] = 1e99

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


        # if log_Koc is zero, assign value of zero
        if self.df_particle.log_Koc[0] == 0:
            self.df_particle['Koc_temperature_correction'] = 0
        elif self.well.schematisation.temp_correction_Koc:
            self.df_particle['Koc_temperature_correction'] = 10 ** self.df_particle.log_Koc * 10 ** (1913 * (1 / (self.df_particle.temperature + 273.15) - 1 / (20 + 273.15)))
        else:
            self.df_particle['Koc_temperature_correction'] = self.df_particle.log_Koc

    def _calculate_state_concentration_in_zone(self):
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

        for i in range(len(self.df_particle)-1):
            if self.df_particle.steady_state_concentration.loc[i+1] is None:

                # if omp is persistent, value at end of zone equal to value incoming to zone
                if self.df_particle.omp_half_life.loc[i+1] == 1e99:
                    self.df_particle.at[i+1, 'steady_state_concentration'] = self.df_particle.steady_state_concentration.loc[i]

                # Column O in Phreatic excel sheet
                # # AH 300 limit only to avoid very small numnbers, makes no difference for other calculations therefore left in
                # Put back in, otherwise there is an error there are too many numbers in the output
                # can't reproduce the error, so take out again
                
                elif (self.df_particle.travel_time.loc[i+1] * self.df_particle.retardation.loc[i+1]
                                                                            / self.df_particle.omp_half_life_temperature_corrected.loc[i+1]) >300:
                    self.df_particle.at[i+1, 'steady_state_concentration'] = 0

                # otherwise, calculate the outcoming concentration from the zone, given the input concentration to the zone.
                # in the case of the vadose zone, the incoming concentration is the initial concentration
                else:
                    self.df_particle.at[i+1, 'steady_state_concentration'] = (self.df_particle.steady_state_concentration.loc[i]
                                                                            / (2 ** (self.df_particle.travel_time.loc[i+1] * self.df_particle.retardation.loc[i+1]
                                                                            / self.df_particle.omp_half_life_temperature_corrected.loc[i+1])))

    def _calculcate_total_breakthrough_travel_time(self):
        ''' Calculate the total time for breakthrough for each flowline at the well

        Returns
        -------
        df_flowline: pandas.dataframe
            Column 'total_breakthrough_travel_time': float

        '''
        self.df_flowline['total_breakthrough_travel_time']  = ""
        self.df_flowline['breakthrough_concentration']  = ""

        for i in range(len(self.df_flowline)):
            flowline_id = i + 1

            df = self.df_particle.loc[self.df_particle['flowline_id'] == flowline_id]
            df.fillna(0)['breakthrough_travel_time']
            self.df_flowline.at[i, 'total_breakthrough_travel_time'] = sum(df.fillna(0)['breakthrough_travel_time'])
            self.df_flowline.at[i, 'breakthrough_concentration'] = df['steady_state_concentration'].iloc[-1]

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
            Column 'substance': string
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
            Column 'temperature': float
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

        self.df_flowline['input_concentration'] = self.well.schematisation.diffuse_input_concentration
        self.df_particle['input_concentration'] = None
        self.df_particle['steady_state_concentration'] = None

        self.df_particle.loc[self.df_particle.zone=='surface', 'input_concentration'] = self.well.schematisation.diffuse_input_concentration
        self.df_particle.loc[self.df_particle.zone=='surface', 'steady_state_concentration'] = self.well.schematisation.diffuse_input_concentration

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
            ind = self.df_particle.flowline_id.iloc[-1]

            if self.well.schematisation.schematisation_type == 'phreatic':
                head = self.well.schematisation._calculate_hydraulic_head_phreatic(distance=distance)
                df_flowline, df_particle = self.well._add_phreatic_point_sources(distance=distance,
                                            depth_point_contamination=depth,
                                            cumulative_fraction_abstracted_water=cumulative_fraction_abstracted_water)

            elif self.well.schematisation.schematisation_type == 'semiconfined':
                bottom_vadose_zone = self.well.schematisation.bottom_vadose_zone_at_boundary

                df_flowline, df_particle = self.well._add_semiconfined_point_sources(distance=distance,
                                        depth_point_contamination=depth,  )

            df_particle['flowline_id'] = df_particle['flowline_id'] + ind

            df_flowline['input_concentration'] = self.well.schematisation.point_input_concentration
            df_particle['input_concentration'] = None
            df_particle['steady_state_concentration'] = None
            df_particle.loc[self.df_particle.zone=='surface', 'input_concentration'] = self.well.schematisation.point_input_concentration
            df_particle.loc[self.df_particle.zone=='surface', 'steady_state_concentration'] = self.well.schematisation.point_input_concentration

            df_flowline['flowline_id'] = df_flowline['flowline_id'] + ind
            df_flowline['flowline_type'] = "point_source"
            df_flowline['flowline_discharge'] = abs(self.well.schematisation.discharge_point_contamination)

            #AH_todo, something here to loop through different point sources if more than one.

            self.df_particle = self.df_particle.append(df_particle)
            self.df_particle.reset_index(drop=True, inplace=True)

            self.df_flowline = self.df_flowline.append(df_flowline)
            self.df_flowline.reset_index(drop=True, inplace=True)

            self.df_flowline['substance'] = self.substance_dict['substance_name']

        self._init_omp()

        self._calculate_Koc_temperature_correction()

        self._calculate_omp_half_life_temperature_correction()

        self._calculate_retardation()

        self._calculate_state_concentration_in_zone()

        self.df_particle['breakthrough_travel_time'] = self.df_particle.retardation * self.df_particle.travel_time

        self._calculcate_total_breakthrough_travel_time()

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
        
        # add the particle release date
        self.df_flowline['particle_release_day'] = (self.start_date - start_date_contamination).days


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
                value_mean = df[df_column].values.mean()
                df[df_column] = df[df_column].fillna(value_mean)
        # Adjust dtype
        df = df.astype({df_column: dtype_})

        return df

    def calc_lambda(self, df_particle, df_flowline, mu1 = 0.149, mu1_std = 0.0932,
                por_eff = 0.33, 
                grainsize = 0.00025,
                alpha0 = 0.001,
                pH0 = 7., 
                const_BM = 1.38e-23,
                temp_water = 10., rho_water = 999.703,
                pathogen_diam = 2.33e-8, v_por = 0.01):

        ''' For more information about the advective microbial removal calculation: 
                BTO2012.015: Ch 6.7 (page 71-74)

        Calculate removal coefficient 'lambda_' [/day].

            lambda_ = k_att + mu_1
            with 'hechtingssnelheidscoëfficiën't k_att [/day] and
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
        # Pathogen diameter 'pathogen_diam' [m]
        # porewater velocity 'v_por' [m/d]
        
        # Check - (example 'E_coli'):
        >> k_att = calc_katt(part_idx = [0], por_eff = [0.33], korrelgrootte = [0.00025],
                bots_eff = 0.001, const_BM = 1.38e-23,
                temp_water = [10.], rho_water = [999.703],
                pathogen_diam = 2.33e-8, v_por = [0.01])
        >> print(k_att)
        >> {0: 0.7993188853572424} # [/day]
        
        Return
        ---------
        df_flowline: pandas.DataFrame
            Column "collision_eff" for each node
            Column "pathogen_diam" for each node

        df_particle: pandas.DataFrame
            Column "relative_distance" for each node
            Column "porewater_velocity" for each node
            Column "k_att" for each node
            Column "lambda" for each node
        '''

        # Create empty column 'relative_distance' in df_particle
        df_particle['relative_distance'] = None  
        # Create empty column 'porewater_velocity' in df_particle
        df_particle['porewater_velocity'] = None     

        # Create empty column 'collision_eff' in df_particle
        df_particle['collision_eff'] = None  
        # Collission efficiency [np.array]
        coll_eff = {}
        for pid in df_flowline.index:
            coll_eff[pid] = alpha0 * 0.9**((df_particle.loc[pid,"pH"].values - pH0)/0.1)

            # Fill df_particle 'collision_eff'
            df_particle.loc[pid,"collision_eff"] = coll_eff[pid]
        
        # Collision efficiency [-]
        df_flowline = self._df_fillna(df_flowline,
                                      df_column = 'collision_eff',
                                      value = alpha0)
        # Pathogen diameter [m]
        df_flowline = self._df_fillna(df_flowline,
                                df_column = 'pathogen_diam',
                                value = pathogen_diam)

        # Water temperature [degrees Celsius]
        df_particle = self._df_fillna(df_particle,
                                df_column = 'temperature',
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

        # Create empty column 'k_att' in df_particle
        df_particle['k_att'] = None
        # Create empty column 'lambda' in df_particle
        df_particle['lambda'] = None        
        
        # if not 'collision_eff' in df_flowline.columns:
        #     df_flowline['collision_eff'] = coll_eff
        # else:
        #     # fill series (if needed with default value)
        #     if df_flowline['collision_eff'].dropna().empty:
        #         df_flowline['collision_eff'] = df_flowline['collision_eff'].fillna(coll_eff)
        #     else: # fill empty rows (with mean value of other records)
        #         coll_eff_mean = df_flowline['collision_eff'].values.mean()
        #         df_flowline['collision_eff'] = df_flowline['collision_eff'].fillna(coll_eff_mean)
   
        # # Pathogen diameter [m]
        # if not 'pathogen_diam' in df_flowline.columns:
        #     df_flowline['pathogen_diam'] = pathogen_diam
        # else:
        #     # fill series (if needed with default value)
        #     if df_flowline['pathogen_diam'].dropna().empty:
        #         df_flowline['pathogen_diam'] = df_flowline['pathogen_diam'].fillna(pathogen_diam)
        #     else: # fill empty rows (with mean value of other records)
        #         pathogen_diam_mean = df_flowline['pathogen_diam'].values.mean()
        #         df_flowline['pathogen_diam'] = df_flowline['pathogen_diam'].fillna(pathogen_diam_mean)

        # # Water temperature
        # if not 'temperature' in df_particle.columns:
        #     df_particle['temperature'] = temp_water
        # else:
        #     # fill series (if needed with default value)
        #     if df_particle['temperature'].dropna().empty:
        #         df_particle['temperature'] = df_particle['temperature'].fillna(temp_water)
        #     else: # fill empty rows (with mean value of other records)
        #         temp_water_mean = df_particle['temperature'].values.mean()
        #         df_particle['temperature'] = df_particle['temperature'].fillna(temp_water_mean)

        # Calculate porewater velocity
        v_por, dist, tdiff = {}, {}, {}   
        for pid in df_flowline.index:


            #Distance array between nodes
            dist[pid] = np.sqrt((df_particle.loc[pid,"xcoord"].diff().values)**2 + \
                                (df_particle.loc[pid,"ycoord"].diff().values)**2 + \
                                (df_particle.loc[pid,"zcoord"].diff().values)**2)
            # Time difference array
            tdiff[pid] = df_particle.loc[pid,"total_travel_time"].diff().values
            # Calculate porewater velocity [m/day] (do not include effective porosity)
            v_por[pid] = abs(dist[pid][1:]/tdiff[pid][1:])

            # Fill column relative_distance in 'df_particle' 
            df_particle.loc[pid,'relative_distance'] = np.array(list(dist[pid][1:]) + [0])
            # df_particle.loc[pid,'relative_distance'][-1] = 0

            # Fill porewater velocity in 'df_particle'
            df_particle.loc[pid,"porewater_velocity"] = np.array(list(v_por[pid]) + [v_por[pid][-1]])
            # In the last pathline row the porewater_velocity is equal to previous velocity
            # df_particle.loc[pid,"porewater_velocity"][-1] = v_por[pid][-1]

            # number of nodes
            n_nodes = len(df_particle.loc[pid,:])
            # #Distance array between nodes
            # dist[pid] = np.zeros((n_nodes-1), dtype = 'float')
            # # Time difference array
            # tdiff[pid] = np.zeros((n_nodes-1), dtype = 'float')
            # for iNode in range(1,n_nodes):
            #     # Calculate internodal distance (m)
            #     dist[pid][iNode-1] = np.sqrt((df_particle.loc[pid,"xcoord"][iNode] - df_particle.loc[pid,"xcoord"][iNode-1])**2 + \
            #         (df_particle.loc[pid,"ycoord"][iNode] - df_particle.loc[pid,"ycoord"][iNode-1])**2 + \
            #         (df_particle.loc[pid,"zcoord"][iNode] - df_particle.loc[pid,"zcoord"][iNode-1])**2)
            #     # Calculate time difference between nodes
            #     tdiff[pid][iNode-1] = (df_particle.loc[pid,"total_travel_time"][iNode] - df_particle.loc[pid,"total_travel_time"][iNode-1])
                

        # Create empty dicts and keys to keep track of arrays per particle id
        k_bots, gamma, As_happ, mu = {}, {}, {}, {}
        D_BM, k_diff, k_att, lambda_ = {}, {}, {}, {}
        # particle ids
        for pid in df_flowline.index:

            # Botsingterm 'k_bots'
            k_bots[pid] = (3/2.)*((1-df_particle.loc[pid,"porosity"].values) / \
                        df_particle.loc[pid,"grainsize"].values) * df_flowline.loc[pid,"collision_eff"]

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
                        (df_particle.loc[pid,"temperature"].values + 42.5)**(3/2)
 
            # Diffusion constant 'D_BM' (Eq.6: BTO2012.015) --> eenheid: m2 s-1
            D_BM[pid] = (const_BM * (df_particle.loc[pid,"temperature"].values + 273.)) / \
                        (3 * np.pi * df_flowline.loc[pid,"pathogen_diam"] * mu[pid])
            # Diffusieconstante 'D_BM' (Eq.6: BTO2012.015) --> eenheid: m2 d-1
            D_BM[pid] *= 86400.

            # Diffusion related attachment term 'k_diff'
            k_diff[pid] = ((D_BM[pid] /
                        (df_particle.loc[pid,"grainsize"].values * df_particle.loc[pid,"porosity"].values * df_particle.loc[pid,"porewater_velocity"].values))**(2/3) * 
                            df_particle.loc[pid,"porewater_velocity"].values)

            # hechtingssnelheidscoëfficiënt k_att [/dag]
            k_att[pid] = k_bots[pid] * 4 * As_happ[pid]**(1/3) * k_diff[pid]
            # removal coefficient 'lambda_' [/day], using the 'mu1' mean.
            lambda_[pid] = k_att[pid] + mu1

            # Fill df_particle 'k_att'
            df_particle.loc[pid,"k_att"] = k_att[pid]
            # Fill df_particle 'lambda'
            df_particle.loc[pid,"lambda"] = lambda_[pid]

        # return (adjusted) df_particle and df_flowline
        return df_particle, df_flowline
        
    def calc_advective_microbial_removal(self,df_particle,df_flowline, 
                                        endpoint_id, trackingdirection = "forward",
                                        mu1 = 0.023, grainsize = 0.00025, alpha0 = 1.E-5, pH0 = 6.8, const_BM = 1.38e-23,
                                        temp_water = 11., rho_water = 999.703, pathogen_diam = 2.33e-8,
                                        conc_start = 1., conc_gw = 0.):
                                        
        ''' Calculate the advective microbial removal along pathlines
            from source to end_point.

            For more information about the advective microbial removal calculation: 
                BTO2012.015: Ch 6.7 (page 71-74)

            Relevant input:
            df_particle: pandas.DataFrame
                Column 'flowline_id'
                Column 'xcoord'
                Column 'ycoord'
                Column 'zcoord'
                Column 'relative_distance'
                Column 'total_travel_time'
                Column 'porosity'
                Column 'lambda'
                Column 'porewater_velocity'

            df_flowline: pandas.DataFrame
                Column 'flowline_id'
                Column 'endpoint_id'
                Column 'flowline_discharge'
                Column 'well_discharge'

            Return 
            df_flowline: pandas.DataFrame
                Column 'input_concentration': float   # (added)

            df_particle: pandas.DataFrame
                Column 'steady_state_concentration': float      # (added)
                Column 'starting_concentration_gw': float                         # (added)
        

            Calculate the steady state concentration along traveled distance per 
            node for each pathline from startpoint to endpoint_id'.
            
            With verwijderingscoëfficiënt 'lambda_' [/d]
            effective porositty 'porosity' [-]
            Starting concentration 'input_concentration' per pathline
            Initial groundwater concentration 'starting_concentration_gw'

            # Algemeen geldt voor de afbraak vanaf maaiveld tot eindpunt (zoals 'lek')
            # C_ = (C_mv-C_gw) * e^-(lambda_/v)*dr + C_gw
            # C_mv = C_maaiveld [aantal/l]; C_gw = C_grondwater [aantal/l]?
            # v = (dr/dt) / por_eff --> poriewatersnelheid [m/d]
            # Bij achtergrondconcentratie C_gw = 0., zoals nu:
            # C_ = C_mv * e^-(lambda_/v)*r

            # Check eindconcentratie van deze tijdstap
            >> _,C_final = calc_pathlineconc(part_idx = [0], C0 = 1., C_gw = 0., por_vals = 0.33,
                        lambda_ = {0: 0.1317345756973715}, xyz_start = [0,0,0], 
                        xyz_eind = [0,0,1], t_start = 0., t_eind = 100.,
                        dist_data = None, time_diff = None,
                        direction = "backward")
            >> print (C_final) [aantal/l]
            >> 1.900378331566034e-06
            
            # Voor berekening per particle path "part_idx"
            # Afgelegde afstanden tussen "xyz_node t" en "xyz_nodel t+1"
            # voor iedere tijdstap (--> path_data)
            dist_data: {0: np.array([1,1.2,1.3]),n+1: np.array([...])}
            # Tijdsduur tussen iedere tijdstap "t" en "t+1"
            # Zie parameter --> time_data
            time_data: {0: np.array([2.2,2.4,1.6]),n+1: np.array([...])}
            # direction: particle tracking direction (forward/backward)

        '''

        # Calculate removal coefficient 'lambda' [/day].
        df_particle, df_flowline = self.calc_lambda(df_particle, df_flowline,
                                                    mu1 = mu1, grainsize = grainsize, alpha0 = alpha0, 
                                                    pH0 = pH0, const_BM = const_BM,
                                                    temp_water = temp_water, 
                                                    rho_water = rho_water, pathogen_diam = pathogen_diam)

        # Starting concentration [-]
        df_flowline = self._df_fillna(df_flowline,
                                df_column = 'starting_concentration',
                                value = conc_start,
                                dtype_ = 'float')

        # Original concentration in groundwater [-]
        df_flowline = self._df_fillna(df_flowline,
                                df_column = 'starting_concentration_gw',
                                value = conc_gw,
                                dtype_ = 'float')

        if 'steady_state_concentration' not in df_particle.columns:
            # Steady state concentration [-]
            df_particle['steady_state_concentration'] = None
        
        for pid in df_flowline.index:
            df_particle.loc[pid,"steady_state_concentration"].iloc[0] = df_flowline.loc[pid,"starting_concentration"]

        # Start dictionaries for relative and steady-state concentrations
        conc_rel = {}
        conc_steady = {}  # concentration corrected for initial concentration
        for pid in df_flowline.index:
            
            # Calculate the relative removal along pathlines
            exp_arg = -((df_particle.loc[pid,"lambda"].values/df_particle.loc[pid,"porewater_velocity"].values) *
                                df_particle.loc[pid,'relative_distance'].values).astype('float')

            conc_rel[pid] = df_flowline.loc[pid,"starting_concentration_gw"] + \
                            (df_flowline.loc[pid,"starting_concentration"] - df_flowline.loc[pid,"starting_concentration_gw"]) * \
                                np.exp(exp_arg)

            if trackingdirection == 'forward':
                conc_steady[pid] = (conc_rel[pid]/df_flowline.loc[pid,"starting_concentration"]).cumprod() * \
                                    df_flowline.loc[pid,"starting_concentration"]
                # Replace 'nan'-values by 0.
                conc_steady[pid]  = np.nan_to_num(conc_steady[pid],0.)
                
                # Index of average final concentration "C_"  idx=-1[C0,...,..,...,Ct-1,C_final]
                idx_final = -1
                                
            elif trackingdirection == 'backward':
                conc_steady[pid] = (conc_rel[pid][::-1]/df_flowline.loc[pid,"starting_concentration"]).cumprod() * \
                                    df_flowline.loc[pid,"starting_concentration"]
                # Replace 'nan'-values by 0.
                conc_steady[pid]  = np.nan_to_num(conc_steady[pid],0.)
                # First value represents concentration of endpoint and v.v.
                # Reverse array--> [C_final,Ct-1,..,...,C0]
                conc_steady[pid] = conc_steady[pid][::-1]
                # Index of average final concentration "C_" idx=0 [C_final,Ct-1,..,...,C0]
                idx_final = 0

            # Add steady_state_concentration to df_particle
            df_particle.loc[pid,'steady_state_concentration'] = conc_steady[pid]

        # Bereken gemiddelde eindconcentratie
        C_final = 0.
        for pid in df_flowline.index:
            if df_flowline.loc[pid,"endpoint_id"] == endpoint_id:
                C_final += (conc_steady[pid][idx_final] * df_flowline.loc[pid,"flowline_discharge"]) / \
                            df_flowline.loc[pid,"well_discharge"]


        # return (adjusted) df_particle and df_flowline
        return df_particle, df_flowline, C_final    

    def compute_pathogen_removal(self):
        #AH_todo
        pass

    def plot_age_distribution(self):
        #AH_todo
        pass

    def plot_logremoval(self):
        #AH_todo
        pass
