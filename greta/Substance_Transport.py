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

#@MartinK alternative to self.schematisaiton.schematisation... this becomes cumbersome the more classes we have

####

#%% ----------------------------------------------------------------------------
# INITIALISATION OF PYTHON e.g. packages, etc.
# ------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from pandas import read_csv
from pandas import read_excel
from tqdm import tqdm  # tqdm gives a progress bar for the simultation
import math
from scipy.special import kn as besselk
import datetime
from datetime import timedelta  

path = os.getcwd()  # path of working directory

class Substance:
    ''' Placeholder class which will later be replaced by the QSAR functionality of AquaPriori. 
    
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

    Attributes
    ----------
    analytical_well: object
        The AnalyticalWell object for the schematisation of the aquifer type.
    omp_inialized: bool
        Boolian indicating whether the Substance object has been initialized
    df_flowline: pandas.DataFrame
        Column 'flowline_id': Integer
        Column 'flowline_type': string
        Column 'discharge': Float
        Column 'particle_release_date': Float
        Column 'input_concentration': float
        Column 'endpoint_id': Integer
        Column 'well_discharge': float
        Column 'recharge_rate': float
        Column 'vertical_resistance_aquitard': float
        Column 'KD': float
        Column 'thickness_full_capillary_fringe': float
        Column 'substance': string
        Column 'moisture_content_vadose_zone': float
        Column 'diameter_borehole': float
        Column 'removal_function': string
        Column 'solid_density_vadose_zone': float
        Column 'solid_density_shallow_aquifer': float
        Column 'solid_density_target_aquifer': float
        Column 'total_breakthrough_travel_time': float
        Column 'breakthrough_concentration': float

    df_particle: pandas.DataFrame
        Column 'flowline_id': int
        Column 'zone': string
        Column 'travel_time_zone': float
        Column 'xcoord': float
        Column 'ycoord': float
        Column 'zcoord': float
        Column 'redox_zone': float
        Column 'temperature': float
        Column 'thickness_layer': float
        Column 'porosity_layer': float
        Column 'dissolved_organic_carbon': float
        Column 'pH': float
        Column 'fraction_organic_carbon': float
        Column 'solid_density_layer': float
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
                analytical_well, 
                substance: Substance):

        '''
        Initialization of the Substanes class, checks for user-defined OMP substance paramters and overrides the database values.

        Parameters
        ----------
        analytical_well: object
            The AnalyticalWell object for the schematisation of the aquifer type.
        substance: object
            The Substance object with the OMP of interest.

        Attributes
        ----------
        @MartinK these are the same as the class, include again?

        '''
        self.analytical_well = analytical_well
        self.omp_inialized = False
        self.df_particle = analytical_well.df_particle
        self.df_flowline = analytical_well.df_flowline
        self.substance = Substance(substance) 


        # AH need to make sure here that the substance passed is the same, e.g. comapre the dictionaries BUT ALSO
        # make sure that user doesn't call one substance in the hydrochemicalschematisation class and another in the concentration class
        # probably only a problem for ourselves, this should be written into a larger "run" class for the model which could avoid this
        if self.substance.substance_name == self.analytical_well.schematisation.substance:
            # Compare the dictionaries and override the default values if the user inputs a value
            # assumes that default dict contains the substance input by the user (we only have three right now though!)
            default_substance_dict = self.substance.substance_dict
            user_substance_dict = self.analytical_well.schematisation.substance_parameters #user input dictionary of values

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
        ''' Initialisation if the Substance is an OMP'''
        if self.omp_inialized:
            pass
        else:
            self.df_particle['omp_half_life'] = self.df_particle['redox_zone'].map(self.substance_dict['omp_half_life'])
            self.df_particle['log_Koc'] = self.substance_dict['log_Koc']
            self.df_particle['pKa'] = self.substance_dict['pKa']

        self.omp_inialized = True


    def _init_pathogen():
        ''' Initialisation if the Substance is a pathogen'''

        pass

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
        if self.analytical_well.schematisation.biodegradation_sorbed_phase:
            self.df_particle['retardation'] = (1 + (1 / (1 + 10 ** (self.df_particle.pH - self.df_particle.pKa)) * self.df_particle.solid_density_layer
                            * (1 - self.df_particle.porosity_layer)
                            * self.df_particle.fraction_organic_carbon * self.df_particle.Koc_temperature_correction)
                    / (self.df_particle.porosity_layer * (1 + (self.df_particle.Koc_temperature_correction * 1 / (1 + 10 ** (self.df_particle.pH - self.df_particle.pKa))
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
        
        #@MartinK this seems a bit roundabout way to access this?
        if self.analytical_well.schematisation.temp_correction_halflife:
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
        elif self.analytical_well.schematisation.temp_correction_Koc:
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
        DOC_inf = self.analytical_well.schematisation.dissolved_organic_carbon_infiltration_water 
        TOC_inf = self.analytical_well.schematisation.total_organic_carbon_infiltration_water 

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
        """ Calculates the concentration in the well of each flowline. Returns 
        the values in 'df_flowline' and 'df_particle' as attributes of the object.
        
        Returns
        -------
        df_flowline: pandas.DataFrame
            Column 'flowline_id': Integer
            Column 'flowline_type': string
            Column 'discharge': Float
            Column 'particle_release_date': Float
            Column 'input_concentration': float
            Column 'endpoint_id': Integer
            Column 'well_discharge': float
            Column 'recharge_rate': float
            Column 'vertical_resistance_aquitard': float
            Column 'KD': float
            Column 'thickness_full_capillary_fringe': float
            Column 'substance': string
            Column 'moisture_content_vadose_zone': float
            Column 'diameter_borehole': float
            Column 'removal_function': string
            Column 'solid_density_vadose_zone': float
            Column 'solid_density_shallow_aquifer': float
            Column 'solid_density_target_aquifer': float
            Column 'total_breakthrough_travel_time': float
                The breakthrough concentration in the well for the OMP taking into account retardation.
            Column 'breakthrough_concentration': float
                The breakthrough concentration in the well for the OMP taking into account sorption
                and biodegradation.

        df_particle: pandas.DataFrame
            Column 'flowline_id': int
            Column 'zone': string
            Column 'travel_time_zone': float
            Column 'xcoord': float
            Column 'ycoord': float
            Column 'zcoord': float
            Column 'redox_zone': float
            Column 'temperature': float
            Column 'thickness_layer': float
            Column 'porosity_layer': float
            Column 'dissolved_organic_carbon': float
            Column 'pH': float
            Column 'fraction_organic_carbon': float
            Column 'solid_density_layer': float
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
        
        self.df_flowline['input_concentration'] = self.analytical_well.schematisation.diffuse_input_concentration
        self.df_particle['input_concentration'] = None
        self.df_particle['steady_state_concentration'] = None

        self.df_particle.loc[self.df_particle.zone=='surface', 'input_concentration'] = self.analytical_well.schematisation.diffuse_input_concentration
        self.df_particle.loc[self.df_particle.zone=='surface', 'steady_state_concentration'] = self.analytical_well.schematisation.diffuse_input_concentration

        if self.analytical_well.schematisation.concentration_point_contamination:
            ''' point contamination '''
            # need to take into account the depth of the point contamination here....
            # need to change the df_particle and df_flowline to only be the flowlines for the point contamination flowline(s)
            # use a single point contamination for now
            # FIRST recalculate the travel times for the contamination, then initialize the class

            #only for a SINGLE point contamination
            distance = self.analytical_well.schematisation.distance_point_contamination_from_well
            depth = self.analytical_well.schematisation.depth_point_contamination
            cumulative_fraction_abstracted_water = (math.pi * self.analytical_well.schematisation.recharge_rate 
                                                        * distance ** 2)/self.analytical_well.schematisation.well_discharge
            ind = self.df_particle.flowline_id.iloc[-1]

            if self.analytical_well.schematisation.schematisation_type == 'phreatic':
                head = self.analytical_well.schematisation.calculate_hydraulic_head_phreatic(distance=distance)
                df_flowline, df_particle = self.analytical_well.add_phreatic_point_sources(distance=distance, 
                                            depth_point_contamination=depth, 
                                            cumulative_fraction_abstracted_water=cumulative_fraction_abstracted_water)

            elif self.analytical_well.schematisation.schematisation_type == 'semiconfined':
                bottom_vadose_zone = self.analytical_well.schematisation.bottom_vadose_zone_at_boundary
                
                df_flowline, df_particle = self.analytical_well.add_semiconfined_point_sources(distance=distance, 
                                        depth_point_contamination=depth,  )

            df_particle['flowline_id'] = df_particle['flowline_id'] + ind
            
            df_flowline['input_concentration'] = self.analytical_well.schematisation.concentration_point_contamination
            df_particle['input_concentration'] = None
            df_particle['steady_state_concentration'] = None
            df_particle.loc[self.df_particle.zone=='surface', 'input_concentration'] = self.analytical_well.schematisation.concentration_point_contamination
            df_particle.loc[self.df_particle.zone=='surface', 'steady_state_concentration'] = self.analytical_well.schematisation.concentration_point_contamination

            df_flowline['flowline_id'] = df_flowline['flowline_id'] + ind
            df_flowline['flowline_type'] = "point_source"
            df_flowline['discharge'] = self.analytical_well.schematisation.discharge_point_contamination

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

        self.df_particle['breakthrough_travel_time'] = self.df_particle.retardation * self.df_particle.travel_time_zone

        self._calculcate_total_breakthrough_travel_time()

    def compute_concentration_in_well_at_date(self):
        ''' Calculates the concentration in the well up to a specific date, 
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

        # reduce the amount of text per line by extracting the following parameters
        compute_contamination_for_date = self.analytical_well.schematisation.compute_contamination_for_date
        start_date_well = self.analytical_well.schematisation.start_date_well
        start_date_contamination = self.analytical_well.schematisation.start_date_contamination
        end_date_contamination = self.analytical_well.schematisation.end_date_contamination

        if start_date_well > start_date_contamination:
            start_date = start_date_well
            back_date_start = start_date_contamination

        elif start_date_well <= start_date_contamination:
            start_date = start_date_contamination
            back_date_start = start_date_well
        
        compute_date = compute_contamination_for_date - start_date
        back_compute_date = start_date - back_date_start

        # calculate the time after which no more contamination
        if end_date_contamination is None:
            pass
        else:
            end_time = end_date_contamination- start_date
            self.df_flowline['end_time_contamination_breakthrough'] = self.df_flowline['total_breakthrough_travel_time'] + end_time.days

        time_array = np.arange(0, compute_date.days+1, 1)
        back_date_array = np.arange(-back_compute_date.days,0, 1)
        time_array = np.append(back_date_array,time_array)
        time_array_dates = pd.date_range(start=back_date_start,end=compute_contamination_for_date)

        #Calculate the concentration in the well, 
        self.df_flowline['concentration_in_well'] = (self.df_flowline['breakthrough_concentration'] 
                            * self.df_flowline['discharge']/ self.df_flowline['well_discharge'])
        df_flowline = self.df_flowline
        well_concentration = []

        #sum the concentration in the well for each timestep
        for i in range(len(time_array)):
            t = time_array[i] 
            if end_date_contamination is None:
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
            Choice of the x-axis as time in years starting at 0, or as the 'Date' since the 
            minimum of 'start_date_well' or 'start_date_contamination'
        as_fraction_input: bool
            If 'True' plots concentration on y-axis as a fraction of the sum of the 
            input concentration (diffuse and point source), [C/C0]
        xlim: array
            The x-axis limits
        ylim: array
            The y-axis limits

        Returns
        -------
        plot/image @MartinK what to put here?
        '''
        
        # reduce the amount of text per line by extracting the following parameters
        concentration_point_contamination = self.analytical_well.schematisation.concentration_point_contamination
        diffuse_input_concentration = self.analytical_well.schematisation.diffuse_input_concentration
        schematisation_type = self.analytical_well.schematisation.schematisation_type
        compute_contamination_for_date = self.analytical_well.schematisation.compute_contamination_for_date
        start_date_well = self.analytical_well.schematisation.start_date_well
        start_date_contamination = self.analytical_well.schematisation.start_date_contamination
        end_date_contamination = self.analytical_well.schematisation.end_date_contamination

        start_date = max(start_date_well,start_date_contamination)
        back_date_start = min(start_date_well,start_date_contamination)
        compute_date = compute_contamination_for_date - start_date

        if concentration_point_contamination is None:
            input_concentration = diffuse_input_concentration 
        else:
            input_concentration = diffuse_input_concentration + concentration_point_contamination

        df_well_concentration = self.compute_concentration_in_well_at_date()

        # as fraction of the input concentration
        if as_fraction_input:
            df_well_concentration[:] = [x / input_concentration for x in df_well_concentration]
            ylabel = 'Fraction of input concentration'
        else: 
            ylabel = 'Concentration'
        
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
        plt.savefig('well_concentration_over_time_'+str(self.substance.substance_name)+'_'+schematisation_type+'.png', dpi=300, bbox_inches='tight')  # save_results_to + '/


    def compute_pathogen_removal(self):
        #AH_todo
        pass

    def plot_age_distribution(self):
        #AH_todo
        pass

    def plot_logremoval(self):
        #AH_todo
        pass
