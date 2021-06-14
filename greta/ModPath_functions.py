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

#%% ----------------------------------------------------------------------------
# INITIALISATION OF PYTHON e.g. packages, etc.
# ------------------------------------------------------------------------------

# %reset -f conda install #reset all variables for each run, -f 'forces' reset, !! 
# only seems to work in Python command window...

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
# from pandas import read_excel
from pandas import read_csv
from pandas import read_excel
import datetime
from tqdm import tqdm  # tqdm gives a progress bar for the simultation
# import pyarrow.parquet as pq
import math
import re # regular expressions
from scipy.special import kn as besselk

path = os.getcwd()  # path of working directory


#%% 

# ------------------------------------------------------------------------------
# Questions
# ------------------------------------------------------------------------------

# 1. 

# ------------------------------------------------------------------------------
# Phreatic and Semi-Confined Aquifer Functions
# ------------------------------------------------------------------------------
# %%

# ########### Defaults ###########
WELL_SCREEN_DIAMETER = .75  # m
BOREHOLE_DIAMETER = .75  # m -> equal to screen diameter to ignore backfilling
TEMPERATURE = 11  # Celcius
K_HOR_AQUIFER = 10  # m/d
VANI_AQUIFER = 1.  # -
K_HOR_CONFINING = .001  # m/d
VANI_CONFINING = 1.  # -
K_HOR_GRAVELPACK = 100  # m/d
VANI_GRAVELPACK = 1.  # -
K_HOR_CLAYSEAL = .001  # m/d
VANI_CLAYSEAL = 1.  # -

DENSITY_AQUIFER = 2650.  # kg/m3
REMOVAL_FUNCTION = 'omp'

DZ_WELL = 0.5  # preferred height of layer [m]
#%%
"""
Parameters
----------
schematisation: string
    'freatic', 'semi-confined', 'riverbankfiltration', 'basinfiltration'
removal_function: string
    'omp' -> get parameters for OMP
    'microbiology' -> get parameters for microbiology
"""


# # ########### INPUT PARAMETERS Aquapriori Bodem "Phreatic OMP" ###########
# schematisation = 'freatic'
# thickness_vadoze_zone = 1.  # m
# thickness_shallow_aquifer = 5.  # m
# thickness_target_aquifer = 10.  # m
# porosity_vadoze_zone = .2  # m3/m3
# porosity_shallow_aquifer = .3  # m3/m3
# porosity_target_aquifer = .25  # m3/m3
# organic_carbon_vadoze_zone = .2  # kg/m3 ??
# organic_carbon_shallow_aquifer = .3  # kg/m3 ??
# organic_carbon_target_aquifer = .25  # kg/m3 ??
# redox_vadoze_zone = 1.  # 1 = (sub)oxic; 2 = anoxic; 3 = deeply anoxic
# redox_shallow_aquifer = 2
# redox_target_aquifer = 3
# well_discharge_m3hour = 20 #m3/h
# recharge_rate = .001
# recharge_conc = 1.
# substance = 'chloridazon'
# vertical_resistance_aquitard   # [d], c_V
# soil_moisture_content           # [m3/m3], θ


# #@basin paramters
# length_basin
# width_basin
# _depth_basin
# horizontal_distance_basin_gallery = horizontal distance between basin bank and drainage gallery [m];
# porosity_recharge_basin 
# groundwater_level_above_saturated_zone = normal maximum rise of watertable above H0 [m];

# from greta.draft_transport_function import HydroChemicalSchematisation as HCS
# HCS_test = HCS()
# print(vars(HCS_test))

# #%%
# class AnalyticalWell():
#     """ Compute travel time distribution using analytical well functions."""
 
#   	def __init__(self):
#     		""" 'unpack/parse' all the variables from the hydrogeochemical schematizization """
#   	  	for key, item for input_dict.items():
 
  
#     def _check_init_freatic():
#        	#check the variables that we need for the individual aquifer types are not NONE aka set by the user
#   			pass
 
#   	def export_to_df(self, what_to_export='all')
#   	    """ Export to dataframe....

#         Parameters
#         ----------
#         what_to_export: String
#         		options: 'all', 'omp_parameters', 'microbial_parameters'
#         """
#   			#delete the unwanted columns depending on what the user asks for here
#   			returns df_flowline, df_particle

#%%  
# the python user will call the function as follows
# well = AnalyticalWell()
# if schematisation == 'freatic':
# 		well.freatic()
# elif schematisation == 'semiconfined':
# 		well.semiconfined()
# else:
#   	raise KeyError('schematisation argument not recognized')
# df_flow, df_particle = well.export_to_df('all')
#%%
# output Alex "phreatic_dict_nogravel.txt" --> saved in "testing dir"
research_dir = os.path.join(path,"..","research")
with open(os.path.join(research_dir,"phreatic_dict_nogravel.txt"),"r") as file_:
    dict_content = file_.read()
    phreatic_scheme = eval(dict_content)

for iKey,iVal in phreatic_scheme.items():
    print(iKey,iVal,"\n")

delr_bounds = []
ncol_bounds = []
for iDict_main in phreatic_scheme:
    
    if iDict_main in ["geo_parameters","recharge_parameters","ibound_parameters",
                      "well_parameters"]:
        
        for iDict_sub in phreatic_scheme[iDict_main]:
            if (iDict_main == "geo_parameters") & (iDict_sub == "vadose"):
                continue
            print(iDict_main,iDict_sub)
            try:
                delr_bounds.append(phreatic_scheme[iDict_main][iDict_sub]["rmin"])
            except Exception as e: print(e)
            try:
                delr_bounds.append(phreatic_scheme[iDict_main][iDict_sub]["rmax"])
            except Exception as e: print(e)
            try:
                ncol_bounds.append([phreatic_scheme[iDict_main][iDict_sub]["ncols"],
                                    phreatic_scheme[iDict_main][iDict_sub]["rmin"],
                                    phreatic_scheme[iDict_main][iDict_sub]["rmax"],
                                    iDict_main,iDict_sub])
            except Exception:
                ncol_bounds.append([1,
                                    phreatic_scheme[iDict_main][iDict_sub]["rmin"],
                                    phreatic_scheme[iDict_main][iDict_sub]["rmax"],
                                    iDict_main,iDict_sub])
print(delr_bounds,"\n",ncol_bounds)
np.unique(np.array(delr_bounds))
# Horizontal boudnaries determined by rmin, rmax
# 
# Refinement1 --> near well 

#%%
'''
phreatic_scheme = {'simulation_parameters': 
    {'schematisation_type': 'phreatic', 
    'computation_method': 'modpath',
    'temp_correction_Koc': True,
    'temp_correction_halflife': True,
    'biodegradation_sorbed_phase': True,
    'compute_thickness_vadose_zone': True,
    'start_date_well': datetime.date(1950, 1, 1),
    'start_date_contamination': datetime.date(1950, 1, 1),
    'compute_contamination_for_date': datetime.date(1950, 4, 11)},
'geo_parameters': 
    {'vadose': 
        {'vadose': True,
        'top': 22,
        'bot': 17,
        'rmin': 0.375,
        'rmax': 1723.5846804982755,
        'porosity': 0.38,
        'moisture_content': 0.15,
        'solid_density': 2.65,
        'f_oc': 0.001,
        'redox': 'anoxic',
        'DOC': 10,
        'pH': 5,
        'T': 11},
    'gravelpack1':
        {'top': 7,
        'bot': -33,
        'rmin': 0.2,
        'rmax': 1.0,
        'hk': 1000,
        'vani': 1},
    'clayseal1':
        {'top': 22,
        'bot': 7,
        'rmin': 0.2,
        'rmax': 2.0,
        'hk': 0.001,
        'vani': 1},
    'mesh_refinement1':
        {'rmin': 0.75,
        'rmax': 40,
        'ncols': 20},
    'mesh_refinement2': 
        {'rmin': 40,
        'rmax': 1723.5846804982755,
        'ncols': 30}
    },
'ibound_parameters':
    {'outer_boundary':
        {'head': 17,
        'top': 7,
        'bot': -33,
        'rmin': 1723.5846804982755,
        'rmax': 1724.5846804982755}
    }, 
'recharge_parameters':
    {'source1':
        {'substance_name': 'benzo(a)pyrene',
        'recharge': 0.0008213552361396304,
        'rmin': 0.75,
        'rmax': 1723.5846804982755,
        'DOC': 0.0,
        'TOC': 0.0,
        'c_in': 0}
    },
'well_parameters':
    {'well1':
        {'Q': 7665.599999999999,
        'top': 7,
        'bot': -33,
        'rmin': 0.0,
        'rmax': 0.2}
    },
'point_parameters':
    {'point1': 
        {'substance_name': 'benzo(a)pyrene',
        'c_in': 100.0,
        'r_start': 0,
        'z_start': 22,
        'q_point': 100.0}
    },
'substance_parameters':
    {'log_Koc': 6.43,
    'pKa': 99,
    'omp_half_life': 
        {'suboxic': 530,
        'anoxic': 2120,
        'deeply_anoxic': 2120}
    },
'bas_parameters': {}
}
'''
#%%
'''
Some clarifications in red below:

Algemeen-> What about temperature? Now it is capital T, I don’t see quickly in the modflow documentation what is wanted here

INTEGER vervangen door FLOAT (behalve nlayer, ncols en andere parameters die echt een integer zijn)
Bijv DOC, TOC, c_in, top, bot, etc.
Kan door aanpassen van default (1 -> 1.0)

Check hoofdlettergebruik. Belangrijk om dit consistent te doen om foutjes bij invoer of code te voorkomen.
top, bot -> prima
pH, DOC, TOC is prima om hoofdletters te doen want dat is een afkorting.
Recharge -> the flopy parameter is “rech” (https://flopy.readthedocs.io/en/latest/_modules/flopy/modflow/mfrch.html). Recharge is geen afkorting, dus moet sowieso met kleine letter.
Khor -> the flopy parameter is “hk” (ff overleg met martinK). Ik heb dit zelf niet goed in de tabel gezet.

Bas_parameters-> Ok then this is an empty dictionary, rest of params come in the Modpath class

-> rmax: mag weg, want wordt nu elders gedefinieerd.

Simulation_paramters 
-> Simulation_paramEters 

-> compute_contamination_for_date, start_date_contamination, end_date_well:
De value moet een “Date” zijn ipv integer. -> This is now a datetime type

Point_parameters
-> voeg even een voorbeeld toe voor Steven. Dat scheelt Steven tijd.

Geo_parameters
Vadoze, layer1, layer2: 
-> de Top en Bot moeten op elkaar aansluiten (mogen geen gaten tussen zitten) en in mASL. Volgens mij gaat hier iets mis.
-> rmin moet gelijk zijn aan de diameter van de diameter_borehole / 2 (0.75 default / 2) -> before this was the diameter_gravelpack, changed to diameter_borehole

Clayseal, gravelpack:
Voeg ajb even een apart voorbeeld toe voor Steven. Dat scheelt Steven tijd.
Dus 2 bestanden voor phreatic (met / zonder filterscreen) en 2 voor confined

Bestand 1:
Default waarden (dus zonder gravelpack)

Bestand 2:
Filterscreen met diameter 0.2 ,
Clayseal ter plaatse van layer1
Gravelpack ter plaatse van layer2

Meshrefinement1 
-> rmin: moet gelijk zijn aan straal van boorgat-> before this was the diameter_gravelpack, changed to diameter_borehole
-> rmax: moet gelijk zijn aan top – bottom van layer2

Meshrefiniment2 
-> rmin: moet gelijk zijn rmax van Meshrefinement1

Recharge_parameters
Rmin -> diameter_borehole / 2
Name -> vervangen door “substance_name”

Substance_parameters -> ok lets discuss, now these params are not passed to the dictionary yet unless user-specified, as the Substance class is not used until the Concentration class. 
Dit bij volgende overleg met martinK bespreken
Optie 1: “substance_name” als key toevoegen (makkelijk als we 1 stof per berekening doen)
Optie 2: nested dictionary van maken, met “substance_name” als key (dan kun je in 1x meerdere stoffen doen)

well_parameters
top -> top layer 1 -> do you mean the top of layer 2? I thought the well was in the target aquifer?
bot -> bottom layer 1 -> same as above, layer2?
rmin -> 0 (midden van put)
rmax -> diameter_filterscreen
'''


class ModPathWell:

    """ Compute travel time distribution using MODFLOW and MODPATH.""" 
    def __init__(self, schematisation: dict): #change schematisation_instance to schematisation
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
        # check the variables that we need for the individual aquifer types are not NONE aka set by the user
        '''check the variables that we need for the individual aquifer types are not NONE aka set by the user'''

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
    def _check_init_semi_confined():
        ''' check the variables that we need for the individual aquifer types are not NONE aka set by the user '''

    def assign_bas_parameters(self, ncols_filterscreen = 1,
                      ncols_gravelpack = 1,
                      ncols_near_well = 20,
                      ncols_far_well = 30,
                      diameter_filterscreen = 0.1,
                      diameter_gravelpack = 0.75,
                      thickness_aquifer = 20.,
                      model_radius = 500.):
        ''' Assign the grid discretization values '''

        # General input dictionary for discretisation:
        self.bas_parameters = {}
        # Number of columns from the well to the outer boundary
        self.bas_parameters["ncols_filterscreen"] = ncols_filterscreen
        self.bas_parameters["ncols_gravelpack"] = ncols_gravelpack
        self.bas_parameters["ncols_near_well"] = ncols_near_well
        self.bas_parameters["ncols_far_well"] = ncols_far_well
        # Filterscreen diameter
        self.bas_parameters["diameter_filterscreen"] = diameter_filterscreen
        # Gravelpack diameter
        self.bas_parameters["diameter_gravelpack"] = diameter_gravelpack
        
        # Thickness in phreatic aquifer is equal to top of upper layer minus
        # bottom of lowest layer () --> use regular expressions? 
        # self.thickness_aquifer = geo_parameters["layer1"]["top"] - geo_parameters["layer2"]["bot"]
        self.bas_parameters["thickness_aquifer"] = thickness_aquifer
        
        # Model radius is equal to the outer model boundary (defined earlier?)
        # self.model_radius = ibound_parameters["outer_boundary"]["rmax"]
        if not hasattr(self,"model_radius"):
            self.bas_parameters["rmax"] = model_radius
        # else:
        #     # Copy model radius from previously defined boundary
        #     self.bas_parameters["rmax"] = self.model_radius
        
def make_radial_discretisation(self, schematisation: dict, dict_keys = None,
                                res_hor_max = 10.):
    ''' Generate distance between columns for axisymmetric model.
    Sets it to self.delr
    - schematisation is of type dict with (sub)dictionaries with keys 'dict_keys'.
    - res_hor_max: maximum horizontal resolution if horizontal resolution 'res_hor'
      is not given in the subdictionary.

    # Assign column data
    - self.delr: column widths of the model columns [np.array]
      rounded to three decimals [mm scale].
    - self.xmid: x-coordinates (middle) of the model columns [np.array]
    - self.ncol: number of model columns

    # Assign row data : The rows have a predefined width [2 rows, 1 m width]
    - self.delc: row widths of the model rows [np.array]
    - self.ymid: y-coordinates (middle) of the model rows [np.array]
    - self.nrow: number of model rows [2]

    Local parameters:
    - refinement_bounds: boundaries where horizontal resolution values may change.
    - col_bounds: Represent the column edges (both xmin and xmax values) 
    '''
    # mesh refinement boundaries
    refinement_bounds = []
    # Detailed model column bounds
    col_bounds = []

    if dict_keys is None:
        dict_keys = [iDict for iDict in schematisation.keys()]
    
    for iDict in dict_keys:
        for iDict_sub in schematisation[iDict]:
            try:
                # rmin
                rmin = schematisation[iDict][iDict_sub]["rmin"]
                refinement_bounds.append(rmin)
            except Exception as e:
                print(e,"continue")
            try:
                # rmax
                rmax = schematisation[iDict][iDict_sub]["rmax"]
                refinement_bounds.append(rmax)
            except Exception as e:
                print(e,"continue")
            # horizontal resolution
            try:
                res_hor = schematisation[iDict][iDict_sub]["res_hor"]
            except Exception as e:
                # Default maximum resolution if "res_hor" is not in subdictionary
                res_hor = res_hor_max 

            num_cols = max(1,math.ceil((rmax-rmin)/res_hor))
            # Determine (in-between) column boundaries
            column_bounds = np.linspace(rmin,rmax,
                                        num = num_cols + 1, endpoint = True)
            col_bounds.extend(list(column_bounds))

    # Only keep unique values for refinement_bounds and col_bounds
    refinement_bounds = np.sort(np.unique(np.array(refinement_bounds)))   
    col_bounds = np.sort(np.unique(np.array(col_bounds)))


    # Create delr array along columns (rounded to 3 decimals)
    self.delr = np.round(np.diff(col_bounds),3)
    # Number of model columns
    self.ncol = len(self.delr)
    # Create xmid array from delr array
    self.xmid = np.empty((self.ncol), dtype= 'float')
    self.xmid[0] = self.delr[0] * 0.5
    for iCol in range(1, self.ncol):
        self.xmid[iCol] = (self.xmid[(iCol - 1)] + \
                            ((self.delr[(iCol)]) + (self.delr[(iCol - 1)])) * 0.5)  

   # Create delc array along rows (rounded to 3 decimals)
    self.delc = np.ones(2,dtype = 'float')
    # Number of model rows
    self.nrow = len(self.delc)

    # Create ymid array from delr array
    self.ymid = np.empty((self.nrow), dtype= 'float')
    self.ymid[0] = self.delc[0] * 0.5
    for iRow in range(1, self.nrow):
        self.ymid[iRow] = (self.ymid[(iRow - 1)] + \
                            ((self.delc[(iRow)]) + (self.delc[(iRow - 1)])) * 0.5)  

    '''
    def make_vertical_discretisation(self):
        """ Generate top's and bot's of MODFLOW layers.
        Sets it to self.tops, self.bots
        """
        self.top = 0.
        lay_bots = {}
        for iGeo in self.geo_parameters:
            lay_bots[iGeo] = self.geo_parameters[iGeo]["bot"]
  			self.bots = None
    
  	def make_discretisation():
        self.grid
 
  	def assign_material(self):  
        self.material  # name of dictionary
        self.wells
   
  
  	def assign_aquifer_parameters():
        self.khor
  		self.VANI
  		#etc.,

  	def assign_wells():
  	    self.WELL_pacakage  # np.array
  		self.Q
  		self.well_id

    	def assign_multi_nodal_well():
  	    self.wellid  # np.array
  		self.....
	def assign_ibound():
        self.ibound  # 3 D grid
  			self.head
    '''

    def assign_ibound(self,filter_top, filter_botm, filter_left, filter_right,zmid, xmid, ibound):

        ''' This fnction is used to assign the constant head cells (ibound --> "-1". 
            These represent a constant well abstraction in combination
            with a net precipitation influx. '''
        
        for iLay in range(self.nlay):
            for iCol in [0]: #range(100):
                if (zmid[iLay] < filter_top) & (zmid[iLay] > filter_botm):
    #                (xmid[iCol] > filter_left) & (xmid[iCol] < filter_right):
                    # Update ibound
                    ibound[iLay, 0, iCol] = -1
                    
        # return ibound
    '''
  	def assign_recharge_boundary():
  		self.recharge  #2 D grid
  
  
    def create_modflow_packages(... , *kwargs):
		""" Generate flopy files."""
		self.dry_cells  #
  		self.other_general_parameters  #
  
    def generate_modflow_files():
		""" Generate input files for modflow computation."""

    def run_modflow():
  		run modflow
  		if condition X is met:
			self.modflow_has_run = True
 
  	def modpath_startpoints_recharge():
				""" Generate startpoints for MODPATH forward computation.
        
        Parameters
        ----------
        xmin_start_particle  # close to well, default = 0
  			xmax_start_particle  # far away from well, default = model radius
  			number_of_particles  # user defined, defualt 100        
        
        Return
        ------
        df_startpoints: pd.DataFrame  # OR DICTIONARY IF THAT IS EASIER -> STEVEN
		        columns: col, row, lay, localx, localy, localz, flowline_id, startpoint_id, discharge
        """
  			init: self.xmin_start_particle = 0
  			self.number_of_particles = ....
	  		self.xmax_start_particle = well_radius if not defined
  
  			# distribute start point particles based on equal distance
				dx = (xmax - xmin) / number_of_particles
	
  			discharge = integrate volume
  
  			self.df_startpoints
  			
  
   	def modpath_endpoints_from_well():
  			# same as above, for wells
 
    def create_modpath_packages(self.df_startpoints):
				self.track_direction  # forward or backward
 
    def generate_modpath_files():
  
    def run_modpath():
 				self.modpath_has_run
  			if condition X is met:
						self.modpath_has_run = True
    '''
    def export_to_df(self, grid_material, what_to_export='all'):
  	    """ Export to dataframe....

        Parameters
        ----------
        what_to_export: String
        		options: 'all', 'omp_parameters', 'microbial_parameters'
        """
        # df_flowline = pd.DataFrame()
        # df_particle = pd.DataFrame()

  		#delete the unwanted columns depending on what the user asks for here
  		# return df_flowline, df_particle
    def phreatic(self):
        self._check_init_phreatic()

        # Make radial discretisation
        self.make_radial_discretisation()

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

'''
#%%  
class ModPathWell_OLD():
    """ Compute travel time distribution using MODFLOW and MODPATH.""" 

    def _check_init_freatic():
        # check the variables that we need for the individual aquifer types are not NONE aka set by the user
  
    def _check_init_semi_confined():
        # check the variables that we need for the individual aquifer types are not NONE aka set by the user

	  def model_size_freatic():
				""" Compute Radius based on discharge of wells and recharge rate.
        
        Returns
        -------
        Radius of model
        """

 		def make_hydrogeological_schematisation()
  			""" make hydrogeological schematisation.
            Sets it to self.radius
        """
        #  convert dict with geological schematisation at well1 location to df 
        df = pd.DataFrame([self.geodict])
  		  # order df from top to bot
  
  			# to do later: give warning if space between top-bottom

  			# identify aquifers and confining layers
  			mask = df['k_vertical'] < 0.01
  
  			# add column to df -> contains well True/False
  
        # find all layers containing well screen
			  # target_layers --> (top 'well' > bottom 'formation') & (bottom 'well' < top 'formation')
  
  			# get aquifer properties (kD, c) in layer containing well -> kD, c_below, c_above
  
  			self.df_schematisation = df
  
    def model_size_semiconfined():
				""" Compute radius based on 5 * spreading distance (labda = sqrt[kDc]).
            Sets it to self.radius
        """
  			KD = .... functie -> 
  			c = .... functie ->  1/ (1/c_above + 1/c_below)
  			labda = sqrt(KD * c)
  			radius = 5 * labda
  			self.radius = radius

    def discretization(radius, df_schematisation):
  			""" Generate discretisation of model.
        
        Parameters
        ----------
        borehole -> 0, 1 or many columns 
        """"
	  	  # Make extra second row to allow modpath
  
  			# horizontal: filterscreen
  
  		  # horizontal: borehole
  
  			# horizontal: tot 1xlengte  filters (of target aquifer)
  
  			# horizontal: tot modelrand
  
  			# verticaal: lagen met putten -> user defined dZ binnen putfilter
  
  			# verticaal: lagen zonder putten -> user defined dZ
  
  			# verticaal: phreatic aquifer -> prevent that cells fall dry
  
  			return del_r, del_c, del_l  #also x_mid, z_mid??

  	def assign_array()
  		'""" for each grid cell -> geo1, geo2, well1, well2, etc."""
  	    self.grid_material  # np.array
 
  	def fixed_head_well_boundary():
        self.grid_ibound

	  def fixed_discharge_well_boundary():
        self.grid_fixed_discharge

  	def recharge_boundary():
  			self.grid2D_recharge
  
  	def fixed_head_model_boundary():

  
	  def no_flow_model_boundary():
  
  
  	def assign_parameters():
  			""" Convert grid_material to parameters."""
  	    self.K_vertical_grid =function(grid_material, dictionary_material_properties)

    def assign_parameters_axisysmetric():

    def create_modflow_packages(... , *kwargs):
				""" Generate flopy files."""
				self.dry_cells  #
  			self.other_general_parameters  #
  
    def generate_modflow_files():
				""" Generate input files for modflow computation."""

    def run_modflow():
  			run modflow
  			if condition X is met:
						self.modflow_has_run = True
 
  	def modpath_startpoints_from_recharge():
				""" Generate startpoints for MODPATH computation.
        
        Parameters
        ----------
        xmin_start_particle  # close to well, default = 0
  			xmax_start_particle  # far away from well, default = model radius
  			number_of_particles  # user defined, defualt 100        
        
        Return
        ------
        df_startpoints: pd.DataFrame  # OR DICTIONARY IF THAT IS EASIER -> STEVEN
		        columns: col, row, lay, localx, localy, localz, flowline_id, startpoint_id, discharge
        """
  			init: self.xmin_start_particle = 0
  			self.number_of_particles = ....
	  		self.xmax_start_particle = well_radius if not defined
  
  			# distribute start point particles based on equal distance
				dx = (xmax - xmin) / number_of_particles
	
  			discharge = integrate volume
  
  			self.df_startpoints
  			
  
   	def modpath_endpoints_from_well():
  			# same as above, for wells
 
    def create_modpath_packages(self.df_startpoints):
				self.track_direction  # forward or backward
 
    def generate_modpath_files():
  
    def run_modpath():
 				self.modpath_has_run
  			if condition X is met:
						self.modpath_has_run = True

  	def export_to_df(self, grid_material, what_to_export='all')
  	    """ Export to dataframe....

        Parameters
        ----------
        what_to_export: String
        		options: 'all', 'omp_parameters', 'microbial_parameters'
        """
  			#delete the unwanted columns depending on what the user asks for here
  			returns df_flowline, df_particle

well = ModpathWell()
if tracking_direction = 'forward'
		well.modpath_startpoints_from_recharge()
elif tracking_direction = 'backward'
		well.modpath_endpoints_from_well()
else
  	raise KeyError('tracking direction argument not recognized')
df_flow, df_particle = well.export_to_df('all')
'''

  
#%%

'''
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
        Column 'start_or_end_point_id': Integer ######################################
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

    def __init__(self, substance: Substance, df_particle, df_flowline, removel_function?):
        self.omp_inialized = False

  
  	def _init_omp()
  		if self.omp_inialized:
	  		self.df_part['Kow'] = self.df_part['redox_zone'].apply(lambda x: substance.get_Kow(x)
   		self.omp_inialized = True


  	def _init_microbiology()


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
       self.df_part...

    def compute_microbiology_removal(self):
                                                             
	  def compute_well_concentration(self, evaluation_time = None)
        """ Returns the concentration in the raw water of each well (as defined by endpoind_id)."""
      	if evaluation_time is None:
          select all flowlines
        else:
          select flowline with break_through_time < evaluation_time                                                  
		  	conc_flowline = concentration at end of selected flowline
				concentration_well = sum (conc_selected_flowline_i * discharge_flowline_i) / sum discharge_all_flowline                                                             


    def plot_concentration(self)
                                           
    def plot_age_distribution(self)

    def plot_logremoval(self)

'''

'''
#%%
class Test():
	def __init__(self):
		self.alex = None
                                                             
	def call_alex(self):                                    
		self.alex = 'called'
	def print_it(self):
		print(self.alex)
                                                    
# %%                       
# the python user will call the function as follows
concentration = Concentration()
if removal_function == 'omp':
		concentration.compute_omp_removal
elif removal_function = 'omp':
		concentration.compute_microbiology_removal
else:
  	raise KeyError('schematisation argument not recognized')
                       
'''