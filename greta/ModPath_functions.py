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
# functions to deal with numpy rec.array data
import numpy.lib.recfunctions as rfn  
import pandas as pd
import os
import sys
import copy

# from pandas import read_excel
from pandas import read_csv
from pandas import read_excel
import datetime
# import pyarrow.parquet as pq
import math
import re # regular expressions
from scipy.special import kn as besselk

try:
    from greta.Analytical_Well import * 
    from greta.Substance_Transport import *
except ModuleNotFoundError as e:
    print(e, ": second try.")
    module_path = os.path.join("..","greta")
    if module_path not in sys.path:
        sys.path.insert(0,module_path)
    from Analytical_Well import * 
    from Substance_Transport import *

    print("Second try to import modules succeeded.")

# path of working directory
path = os.getcwd()  

# run installed version of flopy or add local path 
# (add flopy to requirements.txt: pip install flopy==3.3.1)
try:
    import flopy
    import flopy.utils.binaryfile as bf
except Exception as e:
    flopypth = os.path.abspath(os.path.join('..', '..'))
    sys.path.append(flopypth)
    import flopy
    import flopy.utils.binaryfile as bf

# flopy version
print(flopy.__version__)


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
    def __init__(self, schematisation: object or dict,
                       workspace: str or None = None, modelname: str or None = None,
                       bound_left: str = "xmin", bound_right: str = "xmax",
                       bound_top: str = "top", bound_bot: str = "bot",
                       bound_north: str = "ymin", bound_south: str = "ymax"): 
        ''''unpack/parse' all the variables from the hydrogeochemical schematizization """
       
        Parameters
        ----------
        df_flowline: pandas.DataFrame
            Column 'flowline_id': Integer
            Column 'discharge': Float
                Discharge associated with the flowline (m3/d)
            Column 'particle_release_day': Float
            Column 'input_concentration'
            Column 'endpoint_id': Integer
                ID of Well (or drain) where the flowline ends

        df_particle: pandas.DataFrame
            Column 'flowline_id'
            Column 'travel_time'
            Column 'xcoord'
            Column 'ycoord'
            Column 'zcoord'
            Column 'redox'
            Column 'temperature'
        '''

        # get the non-default parameters
        # self.test_variable = None #AH test variable here to see if errors are caught....
        
        self.workspace = workspace  # workspace
        self.modelname = modelname  # modelname
        # Cbc unit flag (modflow oc-package input)
        iu_cbc = 130
        self.iu_cbc = iu_cbc
        # Expected boundary terms in dicts
        self.bound_left = bound_left
        self.bound_right = bound_right
        self.bound_top = bound_top
        self.bound_bot = bound_bot
        self.bound_north = bound_north
        self.bound_south = bound_south

        # Create output directories
        # Destination root
        self.dstroot = self.workspace + '\\results'  # + '\\results_' + 'yr' + str(yr) + '_sp' + str(per)
        self.dstdircbc = os.path.join(self.dstroot, 'cbc')
        self.dstdirhds = os.path.join(self.dstroot, 'hds')

        # create all directories
        if not os.path.exists(self.dstdirhds):
            os.makedirs(self.dstdirhds)
        if not os.path.exists(self.dstdircbc):
            os.makedirs(self.dstdircbc)

        # Source files
        self.model_hds = os.path.join(self.workspace, self.modelname + '.hds')
        self.model_cbc = os.path.join(self.workspace, self.modelname + '.cbc')
        '''
        Initialize the AnalyticalWell object by making the dictionaries and adding the schematisation (and therefore attributes 
        of the schematisation) to the object.

        Returs
        ------
        Dictionaries (Modflow) with all parameters as an attribute of the function.

        '''
        self.schematisation = schematisation

        if type(self.schematisation) == "dict":
            self.schematisation_dict = self.schematisation
        else:
            #Make dictionaries
            # self.schematisation.make_dictionary()
            self.schematisation_dict = {}
            # Required keys
            self.required_keys = ["simulation_parameters","geo_parameters",
            "ibound_parameters","recharge_parameters",
            "well_parameters","point_parameters", "mesh_refinement",
            "endpoint_id"]
            self._create_schematisation_dict(self.required_keys)
        
    def _create_schematisation_dict(self,required_keys):
        '''
        # Create schematisation_dict 
        # required_keys = ["simulation_parameters","geo_parameters",
        # "ibound_parameters","recharge_parameters",
        # "well_parameters","point_parameters"]
        '''
        for iKey in required_keys:
            if hasattr(self.schematisation,iKey):
                self.schematisation_dict[iKey] = getattr(self.schematisation,iKey)

    def _check_schematisation(self,required_keys):
        
        for req_key in required_keys:
            if not hasattr(self.schematisation,req_key):
                print(f'Error, required variable {req_key} is not defined in schematisation dict.')

    def _check_required_variables(self,required_variables):
        for req_var in required_variables:
            value = getattr(self, req_var)
            print(req_var,value)
            if value is None:
                raise KeyError(f'Error, required variable {req_var} is not defined.')


    def _check_init_phreatic(self):
        # check the variables that we need for the individual aquifer types are not NONE aka set by the user
        '''check the variables that we need for the individual aquifer types are not NONE aka set by the user'''
        
        
        # # Required keys in self.schematisation
        # required_keys = ["simulation_parameters",#repeat for all
        #                       "geo_parameters",
        #                       "recharge_parameters",
        #                       "ibound_parameters",
        #                       "well_parameters",
        #                       "point_parameters",
        #                       "substance_parameters"
        #                     ]

        # Required variables (to run model)
        required_variables = []

        # Check schematisation dictionary
        self._check_schematisation(self.required_keys)
        # Check required variables
        self._check_required_variables(required_variables)
        # # Fill schematisation dictionary
        # self._create_schematisation_dict(self.required_keys)

    def _check_init_semi_confined():
        ''' check the variables that we need for the individual aquifer types are not NONE aka set by the user '''

    def update_property(self, property, value):
        setattr(self, property, value)
        ''' set attribute value to new object attribute/property'''

    def axisym_correction(self, grid: np.array,
                          dtype: str or None = None, theta = 2 * np.pi):

        ''' Correct modflow parameters to correctly calculate axisymmetric flow.
            Adjust the array 'grid' if it already exists, using multiplier theta
            and xmid (nrow == 1 or nrow == 2) or ymid (ncol == 1 or ncol == 2).
        '''
        # Check if model is axisymmetric along rows or columns.
        if dtype is None:
            dtype = grid.dtype
        # Create empty numpy grid    
        grid_axi = np.zeros((self.nlay,self.nrow,self.ncol), dtype = dtype)
        if (self.nrow == 1) | (self.nrow == 2):
            for iCol in range(self.ncol):
                grid_axi[:,:,iCol] = theta * self.xmid[iCol] * grid[:,:,iCol]

        return grid_axi

    def cell_bounds(self,schematisation: dict, dict_key: str = "None",
                        dict_subkey: str = "None",
                        model_type = "axisymmetric"):
        ''' Create cell boundaries of box created using following boundaries:
            self.bound_left: str = "xmin" 
            self.bound_right: str = "xmax"
            self.bound_top: str = "top"
            self.bound_bot: str = "bot"
            self.bound_north: str = "ymin"  # not required for axisymmetric or 2D model
            self.bound_south: str = "ymax"  # not required for axisymmetric or 2D model  
            # Return boundary indices of row and columns plus parameter value
            return layidx_min,layidx_max,rowidx_min,rowidx_max,colidx_min,colidx_max                
        '''
        
        # Loop through schematisation keys (dict_keys)
        iDict = dict_key
        # Use subkeys of schematisation dictionary
        iDict_sub = dict_subkey   

        # coordinate values of boundaries [float]
        left = schematisation[iDict][iDict_sub][self.bound_left]
        right = schematisation[iDict][iDict_sub][self.bound_right]
        try:
            top = schematisation[iDict][iDict_sub][self.bound_top]
        except KeyError:
            print("set top of", iDict, iDict_sub, "to 0.")
            top = 0.
        try: 
            bot = schematisation[iDict][iDict_sub][self.bound_bot]
        except KeyError:
            print("set bottom of", iDict, iDict_sub, "to model bottom.")
            bot = min(self.zmid)

        if not self.model_type == "axisymmetric":
            try:
                north = schematisation[iDict][iDict_sub][self.bound_north]
                south = schematisation[iDict][iDict_sub][self.bound_south]
                # Determine row indices
                rowidx_min = int(np.argwhere((self.ymid <= north) & (self.ymid >= south))[0])
                rowidx_max = int(np.argwhere((self.ymid <= north) & (self.ymid >= south))[-1]) + 1
            except KeyError as e:
                print(e,f"missing {iDict} {iDict_sub}. Continue")
        else:
            north,south,rowidx_min,rowidx_max = None, None,0,1
        try:
            # Determine layer indices
            layidx_min = int(np.argwhere((self.zmid <= top) & (self.zmid >= bot))[0])
            layidx_max = int(np.argwhere((self.zmid <= top) & (self.zmid >= bot))[-1]) + 1
        
        except IndexError as e:
            print(e, iDict,iDict_sub,top,bot, "(top,bot)")
            print("Set layidx_min and layidx_max to None.")
            layidx_min, layidx_max = None, None

        try:
            # Determine column indices
            colidx_min = int(np.argwhere((self.xmid >= left) & (self.xmid <= right))[0])
            colidx_max = int(np.argwhere((self.xmid >= left) & (self.xmid <= right))[-1]) + 1
            # np.where((self.xmid < right) & (self.xmid > left))    
        except IndexError as e:
            print(e, iDict,iDict_sub,left,right, "(left,right)")
            print("Set colidx_min and colidx_max to None.")
            colidx_min, colidx_max = None, None

        return layidx_min,layidx_max,rowidx_min,rowidx_max,colidx_min,colidx_max
            

    def fill_grid(self,schematisation: dict, dict_keys: list or None = None,
                        parameter: str = "None",
                        grid: np.array or None = None,
                        dtype: str or None = 'float'):
        ''' Assign values to 'grid' [np.array] for parameter name 'parameter' [str], 
            using the keys dict_keys [list] in schematisation dictionary 'self.schematisation'.
            'dtype' [grid dtype] --> grid dtype [str] is obtained from grid if initial array is given.
            The schematisation type refers to the model_type["axisymmetric","2D" or "3D"]
            
            Boundaries of the parameter values are to be included in the dictionaries:

            self.bound_left: str = "xmin" 
            self.bound_right: str = "xmax"
            self.bound_top: str = "top"
            self.bound_bot: str = "bot"
            self.bound_north: str = "ymin"
            self.bound_south: str = "ymax"

            The function returns:
            - grid  # grid [np.array] filled with (numeric) values for parameter 'parameter'.

        ''' 
        
        if grid is None:
            # Set dtype
            if dtype is not None:
                dtype = dtype
            else:
                dtype = "float"
            # Create empty np.array
            grid = np.ones((self.nlay,self.nrow,self.ncol), dtype = dtype)
        else:  # np.array is given
            if dtype is not None:
                dtype = grid.dtype
            else:
                dtype = "float"

        if dict_keys is None:
            dict_keys = [iDict for iDict in self.schematisation_dict.keys()]
       
        
        # Loop through schematisation keys (dict_keys)
        for iDict in dict_keys:
            # Loop through subkeys of schematisation dictionary
            for iDict_sub in schematisation[iDict]:   
                
                if parameter in schematisation[iDict][iDict_sub]:



                    # try:
                    # Obtain parameter value
                    parm_val = schematisation[iDict][iDict_sub][parameter]
                    # Obtain cell boundary limits
                    layidx_min,layidx_max,\
                        rowidx_min,rowidx_max,\
                        colidx_min,colidx_max = self.cell_bounds(schematisation,
                                                dict_key = iDict,
                                                dict_subkey = iDict_sub,
                                                model_type = self.model_type)
                    # Fill grid with parameter value 'parm_val'
                    if not None in [layidx_min,layidx_max,colidx_min,colidx_max]:
                        grid[layidx_min: layidx_max,\
                            rowidx_min: rowidx_max,\
                            colidx_min: colidx_max] = parm_val
                    if self.model_type in ["axisymmetric","2D"]:
                        # In 2D model or axisymmetric models an 
                        # inactive row is added to be able to run Modpath successfully.
                        grid[:,1,:] = grid[:,0,:]
  
        # Return the filled grid                
        return grid

    def _assign_cellboundaries(self,schematisation: dict, dict_keys: list or None = None,
                                bound_min: str = "xmin", bound_max: str = "xmax",
                                n_refinement: str = "ncols", ascending: bool = True,
                                res_max: int or float or None = None):
        ''' Function to help creating the grid discretisation inside 'make_discretisation'.
            Determines grid refinement in either the X (ncol), Y (nrow), or Z (nlay) direction.

            Return terms:
            - len_arr        # Length of grid in X,Y or Z-direction [int]
            - cell_sizes     # Distance between cell boundaries in X,Y, or Z-direction [1D np.array]
            - center_points  # Center locations of cells in the direction of grid refinement [1D np.array]
            - bound_list     # Model coordinates of cell boundaries in X, Y or Z-direction [1D np.array]

        '''

        # Keep track of grid boundaries
        bound_list = []

        if dict_keys is None:
            dict_keys = [iDict for iDict in self.schematisation_dict.keys()]

        # Loop through schematisation keys (dict_keys)
        for iDict in dict_keys:
            # Loop through subkeys of schematisation dictionary
            for iDict_sub in schematisation[iDict]:   
                try:
                    # minimum bound
                    val_min = schematisation[iDict][iDict_sub][bound_min]
                except KeyError as e:
                    print(e,f"missing {iDict} {iDict_sub}. Continue")
                    continue

                try:
                    # maximum bound
                    val_max = schematisation[iDict][iDict_sub][bound_max]
                except KeyError as e:
                    print(e,f"missing {iDict} {iDict_sub}. Continue")
                    continue

                try:
                    # number of refinements
                    n_ref = schematisation[iDict][iDict_sub][n_refinement]
                except KeyError:
                    if res_max is None:
                        n_ref = 1
                    else:
                        # Limit the cell resolution using 'res_max' (if not None)
                        n_ref = max(1,math.ceil((val_max-val_min)/res_max))
                    pass  

                # Calculate local resolution [L]
                resolution = (val_max-val_min) / n_ref   

                # Determine (in-between) column boundaries
                boundaries = np.linspace(val_min,val_max,
                                            num = n_ref + 1, endpoint = True)
                bound_list.extend(list(boundaries))

        # Only keep unique values for boundary list 'bound_list' 
        if ascending:  
            bound_list = np.sort(np.unique(np.round(bound_list,3)))
        else:
            bound_list = np.sort(np.unique(np.round(bound_list,3)))[::-1]
        # size of grid cells in the dimension (delr, delc or delv) from zeroeth to n_th cell
        cell_sizes =  abs(np.diff(bound_list))

        # Length of array
        len_arr = len(cell_sizes)

        # Assign center points (xmid | ymid | zmid)
        center_points = np.empty((len_arr), dtype= 'float')
        if ascending: # xmid and ymid arrays are increasing with increasing index number
            center_points[0] = bound_list[0] + cell_sizes[0] * 0.5
            for idx in range(1, len_arr):
                center_points[idx] = center_points[(idx - 1)] + ((cell_sizes[(idx)] + cell_sizes[(idx - 1)]) * 0.5)  
        else: # zmid arrays are decreasing with increasing index number
            center_points[0] = bound_list[0] - cell_sizes[0] * 0.5
            for idx in range(1, len_arr):
                center_points[idx] = center_points[(idx - 1)] - ((cell_sizes[(idx)] + cell_sizes[(idx - 1)]) * 0.5)  

        return len_arr, cell_sizes, center_points, bound_list

        
    def make_discretisation(self, schematisation: dict, dict_keys = None,
                            model_type = 'axisymmetric'):
        ''' Generate spatial grid for model_type choices: 'axisymmetric', '2D' or '3D'.
            
        Parameter 'schematisation' is of type dict with (sub)dictionaries with keys 'dict_keys'.
        The subdictionaries should contain specific keyword arguments to generate the grids.  
        The function indirectly uses "_assign_cellboundaries" to obtain grid data outputs.
        
        # Layer data assignment (model_type: axi-symmetric | 2D | 3D)
        Required keys:
        - "bot"  # bottom of local grid refinement
        - "top"  # top of local grid refinement
        Optional keyword argument(s):
        - "nlayers"  # number of layers within local grid refinement

        - self.delv: layer depths of the model layers [np.array]
                     rounded to two decimals [cm scale].
        - self.zmid: z-coordinates (middle) of the model layers [np.array]
        - self.nlay: number of model layers
        - self.top: model top [float or np.array]
        - self.bot: model bottoms per layer [1D | 2D | 3D np.array of floats]
                
        # Column data assignment (model_type: axi-symmetric | 2D | 3D) 
        Required keys:
        - "xmin"  # left side of local grid refinement
        - "xmax"  # right side of local grid refinement
        Optional keyword argument(s):
        - "ncols" # number of columns within local grid refinement

        # Column data outputs:
        - self.delr: column widths of the model columns [np.array]
            rounded to three decimals [mm scale].
        - self.xmid: x-coordinates (middle) of the model columns [np.array]
        - self.ncol: number of model columns


        # Row data assignment (model_type: 3D)
        # N.B. for modeltype: (axisymmetric | 2D) the rows have a predefined width [2 rows, 1 m width]

        Required keys:
        - "ymin"  # 'northern' side of local grid refinement
        - "ymax"  # 'southern' side of local grid refinement
        Optional keyword argument(s):
        - "nrows" # number of rows within local grid refinement

        # Row data outputs:
        - self.delc: row widths of the model rows [np.array]
        - self.ymid: y-coordinates (middle) of the model rows [np.array]
        - self.nrow: number of model rows
        
        Optional local outputs to return:
        - empty_grid  # Empty numpy array with size (self.nlay,nrow,ncol)
        - lay_bounds  # Upper and lower boundaries of the model grid cells [1D-array]
        - row_bounds  # Left and right boundaries of the model grid cells [1D-array]
        - col_bounds  # North-south boundaries of the grid-cells [1D-array]
        '''

        if schematisation is None:
            schematisation = getattr(self,"schematisation_dict")

        # Assign delv and zmid   
        self.nlay, self.delv, self.zmid, lay_bounds = self._assign_cellboundaries(schematisation = schematisation,
                                                                                  dict_keys = dict_keys,
                                                                                  bound_min = self.bound_bot, bound_max = self.bound_top,
                                                                                  n_refinement = "nlayers", ascending = False)                                   
 
        # Model top
        self.top = max(lay_bounds)
        # Model bottoms
        self.bot = self.top - self.delv.cumsum() 
        # Assign delr and xmid
        self.ncol, self.delr, self.xmid, col_bounds = self._assign_cellboundaries(schematisation = schematisation,
                                                                                  dict_keys = dict_keys,
                                                                    bound_min = self.bound_left, bound_max = self.bound_right,
                                                                    n_refinement = "ncols", ascending = True)

        # Assign delc and ymid
        if model_type in ["axisymmetric","2D"]:
            self.nrow = 2
            self.delc = np.ones((self.nrow),dtype = 'float')
            self.ymid = np.array([0.5,1.5])
            row_bounds = np.array([0.,1.,2.])
        else:
            self.nrow,self.delc,self.ymid,row_bounds = self._assign_cellboundaries(schematisation = schematisation,
                                                                                  dict_keys = dict_keys,
                                                                    bound_min = self.bound_north, bound_max = self.bound_south,
                                                                    n_refinement = "nrows", ascending = True)

        # Create empty model grid
        empty_grid = np.empty((self.nlay,self.nrow,self.ncol), dtype = 'float')

        return empty_grid, lay_bounds, row_bounds, col_bounds

    def set_ibound(self,schematisation: dict,
                        dict_keys: dict or None = None,
                        ibound: np.array or None = None,
                        strt: np.array or None = None):

        ''' This function is used to assign the constant head cells (ibound --> "-1". 
        '''
        if ibound is None:
            print("attribute ibound does not yet exist. Create empty ibound grid")
            ibound = np.ones((self.nlay,self.nrow,self.ncol), dtype = 'int')
                
        if strt is None:
            print("attribute strt does not yet exist. Create empty ibound grid")
            strt = np.zeros((self.nlay,self.nrow,self.ncol), dtype = 'float')
                
        
        
        
        # Relevant dictionary keys
        self.ibound = self.fill_grid(schematisation = self.schematisation_dict,
                                    dict_keys = dict_keys,
                                    parameter = "ibound",
                                    grid = ibound,
                                    dtype = 'int')
        self.strt = self.fill_grid(schematisation = self.schematisation_dict,
                                    dict_keys = dict_keys,
                                    parameter = "head",
                                    grid = strt,
                                    dtype = 'int')
               
    def oc_input(self, spd_oc = {(0, 0): ['save head', 'save budget']}):
        ''' Load OC package parms to model. '''
        self.spd_oc = spd_oc

    def _validate_input(self):
        # controleren op alle benodigde keys van de 
        # schematization class
        # if 'xmin' is not in dictionary.keys():
        #     raise ValueError('xmin is not in dictionayr')
        # dit soort validaties moet eigenlij in de schematization class
        pass

    def assign_material(self,schematisation: dict,
                        dict_keys: dict or None = None):
        ''' Assign grid material using subkeys in schematisation_dict.'''

        if schematisation is None:
            schematisation = getattr(self,"schematisation_dict")

        # dtype of material grid
        dtype = 'object'  # to assign a string without length limitation
        if dict_keys is None:
            dict_keys = [iDict for iDict in self.schematisation_dict.keys()]

        self.material = np.empty((self.nlay,self.nrow,self.ncol), dtype = dtype)
        # Loop through schematisation keys (dict_keys)
        for iDict in dict_keys:
            # Loop through subkeys of schematisation dictionary
            for iDict_sub in schematisation[iDict]:   

                if 'mesh_refinement' not in iDict_sub:

                    # Grid indices
                    try:
                        layidx_min,layidx_max,\
                            rowidx_min,rowidx_max,\
                            colidx_min,colidx_max = self.cell_bounds(schematisation,
                                                    dict_key = iDict,
                                                    dict_subkey = iDict_sub,
                                                    model_type = self.model_type)
                    except Exception as e:
                        print(e,"Continue.")
                        continue
                    # Fill grid with parameter value 'parm_val'
                    if not None in [layidx_min,layidx_max,colidx_min,colidx_max]:

                        self.material[layidx_min: layidx_max,\
                            rowidx_min: rowidx_max,\
                            colidx_min: colidx_max] = iDict_sub
                    
        if self.model_type in ["axisymmetric","2D"]:
                # In 2D model or axisymmetric models an 
                # inactive row is added to be able to run Modpath successfully.
                self.material[:,1,:] = self.material[:,0,:]

    def assign_wellloc(self,schematisation: dict,
                        dict_key: str = "well_parameters",
                        well_names: list or str or None = "None",
                        discharge_parameter = "Q"):

        ''' Determine the location of the pumping wells and the relative discharge per cell.
            
            Boundaries of the well locations are to be included in the dictionaries with 
            keys "dict_keys" [list], based on presence of "discharge_parameter" [str (default = "Q")]:
                       
            bound_left: str = "xmin" 
            bound_right: str = "xmax"
            bound_top: str = "top"
            bound_bot: str = "bot"
            bound_north: str = "ymin"
            bound_south: str = "ymax"

            # Output:
            # Creates a dict with list of tuples representing the well locations
            well_loc = {"well1": [(5,0,0),(6,0,0)],
                        "well2":}
            
            # Calculate cumulative KD per well [m2/day] (to calculate relative discharge per cell)
            KD_well = {"well1": 500.,
                        "well2": 100.,...}
            # Create stress period data [list of lists per stress period]:
            # wel_spd = {0: [[lay,row,col,discharge1],[lay,row,col,discharge2]]

            # Total discharge per day per well
            Qwell_day = = {"well1": -1000..,
                        "well2": -1.,...}
            
        '''

        if (well_names is None) or (well_names == "None"):
            well_names = []
            # Loop through "well_parameters" of schematisation dictionary
            for iKey in schematisation[dict_key]:   
                # Determine if wells exist (with a "discharge_parameter")  
                if discharge_parameter in schematisation[dict_key][iKey]:
                    well_names.append(iKey)
                    
        # Create empty dicts
        well_loc = {}  # well locations
        KD_well = {}   # KD (cumulative) per well
        Qwell_day = {} # Daily flux [m3/d] per well
        # stress period data for well package
        spd_wel = {}
        spd_wel[0] = []
        for iWell in well_names:
            # Daily flux [m3/d]     
            Qwell_day[iWell] = schematisation[dict_key][iWell][discharge_parameter]
            # Calculate boundary indices
            layidx_min,layidx_max,\
                rowidx_min,rowidx_max,\
                colidx_min,colidx_max = self.cell_bounds(schematisation,
                                        dict_key = dict_key,
                                        dict_subkey = iWell,
                                        model_type = self.model_type)


            print("Laymin_max:", (layidx_min,layidx_max),\
                  "Rowmin_max:", rowidx_min,rowidx_max,\
                  "Colmin_max:",colidx_min,colidx_max)
            print("(nlay,nrow,ncol)",(self.nlay,self.nrow,self.ncol))
            
            # Add well locations and stress_period_data
            well_loc[iWell] = []
            KD_well[iWell] = 0.
            for iLay in range(layidx_min,layidx_max):
                for iRow in range(rowidx_min,rowidx_max):
                    for iCol in range(colidx_min, colidx_max):
                        # print(iLay,iRow,iCol)
                        well_loc[iWell].append((iLay,iRow,iCol))
                        # Correct discharge for K_hor near wells and for the possible difference in delv (K * D)
                        KD_well[iWell] += self.hk[iLay,iRow,iCol] * self.delv[iLay]

            print("Well discharge:", Qwell_day[iWell])
            print("KD_well:", iWell, KD_well[iWell])
            # stress period data for well package
            for iLay in range(layidx_min,layidx_max):
                for iRow in range(rowidx_min,rowidx_max):
                   for iCol in range(colidx_min, colidx_max):
                        spd_wel[0].append([iLay, iRow, iCol, Qwell_day[iWell] * \
                                          (self.hk[iLay,iRow,iCol] * self.delv[iLay]) / KD_well[iWell]])
        
        return well_names,well_loc,KD_well, spd_wel, Qwell_day

    ####################
    ### Fill modules ###
    def load_mfobject(self, mf_exe = 'mf2005.exe'):
        self.mf_exe = mf_exe
        
        # Modflow object
        self.mf = flopy.modflow.Modflow(modelname = self.modelname, exe_name= self.mf_exe, model_ws=self.workspace)

    def load_dis(self, per_nr = 0):
        ''' Add dis Package to the MODFLOW model '''
        perlen = self.perlen[per_nr]
        nstp = self.nstp[per_nr]
        steady = self.steady[per_nr]
        
        self.dis = flopy.modflow.ModflowDis(self.mf, self.nlay, self.nrow, self.ncol, 
                                            nper= 1, lenuni = 2, # meters
                                            itmuni = 4, # 3: hours, 4: days
                                            delr= self.delr, delc= self.delc, laycbd= 0, top= self.top,      
                                            botm= self.bot, perlen = perlen, 
                                            nstp= nstp, steady = steady)
                                            # Laycbd --> 0, then no confining lay below

    def spat_ref(self,xll = 0, yll = 0, delr = None, delc = None,
                 units = "meters",lenuni = 1, length_multiplier=1):
        # Define spatial reference of modflow run
        self.mf.sr = flopy.utils.SpatialReference(xll=0, yll=0, delr=self.mf.dis.delr.array, delc=self.mf.dis.delc.array, 
                                    units='meters', lenuni=1, length_multiplier=1)

    def load_bas(self):
        ''' Add basic Package to the MODFLOW model '''
        self.bas = flopy.modflow.ModflowBas(self.mf, ibound = self.ibound, strt = self.strt)
        
    def load_lpf(self):
        ''' Add lpf Package to the MODFLOW model '''
        self.lpf = flopy.modflow.ModflowLpf(self.mf, layavg = 1, ipakcb = self.iu_cbc, hk=self.hk, vka=self.vka, 
                                            ss = self.ss, storagecoefficient = True) 
        # layavg = 1 (--> logarithmic mean); storagecoefficient = True (means: storativity)
        
    def load_pcg(self, hclose=1e-4, rclose = 0.001):
        ''' Add PCG Package to the MODFLOW model '''
        
        self.pcg = flopy.modflow.ModflowPcg(self.mf, hclose=hclose, rclose = rclose)
            
    def load_oc(self, per_nr = 0):
        ''' Add output control Package to the 
        MODFLOW model. '''            

        extension = ['oc', 'hds', 'ddn', 'cbc', 'ibo']
        unitnumber = [14, 51, 52, self.iu_cbc, 0]
        filenames = []
        for iExt in extension:
            filenames.append(self.modelname + "." + iExt)
        self.oc = flopy.modflow.ModflowOc(self.mf, stress_period_data= self.spd_oc,
                                          extension = extension,  # Default extensions
                                          filenames = filenames,
                                          compact = True #,cboufm='(20i5)'
                                          )                                            
        # compact option is a requirement for using Modpath

    def load_rch(self, rech = None):
        ''' Add recharge Package to the MODFLOW model. '''     
                                           
        # Stress period data dict 
        if rech is not None:                       
            self.rech = rech

        self.rch = flopy.modflow.ModflowRch(self.mf, ipakcb = self.iu_cbc, #102, 
                                     rech = self.rech, nrchop = 1)          

    def load_wel(self, spd_wel = None):
        ''' Add well Package to the MODFLOW model. '''     
                                           
        # Stress period data dict                        
        self.spd_wel = spd_wel

        self.wel = flopy.modflow.ModflowWel(self.mf, ipakcb = self.iu_cbc,
                                     stress_period_data= self.spd_wel)
         # If ipakcb != 0: cell budget data is being saved.

    def write_input_mf(self):
        ''' Write package data MODFLOW model. '''
        self.mf.write_input()

    def run_modflowmod(self):         
        ''' Run modflow model '''
        self.success_mf, _ = self.mf.run_model(silent=False)

    # Here are functions to calculate the travel time through vadose zone, shared functions for
    # Analytical and Modflow models
    def _create_radial_distance_array(self, fraction_flux:np.array or None = None):

        ''' Create array of radial distances from the well to a maximum value, radial distance recharge, which
        is the distance from the well needed to recharge the well to meet the pumping demand. '''
        # ah_todo change this to single array of 0.001 to 100
        # right now we have this set up to directly compare with P. Stuyfzand's results
        # in the 'final' version of the model can be a more fine mesh
        # but need to think about how we change the testing

        self.spreading_distance = math.sqrt(self.schematisation.vertical_resistance_aquitard * self.schematisation.KD)
        # AH do not change to model_radius, since the radial distance for recharge is based on the phreatic value for BOTH cases
        self.radial_distance_recharge =  (math.sqrt(self.schematisation.well_discharge
                                                    / (math.pi * self.schematisation.recharge_rate)))

        if fraction_flux is None:
            fraction_flux = np.array([0.00001, 0.0001, 0.001, 0.005])
            fraction_flux = np.append(fraction_flux, np.arange(0.01, 1, 0.01))
            fraction_flux = np.append(fraction_flux, [0.995, 0.9999])
        
        self.fraction_flux = fraction_flux

        radial_distance = self.radial_distance_recharge * np.sqrt(self.fraction_flux)
            
        if self.schematisation_type == 'semiconfined':

            radial_distance = np.append(radial_distance,
                                    [(radial_distance[-1] + ((self.spreading_distance * 3) - radial_distance[-1]) / 3),
                                    (radial_distance[-1] + 2 *
                                    ((self.spreading_distance * 3) - radial_distance[-1]) / 3),
                                        (self.spreading_distance * 3)])
        
        self.radial_distance = radial_distance

    def particle_data(self, partlocs=None, structured=True, particleids=[0],
                                 localx=None, localy=0.5, localz=None,
                                 timeoffset=0.0, drape=None,
                                 pgname = "Recharge", pg_filename = "Recharge.sloc",
                                 trackingdirection = 'forward',
                                 releasedata=0.0):
        ''' Class to create the most basic particle data type (starting location
        input style 1). Input style 1 is the most general input style and provides
        the highest flexibility in customizing starting locations, see flopy docs. '''
        self.partlocs = partlocs # particle starting locations [(lay,row,col),(k,i,j),...]
        # (structured): Boolean defining if a structured (True) or 
        # unstructured particle recarray will be created
        # (particleids) --> if not None: id-vals of particles added
        self.localx = localx
        self.localy = localy
        self.localz = localz
        # Particle id [part group, particle id] - zero based integers
        if not hasattr(self, "pids"):
            print("Create a new particle id dataset dict.")
            self.pids = {}
        if not hasattr(self, "pg"):
            print("Create a new particle group dataset dict.")
            self.pg = {}
            
        self.pids[pgname] = particleids
        #mp7particledata.
        # Particle distribution package - particle allocation
        #modpath.mp7particledata.Part...
        self.pd = flopy.modpath.ParticleData(partlocs = self.partlocs, structured=True,
                                             drape=0, localx= self.localx, 
                                             localy= self.localy, localz= self.localz,
                                             timeoffset = timeoffset, 
                                             particleids = particleids)

    # def _create_radial_distance_particles_recharge(self, recharge_parameters = None,
    #                                             fraction_flux:np.array or None = None):

    #     ''' Create array of radial distances from the well to a maximum value, radial distance recharge, which
    #     is the distance from the well needed to recharge the well to meet the pumping demand. '''

        
    #     if recharge_parameters is None:
    #         recharge_parameters = self.schematisation_dict.get('recharge_parameters')
        
    #     if fraction_flux is None:
    #         fraction_flux = np.array([0.00001, 0.0001, 0.001, 0.005])
    #         fraction_flux = np.append(fraction_flux, np.arange(0.01, 1, 0.01))
    #         fraction_flux = np.append(fraction_flux, [0.995, 0.9999])
        
    #     self.fraction_flux = fraction_flux

    #     ## Particle group data ##
    #     pgroups = list(recharge_parameters.keys())
    #     # xmin and xmax per pg
    #     self.pg_xmin, self.pg_xmax = {}, {}
    #     # Recharge rate per particle group
    #     self.recharge_pg = {}
    #     # Cumulative flux per particle group
    #     self.flux_pg = {}
    #     # Radial distance array per pg
    #     self.radial_distance = {}
    #     for iPG in pgroups:
    #         # xmin
    #         self.pg_xmin[iPG] = recharge_parameters.get(iPG).get("xmin")
    #         # xmax
    #         self.pg_xmax[iPG] = recharge_parameters.get(iPG).get("xmax")
    #         # recharge rate per particle group
    #         self.recharge_pg[iPG] = recharge_parameters.get(iPG).get("recharge")
    #         # cumulative flux per particle group
    #         self.flux_pg[iPG] = math.pi * (self.pg_xmax[iPG]**2 - self.pg_xmin[iPG]**2) * self.recharge_pg[iPG]
        
    #         # Radial distance array per particle group
    #         self.radial_distance[iPG] = self.pg_xmin[iPG] + ((self.pg_xmax[iPG] - self.pg_xmin[iPG]) * np.sqrt(self.fraction_flux))

    def _create_radial_distance_particles_recharge(self, recharge_parameters = None,
                                                   nparticles_cell: int = 1,
                                                   localy=0.5, localz=0.5,
                                                   timeoffset=0.0, drape=0,
                                                   trackingdirection = 'forward',
                                                   releasedata=0.0):
        ''' Class to create the most basic particle data type (starting location
        input style 1). Input style 1 is the most general input style and provides
        the highest flexibility in customizing starting locations, see flopy docs.
        ###########################################################################
        Create array of radial distances from the well to a maximum value, radial distance recharge, which
        is the distance from the well needed to recharge the well to meet the pumping demand.
        '''
        # SUGGESTION SR: --> add term to dictionary 'recharge_parameters' as 'nparticles_cell'
        
        if recharge_parameters is None:
            recharge_parameters = self.schematisation_dict.get('recharge_parameters')

        ## Particle group data ##
        pgroups = list(recharge_parameters.keys())
        # xmin and xmax per pg
        self.pg_xmin, self.pg_xmax = {}, {}
        for iPG in pgroups:
            # xmin
            self.pg_xmin[iPG] = recharge_parameters.get(iPG).get("xmin")
            # xmax
            self.pg_xmax[iPG] = recharge_parameters.get(iPG).get("xmax")

        try:
            # Determine column indices
            colidx_min = int(np.argwhere((self.xmid >= self.pg_xmin[iPG]) & (self.xmid <= self.pg_xmax[iPG]))[0])
            colidx_max = int(np.argwhere((self.xmid >= self.pg_xmin[iPG]) & (self.xmid <= self.pg_xmax[iPG]))[-1])
            # np.where((self.xmid < right) & (self.xmid > left))    
        except IndexError as e:
            print(e,"Set colidx_min and colidx_max to None.")
            colidx_min, colidx_max = None, None            

                
        # particle group filenames
        self.pg_filenames = {iPG: iPG + ".sloc" for iPG in pgroups}
        # particle starting locations [(lay,row,col),(k,i,j),...]
        self.part_locs = {}  
        # Relative location within grid cells (per particle group)
        self.localx = {}
        self.localy = {}
        self.localz = {}
        # Particle id [part group, particle id] - zero based integers
        if not hasattr(self, "pids"):
            print("Create a new particle id dataset dict.")
            self.pids = {}
        if not hasattr(self, "pg"):
            print("Create a new particle group dataset dict.")
            self.pg = {}

        # Particle dta objects
        self.pd = {}
        # Particle counter
        pcount = -1
        for iPG in pgroups:
            # Particles ids (use counter)
            self.pids[iPG] = []
            # Particle location list per particle group
            self.part_locs[iPG] = []
            # Relative location within grid cells (per particle group)
            self.localx[iPG] = []
            self.localy[iPG] = []
            self.localz[iPG] = []
            for iCol in range(colidx_min,colidx_max):
                # Particle row
                p_row = 0
                # Particle layer
                p_lay = 0
                # Particle column
                p_col = iCol

                for iPart_cell in range(nparticles_cell):
                    # Add particle locations (lay,row,col)
                    self.part_locs[iPG].append((p_lay,p_row,p_col))                    
                    # Relative location of the particles in the cells
                    self.localx[iPG].append((iPart_cell + 0.5)/(float(nparticles_cell)))
                    self.localy[iPG].append(localy)
                    self.localz[iPG].append(localz)
                    # particle count
                    pcount += 1 
                    self.pids[iPG].append(pcount)

            # Particle distribution package - particle allocation
            #modpath.mp7particledata.Part...
            self.pd[iPG] = flopy.modpath.ParticleData(partlocs = self.part_locs[iPG], structured=True,
                                                drape=drape, localx= self.localx[iPG], 
                                                localy= self.localy[iPG], localz= self.localz[iPG],
                                                timeoffset = timeoffset, 
                                                particleids = self.pids[iPG])

            # particle group filename
            self.pg_filename = self.pg_filenames[iPG]
            # particle group object
            # modpath.mp7particlegroup.Part.......
            self.pg[iPG] = flopy.modpath.ParticleGroup(particlegroupname=iPG,
                                                                filename=self.pg_filenames[iPG],
                                                                releasedata=releasedata,
                                                                particledata=self.pd[iPG])
            ''' ParticleGroup class to create MODPATH 7 particle group data for
                location input style 1.  '''
            self.trackingdirection = trackingdirection

            
        

        ## OLD analytical code ##
        # self.spreading_distance = math.sqrt(self.schematisation.vertical_resistance_aquitard * self.schematisation.KD)
        # # AH do not change to model_radius, since the radial distance for recharge is based on the phreatic value for BOTH cases
        # self.radial_distance_recharge =  (math.sqrt(self.schematisation.well_discharge
        #                                             / (math.pi * self.schematisation.recharge_rate)))

        # if fraction_flux is None:
        #     fraction_flux = np.array([0.00001, 0.0001, 0.001, 0.005])
        #     fraction_flux = np.append(fraction_flux, np.arange(0.01, 1, 0.01))
        #     fraction_flux = np.append(fraction_flux, [0.995, 0.9999])
        
        # self.fraction_flux = fraction_flux

        # radial_distance = self.radial_distance_recharge * np.sqrt(self.fraction_flux)
            
        # if self.schematisation_type == 'semiconfined':

        #     radial_distance = np.append(radial_distance,
        #                             [(radial_distance[-1] + ((self.spreading_distance * 3) - radial_distance[-1]) / 3),
        #                             (radial_distance[-1] + 2 *
        #                             ((self.spreading_distance * 3) - radial_distance[-1]) / 3),
        #                                 (self.spreading_distance * 3)])
        
        # self.radial_distance = radial_distance

    ####################################
    ### ModPath 7 input and analyses ###
    ####################################
    ####### Fill Modpath modules #######
    def load_MP7object(self, mp_exe = 'mpath7', mf_model = None,
                      headfilename = None, budgetfilename = None):
        ''' Load modpath7 model, using input from modflow model 'mf_model' (DIS),
            results from head calculations, and cell-by-cell data.'''
        self.mp_exe = mp_exe
        if mf_model is not None:
            self.mf = mf_model
        if headfilename is None:  # Filename of the MODFLOW output head file.
            self.headfile = self.model_hds # os.path.join(self.workspace, 
        else:
            self.headfile = headfilename
        if budgetfilename is None: # Filename of the MODFLOW output cell-by-cell budget file.
            self.cbcfile = self.model_cbc
        else:     
            self.cbcfile = budgetfilename
            
        # Modpath object
        self.mp7 = flopy.modpath.Modpath7(modelname = self.modelname + "_mp", model_ws=self.workspace,
                                             exe_name= self.mp_exe, flowmodel = self.mf)#,
                                             #headfilename = self.headfile, budgetfilename = self.cbcfile)

    def mpbas_input(self, prsity = 0.3, defaultiface = {'RECHARGE': 6}):
        ''' Read model prsity and iface values into object. '''
        self.prsity = prsity
        self.defaultiface = defaultiface
        
    def load_mpbas(self):
        ''' Add BAS Package to the ModPath model '''
        
        self.mpbas = flopy.modpath.Modpath7Bas(model = self.mp7, porosity = self.porosity,
                                                      defaultiface = self.defaultiface)

    def modpath_simulation(self,mp_model = None, trackingdirection = 'backward',
                           simulationtype = 'combined', stoptimeoption = 'extend',
                           particlegroups = None, zones = None):        
        ''' input MODPATH Simulation File Package Class, see flopy docs. '''
        if mp_model is not None:
            self.mp7 = mp_model
            
        if self.trackingdirection is None:
            self.trackingdirection = trackingdirection
            
            
        # Zones are not read in detail
        if zones is None:
            self.zones = [1] * self.nlay
        else:
            self.zones = zones
        
        self.mp7sim = flopy.modpath.Modpath7Sim(model = self.mp7, mpnamefilename=None, 
                                         listingfilename=None, endpointfilename=None,
                                         pathlinefilename=None, timeseriesfilename=None,
                                         tracefilename=None, simulationtype = simulationtype,
                                         trackingdirection = self.trackingdirection, 
                                         weaksinkoption='stop_at', weaksourceoption='pass_through',
                                         budgetoutputoption='summary', traceparticledata=None, #[0,0], #self.pid, 
                                         budgetcellnumbers=None, referencetime=0.,
                                         stoptimeoption = stoptimeoption, stoptime=None,
                                         timepointdata=None, zonedataoption='off',  # timepointdata=[100*24,1/24.]
                                         stopzone='off', zones=self.zones, retardationfactoroption='off',
                                         retardation=1.0, particlegroups=particlegroups,
                                         extension='mpsim')
        # stoptimeoption: "extend" --> particle simulation continues beyond time specified in modflow BAS package
        # timepointdata = [100*24,1/24.] # max 100 days of data

    def write_input_mp(self):
        ''' Write package data ModPath model. '''
        self.mp7.write_input()

    def run_ModPathmod(self):         
        ''' Run ModPath model '''
        # self.mp7.run_model(silent=False)
        self.success_mp,_ = self.mp7.run_model(silent=False)

### Functie: Check xmin, xmax,
    '''

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
  	    self.WELL_package  # np.array
  		self.Q
  		self.well_id

    	def assign_multi_nodal_well():
  	    self.wellid  # np.array
  		self.....
	def assign_ibound():
        self.ibound  # 3 D grid
  			self.head
    '''


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
    
    def export_to_df(self, grid_material, what_to_export='all'):
  	    """ Export to dataframe....

        Parameters
        ----------
        what_to_export: String
        		options: 'all', 'omp', 'microbial_parameters'
        """
        # df_flowline = pd.DataFrame()
        # df_particle = pd.DataFrame()

  		#delete the unwanted columns depending on what the user asks for here
  		# return df_flowline, df_particle
    '''
    def _export_to_df(self,
        df_output,
        distance,
        total_travel_time,
        travel_time_unsaturated,
        travel_time_shallow_aquifer,
        travel_time_target_aquifer,
        discharge_point_contamination=None, #AH_todo do we need this?
        ):
        """ Makes 'df_flowline' and 'df_particle'
        # MK: attributes, dont make thema rguments
        #AH_todo review the logic of how the diffusion/point sources are calculated with @martinK

            Parameters
            ----------
            df_output: pandas.DataFrame
                Column 'total_travel_time': float
                Column 'travel_time_unsaturated': float
                Column 'travel_time_shallow_aquifer': float
                Column 'travel_time_target_aquifer': float
                Column 'radial_distance': float
                Column 'head': float
                Column 'cumulative_fraction_abstracted_water': float
                Column 'flowline_discharge': float
            distance: array
                Array of distance(s) [m] from the well.
                For diffuse sources 'distance' is the 'radial_distance' array.
                For point sources 'distance' is 'distance_point_contamination_from_well'.
            total_travel_time: array
                Sum of the unsaturated, shallow and target aquifer travel times for each
                point in the given distance array, [days].
            travel_time_unsaturated: array
                Travel time in the unsaturated zone for each point in the given distance array returned as
                attrubute of the function, [days].
            travel_time_shallow_aquifer: array
                Travel time in the shallow aquifer for each point in the given distance array, [days].
            travel_time_target_aquifer: array
                Travel time in the target aquifer for each point in the given distance array, [days].

            Returns
            -------
            df_flowline: pandas.DataFrame
                Column 'flowline_id': Integer
                Column 'flowline_type': string
                    Described the type of contamination associated with the flowline,
                    either 'diffuse_source' or 'point_source'.
                Column 'discharge': Float
                    Discharge associated with the flowline, [m3/d].
                Column 'particle_release_date': Float
                Column 'input_concentration'
                Column 'endpoint_id': Integer
                    ID of Well (or drain) where the flowline ends.
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

            df_particle: pandas.DataFrame
                Column 'flowline_id': int
                Column 'zone': string
                    Zone in the aquifer, ground surface 'surface', 'vadose_zone', 'shallow_aquifer' or 'target_aquifer'
                Column 'travel_time_zone': float
                    Travel time in the respective aquifer zone given in column 'zone.
                Column 'xcoord': float
                Column 'ycoord': float
                Column 'zcoord': float
                Column 'redox_zone': float
                    'suboxic', 'anoxic', deeply_anoxic'
                Column 'temperature': float
                    Of the respective aquifer zone.
                Column 'thickness_layer': float
                Column 'porosity_layer': float
                Column 'dissolved_organic_carbon': float
                Column 'pH': float
                Column 'fraction_organic_carbon': float
                Column 'solid_density_layer': float
            """
        # MK:  this is something that I could imagine that you use it as an argument to the
        # function.
        #AH_todo, @MartinK, what is the advantage over using the attribute?
        what_to_export = self.schematisation.what_to_export


    def load_mfmodel(self, fname_nam):
        ''' Load modflow model from namefile. '''
        self.mf = flopy.modflow.Modflow.load(fname_nam)

    def mfmodelrun(self):
        ''' Main model run code. '''
        print ("Run model:",self.workspace, self.modelname +"\n")

        # Start model run
        # Load modflow packages
        for iper in range(self.nper):
            # If the starting heads have been changed during the first model run:
            if iper != 0:

                try:
                    # Load head data
                    _, self.strt = self.load_hdsobj(fname = self.model_hds,time = -1)
                except Exception as e:
                    print (e,"Error loading strt_fw data\nstrt_fw set to '0'.")
                    self.strt = 0.
            
            # Open modflow object
            self.load_mfobject(mf_exe = 'mf2005.exe')
            # Load packages
            self.load_dis(per_nr = iper)   
            self.load_bas()   
            self.load_lpf()   
            self.load_pcg()   
            self.load_oc(per_nr = iper)
            if len(self.well_names) > 0:
                try:                        
                    self.load_wel(spd_wel = self.spd_wel)
                except Exception as e:
                    print(e, "error loading well package.")
            try:
                self.load_rch(rech = self.recharge)
            except Exception as e:
                print(e, "no recharge assigned.")

            # Assign spatial reference:
            # this is a fix until the grid object is available                    
            self.spat_ref(xll = self.xll, yll = self.yll, delr = self.delr,
                     delc = self.delc, units = "meters",
                     lenuni = 1, length_multiplier=1)

            # Write modflow, mt3dms and seawat input
            try: 

                self.write_input_mf()
            finally:
                pass
            
            # Finally removes the given exceptions from memory
        
            # Try to delete the previous output files, to prevent accidental use of older files
            try:  
                os.remove(self.model_hds)
            except FileNotFoundError:
                pass
            try:  
                os.remove(self.model_cbc)
            except FileNotFoundError:
                pass
            
            # Run modflow model
            try:
                self.run_modflowmod()
            except Exception as e:
                self.success_mf = False
                print(e, self.success_mf)
            # print(self.success_mf, self.buff)

            # Model run completed succesfully
            print("Model run", self.workspace, self.modelname, "completed without errors:", self.success_mf)

    def MP7modelrun(self, mf_namfile = None, mp_exe = None):
        ''' Modpath model run.'''
        print ("Run modpath:",self.workspace, self.modelname +"\n")

        # Copy mfmodel to modpath section
        if mf_namfile is None:  
            if hasattr(self, "mf"):
                print("Modflow model 'mf' already exists in object.")
                pass
            else:
                # Load mfmodel assuming namfile exists in same folder
                self.mf_namfile = os.path.join(self.workspace, self.modelname + '.nam')
                self.load_mfmodel(self.mf_namfile)
        else:
            if hasattr(self, "mf"):
                # modflow model "mf" already exists in object
                pass
            else:
                print("Load mf model using given namfile location:", mf_namfile)
                # mf_namefile
                self.mf_namfile = mf_namfile
                self.load_mfmodel(self.mf_namfile)

        
        # Load the MP7 module/object
        if mp_exe is None:
            self.mp_exe = "mpath7"
        else:
            self.mp_exe = mp_exe

        # Load mp7 object into Aximodel class
        self.load_MP7object(mp_exe = self.mp_exe, mf_model = self.mf)#,

        # Load modpath basic package
        self.load_mpbas()
        
        # Select particle groups as input to the model
        self.particlegroups = []
        for iPG in self.pg:
            self.particlegroups.append(self.pg[iPG])
        # Only for unstruct grids:
        #- Load modpath unstructured grid
        #- Load modpath time discretization
        
        # Write modpath simulation File:
        self.modpath_simulation(mp_model = None,# trackingdirection = None,
                           simulationtype = 'combined', stoptimeoption = 'extend',
                           particlegroups = self.particlegroups) 

        # Write files and execute the modpath model
        self.write_input_mp()

        try:
            self.run_ModPathmod()
            # Model run completed succesfully
            print("ModPath run", self.workspace, self.modelname, "completed succesfully.")
        except Exception as e:
            self.success_mp = False
            print(e, "ModPath run", self.workspace, self.modelname, "failed.")           

    def phreatic(self):
        ''' Modflow & modpath calculation using subsurface schematisation type 'phreatic'. '''
        self._check_init_phreatic()
        # Model_type axisymmetric
        self.model_type = "axisymmetric"

        # Make radial discretisation
        # Use dictionary keys from schematisation
        dict_keys = ["geo_parameters","recharge_parameters","ibound_parameters",
                      "well_parameters","mesh_refinement"]
        self.make_discretisation(schematisation = self.schematisation_dict, dict_keys = dict_keys,
                            model_type = self.model_type)

        # Set ibound grid and starting head
        dict_keys = ["ibound_parameters"]
        self.strt = np.zeros((self.nlay,self.nrow,self.ncol), dtype = 'float')
        self.set_ibound(schematisation = self.schematisation_dict,
                            dict_keys = dict_keys,
                            ibound = None,
                            strt = self.strt)


        # list active packages
        active_packages = []
        # Parameter requirement
        package_parms = {"BAS": "ibound"}

        ''' ### lpf package input parms ###
            
            hk: Horizontal conductivity
            vka: Vertical conductivity (-> if layvka = 0 (default))
            ss: Specific storage (1/m)
            storativity: # If True (default): Ss stands for storativity [-] instead of specific storage [1/m]
            layavg: Layer average of hydraulic conductivity
                layavg: 0 --> harmonic mean 
                layavg: 1 --> logarithmic mean (is used in axisymmetric models).
        '''

        if self.model_type == "axisymmetric":
            self.layavg = 1
        else:  
            self.layavg = 0
        # geohydrological parameter names # list of schematisation_dict terms & dtype
        self.geoparm_names = {"moisture_content": [["geo_parameters"],"float"],
                        "hk": [["geo_parameters"],"float"],
                        "vani": [["geo_parameters"],"float"],
                        "porosity": [["geo_parameters"],"float"],
                        "solid_density": [["geo_parameters"],"float"],
                        "fraction_organic_carbon": [["geo_parameters"],"float"],
                        "redox": [["geo_parameters"],"object"],
                        "dissolved_organic_carbon": [["geo_parameters"],"float"],
                        "pH": [["geo_parameters"],"float"],
                        "temperature": [["geo_parameters"],"float"]
                        }
                      

        for iParm, dict_keys in self.geoparm_names.items():
            # Temporary value grid
            unitgrid = np.ones((self.nlay,self.nrow,self.ncol), dtype = dict_keys[1])
            # Fill grid with new values
            grid = self.fill_grid(schematisation = self.schematisation_dict,
                                dict_keys = dict_keys[0],
                                parameter = iParm,
                                grid = unitgrid,
                                dtype = dict_keys[1])
            self.update_property(property = iParm, value = grid)

        # Assign material grid
        self.assign_material(schematisation = self.schematisation_dict,
                            dict_keys = ["geo_parameters"])

        # Create (uncorrected) array for kv ("vka"), using "kh" and "vani" (vertical anisotropy)
        self.hk[self.material == None] = 999.
        self.vani[self.material == None] = 999.
        # Vertical conductivity [m/d]
        self.vka = self.hk / self.vani
        # and for 'storativity' ("ss": specific storage)
        self.ss = np.ones((self.nlay,self.nrow,self.ncol), dtype = 'float') * 1.E-6
        # Axisymmetric flow properties
        axisym_parms = ["hk","vka","ss"]
        if self.model_type == "axisymmetric":
            for iParm in axisym_parms:
                grid_uncorr = getattr(self,iParm)
                grid_axi = self.axisym_correction(grid = grid_uncorr)
                # Update attribute
                self.update_property(property = iParm, value = grid_axi)

        # Create input oc package
        self.oc_input(spd_oc = {(0, 0): ['save head', 'save budget']})

        # Create recharge package
        rech_parmnames = {"recharge": [["recharge_parameters"],"float"]}
        for iParm, dict_keys in rech_parmnames.items():
            # Temporary value
            grid = self.fill_grid(schematisation = self.schematisation_dict,
                            dict_keys = dict_keys[0],
                            parameter = iParm,
                            grid = None,
                            dtype = dict_keys[1])
            
            if self.model_type == "axisymmetric":
                grid_axi = self.axisym_correction(grid = grid)[0,:,:]
                # Update attribute
                self.update_property(property = iParm, value = grid_axi)
            else:
                # Update attribute
                self.update_property(property = iParm, value = grid[0,:,:])

        # Well input
        # !!! Obtain node numbers of well locations (use indices) !!!
        # leakage discharge from well
        well_names = [iWell for iWell in self.schematisation_dict["well_parameters"] if \
                        "Q" in self.schematisation_dict["well_parameters"][iWell].keys()]
       
        self.well_names,\
            self.well_loc,\
                self.KD_well,\
                    self.spd_wel,\
                        self.Qwell_day = self.assign_wellloc(schematisation = self.schematisation_dict,
                                                        dict_key = "well_parameters",
                                                        well_names = None,
                                                        discharge_parameter = "Q")

        # Load wel parms to model
        # self.wel_input(spd_wel = self.spd_wel)

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

    def load_hdsobj(self, fname = None, time = None):
        
        ''' Return head data from file.
            If time = -1 --> return final time and head grid,
            elif time = 'all' --> return all time values and head grids,
            elif time = [1.,2.,time_n]--> Return head grids for prespecified times. '''
        
        try:
            # Read binary concentration file
            hdsobj = bf.HeadFile(fname, precision = 'single', verbose = False)
            times = hdsobj.get_times()
            head_dat = {}
                
            if time == -1:
                head_dat = hdsobj.get_data(totim = times[-1])
            elif time == 'all':
                for iTime in times:
                    head_dat[iTime] = hdsobj.get_data(totim = iTime)
            else:
                try:
                    for iTime in time:
                        head_dat[iTime] = hdsobj.get_data(totim = iTime)
                except Exception as e:
                    print ("time values are not in saved list of times")
        except Exception as e:
            print ("Error loading head data\nhead_dat set to '0'.")
            head_dat = 0.
            pass
        finally:
            hdsobj.close()

        return times, head_dat

    def read_binarycbc(fname):
        ''' Read binary cell budget file (fname). 
            This is modflow output.'''
        cbcobj = bf.CellBudgetFile(fname)
        print(cbcobj.list_records())
        
        frf = cbcobj.get_data(text='FLOW RIGHT FACE')[0]
        flf = cbcobj.get_data(text='FLOW LOWER FACE')[0]
        cbcobj.close()
        
        return frf, flf 

    def get_nodes(self,locs):
        ''' Obtain/return model node index (int) belonging to
            layer 'iLay', row 'iRow' and column 'iCol'.
            for point locations 'locs'. '''
        nodes = []
        for iLay,iRow,iCol in locs:
            nodes.append(iLay * self.nrow * self.ncol + iRow * self.ncol + iCol)
        return nodes

    def calc_node_indices(self,xyz_nodes, particle_list: list or None = None):
        ''' Obtain/return layer,row,column idx ("node_indices") as dict of np.arrays
            corresponding to xyz-coördinates of type dict per tracked particle (as key)
            with nodes "xyz_nodes" (np.array). 
            
            Method requires center points xmid, ymid, zmid to be predefined.
            
        '''

        if particle_list is None:
            particle_nodes = list(xyz_nodes.keys())
        else:  # requires check if all indices occur in 'xyz_nodes'
            particle_nodes = [idx for idx in particle_list if idx in xyz_nodes.keys()]
            if len(particle_nodes) < particle_list:
                print("Warning: particles do not match 'xyz_nodes' indices.\n",
                "Function uses all 'xyz_nodes' keys instead (=default).")
                particle_nodes = list(xyz_nodes.keys())

        # Create dict for Layer, row, column indices per particle
        node_indices = {}
        for iPart in particle_nodes:
            # Number of nodes
            nr_nodes = xyz_nodes[iPart].shape[0]
    #                                    print(iPart)
            node_indices[iPart] = [(np.argmin(abs(self.zmid-xyz_nodes[iPart][iNode][2])), \
                                        np.argmin(abs(self.ymid-xyz_nodes[iPart][iNode][1])), \
                                        np.argmin(abs(self.xmid-xyz_nodes[iPart][iNode][0]))) \
                                                        for iNode in range(nr_nodes)]


        return node_indices

    def read_pathlinedata(self,fpth, nodes):
        ''' read pathlinedata from file fpth (extension: '*.mppth'), 
            given the particle release node index 'nodes' obtained from 
            a tuple or list of tuples (iLay,iRow,iCol).
            Cell node indices can be obtained using the method get_nodes((iLay,iRow,iCol))
            
            Return xyz_nodes, dist, tdiff, dist_tot, time_tot, pth_data 
            
            with xyz_nodes, the xyz coördinates of each node per tracked particle
            dist (L): the distance between each node per tracked particle
            tdiff (T): the travel time (difference) between each node per particle
            dist_tot: total distance cvered by each tracked particle
            time_tot: total duration between release and ending of each particle
            pth_data: complete recarray is returned to see what's in the file fpth.
            
                '''
        pth_object = flopy.utils.PathlineFile(fpth)
        # Raw pathline data
        pth_data = pth_object.get_destination_pathline_data(nodes)
        time, xyz_nodes, txyz, dist, tdiff = {}, {}, {}, {}, {}
        # number of particles within file
        npart = len(pth_data)
        # Numpy arrays storing the total distance traveled for each particle
        dist_tot = np.zeros((npart), dtype = 'float')
        # Numpy arrays storing the total travel time for each particle
        time_tot = np.zeros((npart), dtype = 'float')
        for idx,iPart in enumerate(pth_data):
            # n_nodes (ruw)
            n_nodes_raw = len(pth_data[idx]["x"])
            # XYZ data
            txyz[idx] = np.empty((n_nodes_raw,4), dtype = 'float')
            txyz[idx] = np.array([pth_data[idx]["time"],
                                pth_data[idx]["x"],
                                pth_data[idx]["y"],
                                pth_data[idx]["z"],
                                ]).T
            # Remove identical data (based on identical times)
            txyz[idx] = np.unique(txyz[idx], axis = 0)
            # time data
            time[idx] = txyz[idx][:,0]
            # xyz data
            xyz_nodes[idx] = txyz[idx][:,1:]
            
            # Determine number of remaining nodes
            n_nodes = min(xyz_nodes[idx].shape[0],time[idx].shape[0])
            if xyz_nodes[idx].shape[0] != time[idx].shape[0]:
                print(xyz_nodes[idx].shape[0],time[idx].shape[0])
            #Distance array between nodes
            dist[idx] = np.zeros((n_nodes-1), dtype = 'float')
            # Time difference array
            tdiff[idx] = np.zeros((n_nodes-1), dtype = 'float')
            for iNode in range(1,n_nodes):
                # Calculate internodal distance (m)
                dist[idx][iNode-1] = np.sqrt((xyz_nodes[idx][iNode][0] - xyz_nodes[idx][iNode-1][0])**2 + \
                    (xyz_nodes[idx][iNode][1] - xyz_nodes[idx][iNode-1][1])**2 + \
                    (xyz_nodes[idx][iNode][2] - xyz_nodes[idx][iNode-1][2])**2)
                # Calculate time difference between nodes
                tdiff[idx][iNode-1] = (time[idx][iNode] - time[idx][iNode - 1])
                
            # Total distance covered per particle
            dist_tot[idx] = dist[idx].sum()
            # Total time covered per particle
            time_tot[idx] = time[idx][-1]
                        
        return xyz_nodes, dist, tdiff, dist_tot, time_tot, pth_data

    # Make df_particle
    def fill_df_particle (df,
                        xyz_nodes,
                        dist,
                        tdiff,
                        dist_tot,
                        time_tot
                        # distance,
                        # travel_time_unsaturated,
                        # travel_time_shallow_aquifer,
                        # travel_time_target_aquifer,
                        # total_travel_time
                        ):
        '''Fill the df_particle by flowline_id

        Parameters
        ----------
        df: pandas.dataframe
            Holder for df_particle, filled by the flowline_id
        distance: array
            Array of distance(s) [m] from the well.
            For diffuse sources 'distance' is the 'radial_distance' array.
            For point sources 'distance' is 'distance_point_contamination_from_well'.
        travel_time_unsaturated: array
            Travel time in the unsaturated zone for each point in the given distance array returned as
            attrubute of the function, [days].
        travel_time_shallow_aquifer: array
            Travel time in the shallow aquifer for each point in the given distance array, [days].
        travel_time_target_aquifer: array
            Travel time in the target aquifer for each point in the given distance array, [days].
        total_travel_time: array
            Sum of the unsaturated, shallow and target aquifer travel times for each
            point in the given distance array, [days].

        Returns
        -------
        df: pandas.dataframe
            Holder for df_particle, filled by the flowline_id

        flowline_id: ID of the flowline
        zone
        travel_time
        total_travel_time
        xcoord
        ycoord
        zcoord
        redox
        temperature
        travel_distance
        porosity
        dissolved_organic_carbon
        pH
        fraction_organic_carbon
        solid_density

        

        df.loc[0] = [flowline_id,
                    "surface",
                    0,
                    0,
                    distance,
                    self.schematisation.model_width,
                        self.schematisation.ground_surface,
                        None,
                        self.schematisation.temperature,
                        None,
                        None,
                        None,
                        None,
                        None,
                        None]

        df.loc[1] = [flowline_id,
                        "vadose_zone",
                        travel_time_unsaturated,
                        travel_time_unsaturated,
                        distance,
                        self.schematisation.model_width,
                        self.schematisation.bottom_vadose_zone_at_boundary, # @MartinvdS should this be the thickness_vadose_zone_drawdown??
                        self.schematisation.redox_vadose_zone,
                        self.schematisation.temperature,
                        self.schematisation.thickness_vadose_zone_at_boundary,
                        self.schematisation.porosity_vadose_zone,
                        self.schematisation.dissolved_organic_carbon_vadose_zone,
                        self.schematisation.pH_vadose_zone,
                        self.schematisation.fraction_organic_carbon_vadose_zone,
                        self.schematisation.solid_density_vadose_zone,
                        ]

        df.loc[2] = [flowline_id,
                    "shallow_aquifer",
                        travel_time_shallow_aquifer,
                        travel_time_unsaturated + travel_time_shallow_aquifer,
                        distance,
                        self.schematisation.model_width,
                        self.schematisation.bottom_shallow_aquifer,
                        self.schematisation.redox_shallow_aquifer,
                        self.schematisation.temperature,
                        self.schematisation.thickness_shallow_aquifer, # @MartinvdS does this need to account for the vadose zone drawdown?
                        self.schematisation.porosity_shallow_aquifer,
                        self.schematisation.dissolved_organic_carbon_shallow_aquifer,
                        self.schematisation.pH_shallow_aquifer,
                        self.schematisation.fraction_organic_carbon_shallow_aquifer,
                        self.schematisation.solid_density_shallow_aquifer,
                        ]

        df.loc[3] = [flowline_id,
                    "target_aquifer",
                        travel_time_target_aquifer,
                        total_travel_time,
                        self.schematisation.diameter_borehole/2, #at the well
                        self.schematisation.model_width,
                        self.schematisation.bottom_target_aquifer,
                        self.schematisation.redox_target_aquifer,
                        self.schematisation.temperature,
                        self.schematisation.thickness_target_aquifer,
                        self.schematisation.porosity_target_aquifer,
                        self.schematisation.dissolved_organic_carbon_target_aquifer,
                        self.schematisation.pH_target_aquifer,
                        self.schematisation.fraction_organic_carbon_target_aquifer,
                        self.schematisation.solid_density_target_aquifer,
                        ]
        return df
        '''

    def run_model(self,
                    # simulation_parameters: dict or None = None,
                    xll = 0., yll = 0., perlen:dict or float or int = 365.*50, 
                    nstp:dict or int = 1, nper:int = 1,
                    steady:dict or bool = True,
                    run_mfmodel = True, run_mpmodel = True):
        ''' Run the combined modflow and modpath model using one 
            of four possible schematisation types:
            - "Phreatic"
            - "Semi-confined"
            - "Recharge basin (BAR)"
            - "River bank filtration (RBF)"
            Currently (13-7-2021) only the Phreatic schematisation is supported.'''

        # print(self.schematisation)
        # Run modflow model (T/F)
        self.run_mfmodel = run_mfmodel
        # Run modpath model (T/F)
        self.run_mpmodel = run_mpmodel
        # simulation parameters
        self.simulation_parameters = self.schematisation_dict.get('simulation_parameters')
        # if simulation_parameters is None:
        #     self.simulation_parameters = self.schematisation_dict.simulation_parameters
        # else:
        #     self.simulation_parameters = simulation_parameters

        # Simulation parameters
        # Dict with stress period lengths
        if type(perlen) != "dict":
            self.perlen = {0: perlen}    
        else:
            self.perlen = perlen
        # Nr of time periods per stress period (int)
        if type(nstp) != "dict":
            self.nstp = {0: nstp} 
        else:
            self.nstp = nstp 

        # Nr of stress periods  (int)
        self.nper = nper

        # Steady state model run (True/False)
        if type(steady) != "dict":
            self.steady = {0: steady} 
        else:    
            self.steady = steady  

        # Define reference lowerleft
        self.xll = xll
        self.yll = yll

        # Type of scenario
        try:
            self.schematisation_type = self.simulation_parameters.get("schematisation_type")
        except KeyError as e:
            print(e)

        if self.schematisation_type == "phreatic":
            print("Run phreatic model")
            self.phreatic()
            if self.run_mfmodel:
                # Run modflow model
                self.mfmodelrun()
        # elif self.schematisation_type == "Semi-confined":

        # Modpath simulation
        if self.run_mpmodel:

            # Create radial distance array with particle locations
            # self._create_radial_distance_array() # Analytische fluxverdeling
            self._create_radial_distance_particles_recharge(recharge_parameters = None,
                                                            nparticles_cell = 1,
                                                            localy=0.5, localz=0.5,
                                                            timeoffset=0.0, drape=0,
                                                            trackingdirection = 'forward',
                                                            releasedata=0.0)  

            # Default flux interfaces
            defaultiface = {'RECHARGE': 6, 'ET': 6}
            
            
            ## Model mpbas input ##
            self.mpbas_input(prsity = self.porosity, defaultiface = defaultiface)

            # Run modpath model
            self.MP7modelrun()

            # self.success_mp = True
    
            print("modelrun of type", self.schematisation_type, "completed.")

            ## Post-processing: ##
            # ##cbc
            # cbcfile = os.path.join(self.workspace, self.modelname + '.cbc') # r"r:\P402045_014\microbiologisch_risico_lekke_peilbuis\python\MP7_V5_20200429\GHscen_A_onvzone_afwezig_lek_klein_lekdiepte0_5\GHscen_A_onvzone_afwezig_lek_klein_lekdiepte0_5.cbc"
            # cbb = flopy.utils.binaryfile.CellBudgetFile(filename = cbcfile,
            #                                             precision='single',
            #                                             verbose=False)
            
            # # Flux per particle group, per path
            # fluxnode_frf = {}
            # frf = cbb.get_data(kstpkper=(0,0), text='FLOW RIGHT FACE')
            # for iPG in self.pg:
            #     fluxnode_frf[iPG] = {}
            #     for id_,iNode in enumerate(self.part_locs[iPG]):
            #         print(id_,iNode, "right face flux:", frf[0][iNode])
            #         fluxnode_frf[iPG][iNode] = frf[0][iNode]
            # # close budget object afterwards
            # cbb.close()

            # Empty output dicts
            self.xyz_nodes = {}
            self.time_diff = {}
            self.particle_data = {}
            # dataframes of particle data (dict)
            df_particle = {}
            # list the particle_data dataframes
            df_particle_list = []
            # Pathline output file
            self.mppth = os.path.join(self.workspace, self.modelname + '_mp.mppth')

            for iPG in self.part_locs:  # use endpoint_id dict or list
                # xyz_locs
                self.xyz_nodes[iPG] = {}
                # Save flow duration (time_diff) of pathlines
                self.time_diff[iPG] = {}
                # group recarray per particle location
                self.particle_data[iPG] = {}


                # Nodes to retrieve
                print("Particle group",iPG, "nr of nodes:", str(len(self.part_locs.get(iPG))))
                for id_,iNode in enumerate(self.part_locs.get(iPG)):
                    # print(id_,iNode)
                    try:
                        self.nodes = self.get_nodes(iNode)
                    except Exception:
                        self.nodes = self.get_nodes([iNode])
                    print(id_,iNode, self.nodes)

                    # nodes_rel[ = flopy.utils.ra_slice(m.wel.stress_period_data[0], ['k', 'i', 'j'])
                    #     nodes_well = prf.get_nodes(locs = wel_locs, nrow = nrow, ncol = ncol) # m.dis.get_node(wel_locs.tolist())

                    # Read pathline data
                    self.xyz_nodes[iPG][iNode], \
                    self.dist_data,  \
                    self.time_diff[iPG][iNode],  \
                    self.dist_tot,   \
                    self.time_tot,   \
                    self.pth_data =  \
                                self.read_pathlinedata(fpth = self.mppth,
                                                    nodes = self.nodes) #pg_nodes[iGroup])  
                    
                    '''                             
                    xyz_points, \    # XYZ data 
                    dist_data,  \    # Distance araay between nodes
                    time_diff,  \    # Time difference array
                    dist_tot,   \    # Total distance covered per particle
                    time_tot,   \    # Total time covered per particle
                    pth_data =  \    # Raw pathline data
                    '''        
                    # col, lay, row index
                    node_indices = self.calc_node_indices(xyz_nodes = self.xyz_nodes[iPG][iNode])
                    # Particle indices
                    part_idx = [iPart for iPart in self.xyz_nodes[iPG][iNode]] #node_indices]
                    # Test array travel times (days)
                    tot_time_arr = {iPart: self.time_diff[iPG][iNode][iPart].sum() for iPart in part_idx}
                    # Fill recarray
                    self.particle_data[iPG][iNode] = copy.deepcopy(self.pth_data)
                    # Create rec.arrays for porosity, pH, T, etc. to append to particle_data
                    parm_list = ["porosity","solid_density","fraction_organic_carbon",
                                "redox","dissolved_organic_carbon", "pH","temperature","material"]
                    for iPart in part_idx:

                        # Loop through rec.arrays of pth_data using part_idx 
                        for iParm in parm_list:


                        
                            # Material property array
                            material_property_arr = getattr(self,iParm)
                            # dtype of array
                            mat_dtype = material_property_arr.dtype
                            if mat_dtype == 'object':
                                dtype_ = '|S20'
                            else:
                                dtype_ = mat_dtype
                                
                            # Numpy array values
                            parm_values = np.array([material_property_arr[idx[0],idx[1],idx[2]] for idx in node_indices[iPart]],
                                                    dtype = [(iParm,mat_dtype)])
                            # Numpy recarray of material property
                            parm_values_rec = parm_values.view(np.recarray)
                            # Append recarray to particle_data (dict of dicts of np.recarray)
                            self.particle_data[iPG][iNode][iPart] = rfn.rec_append_fields(self.particle_data[iPG][iNode][iPart], iParm, parm_values_rec, 
                                                                        dtypes=dtype_)
                            # self.particle_data[iPG][iNode][iPart] = np.append(self.particle_data[iPG][iNode][iPart],parm_values_rec)
            
                        # Export rec.arrays as pd dataframe
                        df_particle[f"pg: {iPG} node: {iNode} particle: {iPart}"] = pd.DataFrame.from_records(data = self.particle_data[iPG][iNode][iPart],
                                                                                    index = "particleid",
                                                                                    exclude = None).iloc[:-1,:]

                        particle_subset_fname = os.path.join(self.dstroot,f"pg{iPG}_node{iNode}_particle{iPart}.csv")
                        df_particle[f"pg: {iPG} node: {iNode} particle: {iPart}"].to_csv(particle_subset_fname)
                        # Append dataframes
                        df_particle_list.append(df_particle[f"pg: {iPG} node: {iNode} particle: {iPart}"])

            # Concatenate df_particle dataframes ('combined')
            df_particle_comb = pd.concat(df_particle_list, axis = 0, ignore_index = False)

            # df_particle file name
            particle_fname = os.path.join(self.dstroot,self.schematisation_type + "_df_particle.csv")
            df_particle_comb.to_csv(particle_fname)
            
            print("Post-processing modpathrun of type", self.schematisation_type, "completed.")


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
        		options: 'all', 'omp', 'microbial_parameters'
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
class SubstanceTransport():
    """ Returns concentration in a groundwater well for a given Organic Micro Pollutant or microbial species.

    Parameters
    ----------
    df_flowline: pandas.DataFrame
        Column 'flowline_id': Integer
        Column 'discharge': Float
            Discharge associated with the flowline (m3/d)
        Column 'particle_release_day': Float
        Column 'input_concentration'
        Column 'start_or_end_point_id': Integer ######################################
        		ID of Well (or drain) where the flowline ends

    df_particle: pandas.DataFrame
        Column 'flowline_id'
        Column 'travel_time'
        Column 'xcoord'
        Column 'ycoord'
        Column 'zcoord'
        Column 'redox'
        Column 'temperature'
        Column 'Kow'  # only necessary for OMP
        Column 'Labda'  # only necessary for pathogen

    Returns
    -------    
		

    """

    def __init__(self, substance: Substance, df_particle, df_flowline, removel_function?):
        self.omp_inialized = False

  
  	def _init_omp()
  		if self.omp_inialized:
	  		self.df_part['Kow'] = self.df_part['redox'].apply(lambda x: substance.get_Kow(x)
   		self.omp_inialized = True


  	def _init_pathogen()


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

    def compute_pathogen_removal(self):
                                                             
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
concentration = SubstanceTransport()
if removal_function == 'omp':
		concentration.compute_omp_removal
elif removal_function = 'omp':
		concentration.compute_pathogen_removal
else:
  	raise KeyError('schematisation argument not recognized')
                       
'''

#%%
if __name__ == "__main__":
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
        'pathogen' -> get parameters for pathogen
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
    #         		options: 'all', 'omp', 'microbial_parameters'
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
    # # output Alex "phreatic_dict_nogravel.txt" --> saved in "testing dir"
    # research_dir = os.path.join(path,"..","research")
    # with open(os.path.join(research_dir,"phreatic_dict_nogravel.txt"),"r") as file_:
    #     dict_content = file_.read()
    #     phreatic_scheme = eval(dict_content)

    # check_schematisation = False # Check input dict (T/F)
    # if check_schematisation:
    #     for iKey,iVal in phreatic_scheme.items():
    #         print(iKey,iVal,"\n")

