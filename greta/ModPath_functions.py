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

# Plotting modules
import matplotlib.pyplot as plt
import matplotlib 
import matplotlib.colors as colors

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

from greta.Analytical_Well import AnalyticalWell

# try:
#     from greta.Analytical_Well import * 
#     from greta.Substance_Transport import *
# except ModuleNotFoundError as e:
#     print(e, ": second try.")
#     module_path = os.path.join("..","greta")
#     if module_path not in sys.path:
#         sys.path.insert(0,module_path)
#     from Analytical_Well import * 
#     from Substance_Transport import *

#     print("Second try to import modules succeeded.")

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
    def __init__(self, schematisation: AnalyticalWell,
                       workspace: str or None = None, modelname: str or None = None,
                       bound_left: str = "xmin", bound_right: str = "xmax",
                       bound_top: str = "top", bound_bot: str = "bot",
                       bound_north: str = "ymin", bound_south: str = "ymax",
                       trackingdirection = "forward"): 
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
        # Direction of calculating flow along pathlines (in modpath)
        self.trackingdirection = trackingdirection


        # Create output directories
        # Destination root
        self.dstroot = self.workspace + '\\results'

        # Create workspace(s)
        if not os.path.exists(self.dstroot):
            os.makedirs(self.dstroot)

        # Source files
        self.model_hds = os.path.join(self.workspace, self.modelname + '.hds')
        self.model_cbc = os.path.join(self.workspace, self.modelname + '.cbc')
        '''
        Initialize the AnalyticalWell object by making the dictionaries and adding the schematisation (and therefore attributes 
        of the schematisation) to the object.

        Returns
        ------
        Dictionaries (Modflow) with all parameters as an attribute of the function.

        '''
        self.schematisation = schematisation
        # Required keys
        self.required_keys = ["simulation_parameters","geo_parameters",
        "ibound_parameters","recharge_parameters",
        "well_parameters","concentration_boundary_parameters","point_parameters", "mesh_refinement",
        "endpoint_id"]

        if type(self.schematisation) == dict:
            self.schematisation_dict = self.schematisation
        else:
            #Make dictionaries
            # self.schematisation.make_dictionary()
            self.schematisation_dict = {}
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

    def _check_parameters(self,required_parameters):
        for req_var in required_parameters:
            value = getattr(self, req_var)
            print(req_var,value)
            if value is None:
                raise KeyError(f'Error, required variable {req_var} is not defined.')

    def _check_packages(self, schematisation: dict,
                            package_list: list = ["ibound"]):
        ''' Determine whether required packages and corresponding dictionaries
            occur in the schematisation dictionary 'schematisation'.'''

        for iPackage in package_list:
            if schematisation.get(iPackage + '_parameters', default=None) is None: 
                raise KeyError(f'Error, required package parameters for {iPackage} are not defined.')
          

    ## Tijdelijk overslaan (*codeSteven_20211021.xlsx)
    def _check_init_schematisation(self, required_keys,required_parameters):
        '''check the parameters that we need for the individual aquifer types are not NONE aka set by the user'''
        
        # Check schematisation dictionary
        self._check_schematisation(required_keys)
        # Check required parameters
        self._check_parameters(required_parameters)

    def _check_schematisation_input(self):
        ''' Modflow & modpath calculation using subsurface schematisation type 'phreatic'. '''
        
        # list active packages
        active_packages = ["DIS","BAS","LPF","OC","PCG"]
        # Parameter requirement
        package_parms = {"BAS": "ibound"}

        # Required keys in self.schematisation
        required_keys = ["simulation_parameters","geo_parameters",
        "ibound_parameters","recharge_parameters","well_parameters",
        "point_parameters", "mesh_refinement",
        "endpoint_id"]

        # Required parameters (to run model)
        required_parameters = []

        # Check input for schematisation type
        self._check_init_schematisation(required_keys = required_keys, required_parameters = required_parameters)  


    def _update_property(self, property, value):
        setattr(self, property, value)
        ''' set attribute value to new object attribute/property'''

    def extract_vadose_zone_parameters(self,schematisation: dict,
                                        dict_key = "geo_parameters",
                                        parameter: str = "vadose",
                                        remove_keys = True):
        ''' Extract vadose_zone interval from geo_parameters dictionary 'dict_key',
        if 'parameter' value in the 'dict_key' (e.g. "vadose") is True.
        remove_keys: --> remove the key from schematisation[dict_key] if
        it is a vadose_zone (T/F). '''
        vadose_parameters = {}
        for iDict_sub in schematisation[dict_key]:
            if parameter in schematisation[dict_key][iDict_sub]:
                if schematisation[dict_key][iDict_sub][parameter] == True:  # is geological layer a vadose zone? T/F
                    vadose_parameters[iDict_sub] = schematisation[dict_key][iDict_sub]

        # Add vadose_parameters as key for schematisation dictionary
        schematisation["vadose_parameters"] = vadose_parameters

        # Remove the keys added as vadose_parameters (vadose_zone)
        if remove_keys:
            for iKey in vadose_parameters:
                schematisation[dict_key].pop(iKey)

        return schematisation

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
                    n_ref = max(1,schematisation[iDict][iDict_sub][n_refinement])
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
            top = self.top
        # CHECK FOR ERRORS IF KEYWORD "bot" is not given
        try: 
            bot = schematisation[iDict][iDict_sub][self.bound_bot]
        except KeyError:
            print("set bottom of", iDict, iDict_sub, "to model bottom.")
            bot = min(self.bot)

        if not self.model_type == "axisymmetric":
            try:
                north = schematisation[iDict][iDict_sub][self.bound_north]
                south = schematisation[iDict][iDict_sub][self.bound_south]
                # Determine row indices
                rowidx_min = int(np.argwhere((self.ymid < north) & (self.ymid >= south))[0])
                rowidx_max = int(np.argwhere((self.ymid < north) & (self.ymid >= south))[-1]) + 1
            except KeyError as e:
                print(e,f"missing {iDict} {iDict_sub}. Continue")
        else:
            north,south,rowidx_min,rowidx_max = None, None,0,1
        try:
            # Determine layer indices
            layidx_min = int(np.argwhere((top >= self.zmid) & (bot < self.zmid))[0])
            layidx_max = int(np.argwhere((top >= self.zmid) & (bot < self.zmid))[-1]) + 1
        
        except IndexError as e:
            print(e, iDict,iDict_sub,top,bot, "(top,bot)")
            print("Set layidx_min and layidx_max to None.")
            layidx_min, layidx_max = None, None

        try:
            # Determine column indices
            colidx_min = int(np.argwhere((left <= self.xmid) & (right > self.xmid))[0])
            colidx_max = int(np.argwhere((left <= self.xmid) & (right > self.xmid))[-1]) + 1
            # np.where((self.xmid < right) & (self.xmid > left))    
        except IndexError as e:
            print(e, iDict,iDict_sub,left,right, "(left,right)")
            print("Set colidx_min and colidx_max to None.")
            colidx_min, colidx_max = None, None

        return layidx_min,layidx_max,rowidx_min,rowidx_max,colidx_min,colidx_max

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

                # if 'mesh_refinement' not in iDict_sub:

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
                        colidx_min: colidx_max][self.material[layidx_min: layidx_max,\
                        rowidx_min: rowidx_max,\
                        colidx_min: colidx_max] == None] = iDict_sub
                    
        if self.model_type in ["axisymmetric","2D"]:
                # In 2D model or axisymmetric models an 
                # inactive row is added to be able to run Modpath successfully.
                self.material[:,1,:] = self.material[:,0,:]


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

      

    def assign_wellloc(self,schematisation: dict,
                        dict_key: str = "well_parameters",
                        well_names: list or str or None = "None",
                        discharge_parameter = "well_discharge"):

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
            # Daily flux [m3/d]   (negative value)  
            Qwell_day[iWell] = schematisation[dict_key][iWell][discharge_parameter]
            # Calculate boundary indices
            layidx_min,layidx_max,\
                rowidx_min,rowidx_max,\
                colidx_min,colidx_max = self.cell_bounds(schematisation,
                                        dict_key = dict_key,
                                        dict_subkey = iWell,
                                        model_type = self.model_type)

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

            # stress period data for well package
            for iLay in range(layidx_min,layidx_max):
                for iRow in range(rowidx_min,rowidx_max):
                   for iCol in range(colidx_min, colidx_max):
                        spd_wel[0].append([iLay, iRow, iCol, Qwell_day[iWell] * \
                                          (self.hk[iLay,iRow,iCol] * self.delv[iLay]) / KD_well[iWell]])
        
        return well_names,well_loc,KD_well, spd_wel, Qwell_day

    ####################
    ### Fill modules ###
    def create_modflow_input(self, **kwargs):
        """ Generate flopy files."""
            
        # Use dictionary keys from schematisation
        dict_keys = ["geo_parameters","recharge_parameters","ibound_parameters",
                      "well_parameters","mesh_refinement"]
        self.make_discretisation(schematisation = self.schematisation_dict, dict_keys = dict_keys,
                            model_type = self.model_type)

        # Set ibound grid and starting head
        dict_keys_ibound = ["ibound_parameters"]

        self.ibound = np.ones((self.nlay,self.nrow,self.ncol), dtype = 'int')
        self.strt = np.zeros((self.nlay,self.nrow,self.ncol), dtype = 'float')

        # Relevant dictionary keys
        self.ibound = self.fill_grid(schematisation = self.schematisation_dict,
                                    dict_keys = dict_keys_ibound,
                                    parameter = "ibound",
                                    grid = self.ibound,
                                    dtype = 'int')
        self.strt = self.fill_grid(schematisation = self.schematisation_dict,
                                    dict_keys = dict_keys_ibound,
                                    parameter = "head",
                                    grid = self.strt,
                                    dtype = 'int')      

        # Set ibound of additional rows equal to zero
        self.ibound[:,1:,:] = 0

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
        self.geoparm_names = {"moisture_content": [["geo_parameters"],"float",0.35],
                        "hk": [["geo_parameters"],"float",999.],
                        "vani": [["geo_parameters"],"float",999.],
                        "porosity": [["geo_parameters"],"float",1.],
                        "solid_density": [["geo_parameters"],"float",2.5],
                        "fraction_organic_carbon": [["geo_parameters"],"float",0.0],
                        "redox": [["geo_parameters"],"object",0],
                        "dissolved_organic_carbon": [["geo_parameters"],"float",0],
                        "pH": [["geo_parameters"],"float",7],
                        "temperature": [["geo_parameters"],"float",11]
                        }
                      

        for iParm, dict_keys in self.geoparm_names.items():
            # Temporary value grid
            unitgrid = np.zeros((self.nlay,self.nrow,self.ncol), dtype = dict_keys[1]) + dict_keys[2]
            # Fill grid with new values
            grid = self.fill_grid(schematisation = self.schematisation_dict,
                                dict_keys = dict_keys[0],
                                parameter = iParm,
                                grid = unitgrid,
                                dtype = dict_keys[1])
            self._update_property(property = iParm, value = grid)

        # Assign material grid
        self.assign_material(schematisation = self.schematisation_dict,
                            dict_keys = ["geo_parameters","ibound_parameters","well_parameters"])

        # Create (uncorrected) array for kv ("vka"), using "kh" and "vani" (vertical anisotropy)
        for iLay in range(self.nlay):
            for iRow in range(self.nrow):
                for iCol in range(self.ncol):
                    # Check if material occurs in "geo_parameters"
                    if self.material[iLay,iRow,iCol] not in self.schematisation_dict["geo_parameters"].keys():
                        self.hk[iLay,iRow,iCol] = 999.
                        self.vani[iLay,iRow,iCol] = 999.
                        self.redox[iLay,iRow,iCol] = "anoxic"

        # Vertical conductivity [m/d]
        self.vka = self.hk / self.vani
        # and for 'storativity' ("ss": specific storage)
        self.ss = np.ones((self.nlay,self.nrow,self.ncol), dtype = 'float') * 1.E-6
        # Uncorrected porosity parameter (for removal calculation)
        self.prsity_uncorr = copy.deepcopy(self.porosity)

        # Axisymmetric flow properties
        axisym_parms = ["hk","vka","ss","porosity"]
        if self.model_type == "axisymmetric":
            for iParm in axisym_parms:
                grid_uncorr = getattr(self,iParm)
                grid_axi = self.axisym_correction(grid = grid_uncorr)
                # Update attribute
                self._update_property(property = iParm, value = grid_axi)

        ### Outstanding issue: how to deal with multiple recharge sources 
        # Create recharge package
        rech_parmnames = {"recharge": [["recharge_parameters"],"float"]}
        for iParm, dict_keys in rech_parmnames.items():
            # Temporary value
            grid_temp = np.zeros((self.nlay,self.nrow,self.ncol), dtype = 'float')
            grid = self.fill_grid(schematisation = self.schematisation_dict,
                            dict_keys = dict_keys[0],
                            parameter = iParm,
                            grid = grid_temp,
                            dtype = dict_keys[1])
            
            if self.model_type == "axisymmetric":
                rech_grid = self.axisym_correction(grid = grid)[0,:,:]
            else:
                rech_grid = grid[0,:,:]

            # Update attribute (recharge)
            self._update_property(property = iParm, value = rech_grid)

        # Well input
        # !!! Obtain node numbers of well locations (use indices) !!!
        # leakage discharge from well
        well_names = [iWell for iWell in self.schematisation_dict["well_parameters"] if \
                        "well_discharge" in self.schematisation_dict["well_parameters"][iWell].keys()]
       
        self.well_names,\
            self.well_loc,\
                self.KD_well,\
                    self.spd_wel,\
                        self.Qwell_day = self.assign_wellloc(schematisation = self.schematisation_dict,
                                                        dict_key = "well_parameters",
                                                        well_names = None,
                                                        discharge_parameter = "well_discharge")



    def create_mfobject(self, mf_exe = 'mf2005.exe', fname_nam = None):
        ''' create mf-object (if fname_nam equals None).
            Otherwise: load mf-model using nam-file 'fname_nam'.
            '''
        self.mf_exe = mf_exe
        if fname_nam is None:
            # Create new Modflow object
            self.mf = flopy.modflow.Modflow(modelname = self.modelname, exe_name= self.mf_exe, model_ws=self.workspace)
        else:
            # Load Modflow model
            self.mf = flopy.modflow.Modflow.load(fname_nam)

    def create_dis(self, per_nr = 0):
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

    def create_bas(self):
        ''' Add basic Package to the MODFLOW model '''
        self.bas = flopy.modflow.ModflowBas(self.mf, ibound = self.ibound, strt = self.strt)
        
    def create_lpf(self):
        ''' Add lpf Package to the MODFLOW model '''
        self.lpf = flopy.modflow.ModflowLpf(self.mf, layavg = 1, ipakcb = self.iu_cbc, hk=self.hk, vka=self.vka, 
                                            ss = self.ss, storagecoefficient = True) 
        # layavg = 1 (--> logarithmic mean); storagecoefficient = True (means: storativity)
        
    def create_pcg(self, hclose=1e-4, rclose = 0.001):
        ''' Add PCG Package to the MODFLOW model '''
        
        self.pcg = flopy.modflow.ModflowPcg(self.mf, hclose=hclose, rclose = rclose)
            
    def create_oc(self, spd_oc = {(0, 0): ['save head', 'save budget']},
                per_nr = 0):
        ''' Add output control Package to the 
        MODFLOW model. '''    

        # Extensions of  output files
        extension = ['oc', 'hds', 'ddn', 'cbc', 'ibo']
        # Unit numbers (of modflow packages and output)
        unitnumber = [14, 51, 52, self.iu_cbc, 0]
        filenames = []
        for iExt in extension:
            filenames.append(self.modelname + "." + iExt)

        self.oc = flopy.modflow.ModflowOc(self.mf, stress_period_data= spd_oc,
                                          extension = extension,  # Default extensions
                                          filenames = filenames,
                                          compact = True #,cboufm='(20i5)'
                                          )                                            
        # compact option is a requirement for using Modpath

    def create_rch(self, rech = None):
        ''' Add recharge Package to the MODFLOW model. '''     
                                           
        # Stress period data dict 
        if rech is not None:                       
            self.rech = rech

        self.rch = flopy.modflow.ModflowRch(self.mf, ipakcb = self.iu_cbc, #102, 
                                     rech = self.rech, nrchop = 1)          

    def create_wel(self, spd_wel = None):
        ''' Add well Package to the MODFLOW model. '''     
                                           
        # Stress period data dict                        
        self.spd_wel = spd_wel

        self.wel = flopy.modflow.ModflowWel(self.mf, ipakcb = self.iu_cbc,
                                     stress_period_data= self.spd_wel)
         # If ipakcb != 0: cell budget data is being saved.

    def create_modflow_packages(self, **kwargs):
        ''' Create modflow packages used in the model run.'''

        # Start model run
        # Load modflow packages
        for iper in range(self.nper):
            # If the starting heads have been changed during the first model run:
            if iper != 0:

                try:
                    # Load head data
                    _, self.strt = self.read_hdsobj(fname = self.model_hds,time = -1)
                except Exception as e:
                    print (e,"Error loading strt_fw data\nstrt_fw set to '0'.")
                    self.strt = 0.
            
            # Open modflow object
            self.create_mfobject(mf_exe = 'mf2005.exe')
            # Load packages
            self.create_dis(per_nr = iper)   
            self.create_bas()   
            self.create_lpf()   
            self.create_pcg()   
            # Output control stress period data
            self.spd_oc = {(0, 0): ['save head', 'save budget']}
            self.create_oc(spd_oc = self.spd_oc)
            if "well_parameters" in kwargs.keys():
                if len(kwargs["well_parameters"]) > 0:
                    try:                        
                        self.create_wel(spd_wel = self.spd_wel)
                    except Exception as e:
                        print(e, "error loading well package.")
            if "recharge_parameters" in kwargs.keys():
                if len(kwargs["recharge_parameters"]) > 0:
                    try:
                        self.create_rch(rech = self.recharge)
                    except Exception as e:
                        print(e, "no recharge assigned.")

    def generate_modflow_files(self):
        ''' Write package data MODFLOW model. '''
        self.mf.write_input()

        # Try to delete the previous output files, to prevent accidental use of older files
        try:  
            os.remove(self.model_hds)
        except FileNotFoundError:
            pass
        try:  
            os.remove(self.model_cbc)
        except FileNotFoundError:
            pass

    def run_modflowmod(self):         
        ''' Run modflow model '''
        self.success_mf, _ = self.mf.run_model(silent=False)

    def mfmodelrun(self):
        ''' Main model run code. '''
        print ("Run model:",self.workspace, self.modelname +"\n")
        # Create modflow packages
        self.create_modflow_packages(**self.schematisation_dict)
        # Generate modflow files
        self.generate_modflow_files()

        # Run modflow model
        try:
            self.run_modflowmod()
        except Exception as e:
            self.success_mf = False
            print(e, self.success_mf)
        # print(self.success_mf, self.buff)

        # Model run completed succesfully
        print("Model run", self.workspace, self.modelname, "completed without errors:", self.success_mf)

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

    def _create_diffuse_particles(self, recharge_parameters = None,
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
        else:
            recharge_parameters = self.schematisation_dict.get(recharge_parameters)
        ## Particle group data ##
        if recharge_parameters is None:
            pgroups = []
        else:
            pgroups = list(recharge_parameters.keys())

        # xmin and xmax per pg
        if not hasattr(self, "pg_xmin"):
            self.pg_xmin = {}
        if not hasattr(self, "pg_xmax"):
            self.pg_xmax = {}
        # ymin and ymax per pg
        if not hasattr(self, "pg_ymin"):
            self.pg_ymin = {}
        if not hasattr(self, "pg_ymax"):
            self.pg_ymax = {}
        # zmin and zmax per pg
        if not hasattr(self, "pg_zmin"):
            self.pg_zmin = {}
        if not hasattr(self, "pg_zmax"):
            self.pg_zmax = {}

        # particle group filenames
        if not hasattr(self, "pg_filenames"):
            self.pg_filenames = {}
        for iPG in pgroups:
            self.pg_filenames[iPG] = iPG + ".sloc"

    	# flowline_type: point_source (or diffuse_source)
        if not hasattr(self, "flowline_type"):
            self.flowline_type = {}
        if not hasattr(self, "point_discharge"):
            self.point_discharge = {}      

        # particle starting locations [(lay,row,col),(k,i,j),...]
        if not hasattr(self, "part_locs"):
            self.part_locs = {}  
        # Relative location within grid cells (per particle group)
        if not hasattr(self, "localx"):
            self.localx = {}  
        if not hasattr(self, "localy"):
            self.localy = {}  
        if not hasattr(self, "localz"):
            self.localz = {}  

        # Particle id [part group, particle id] - zero based integers
        if not hasattr(self, "pids"):
            print("Create a new particle id dataset dict.")
            self.pids = {}
        if not hasattr(self, "pg"):
            print("Create a new particle group dataset dict.")
            self.pg = {}
        # Particle group nodes
        if not hasattr(self, "pg_nodes"):
            print("Create a new particle group nodes dataset dict.")
            self.pg_nodes = {}

        # Particle data objects
        if not hasattr(self, "pd"):
            print("Create a new particle dataset dict.")
            self.pd = {}


        # Particle counter
        if not hasattr(self, "pcount"):
            print("Create a new particle counter.")
            self.pcount = -1
            
        for iPG in pgroups:

            if len(iPG) == 0:
                # Particle group not entering via recharge
                continue

            # xmin
            self.pg_xmin[iPG] = recharge_parameters.get(iPG).get("xmin")
            # xmax
            self.pg_xmax[iPG] = recharge_parameters.get(iPG).get("xmax")

            # Only if both ymin and ymax are given
            if ("ymin" in recharge_parameters.get(iPG)) & \
                ("ymax" in recharge_parameters.get(iPG)):

                # ymin
                self.pg_ymin[iPG] = recharge_parameters.get(iPG).get("ymin")
                # ymax
                self.pg_ymax[iPG] = recharge_parameters.get(iPG).get("ymax")
            else:
                # ymin
                self.pg_ymin[iPG] = self.ymid[0]
                # ymax
                self.pg_ymax[iPG] = self.ymid[0]



            # zmin
            if "zmin" in recharge_parameters.get(iPG):
                self.pg_zmin[iPG] = recharge_parameters.get(iPG).get("zmin")
            else:
                if "zmax" in recharge_parameters.get(iPG):
                    self.pg_zmin[iPG] = recharge_parameters.get(iPG).get("zmax")
                else:
                    self.pg_zmin[iPG] = self.top
            # zmax
            if "zmax" in recharge_parameters.get(iPG):
                self.pg_zmax[iPG] = recharge_parameters.get(iPG).get("zmax")
            else:
                if "zmin" in recharge_parameters.get(iPG):
                    self.pg_zmax[iPG] = recharge_parameters.get(iPG).get("zmin")
                else:
                    self.pg_zmax[iPG] = self.pg_zmin[iPG]

            try:
                # Determine column indices
                colidx_min = int(np.argwhere((self.xmid + 0.5 * self.delr >= self.pg_xmin[iPG]) & \
                    (self.xmid - 0.5 * self.delr <= self.pg_xmax[iPG]))[0])
                colidx_max = int(np.argwhere((self.xmid + 0.5 * self.delr >= self.pg_xmin[iPG]) & \
                    (self.xmid - 0.5 * self.delr <= self.pg_xmax[iPG]))[-1])
                # np.where((self.xmid < right) & (self.xmid > left))    
            except IndexError as e:
                print(e,"Set colidx_min and colidx_max to None.")
                colidx_min, colidx_max = None, None            

            try:
                # Determine row indices
                rowidx_min = int(np.argwhere((self.ymid + 0.5 * self.delc >= self.pg_ymin[iPG]) & \
                    (self.ymid - 0.5 * self.delc <= self.pg_ymax[iPG]))[0])
                rowidx_max = int(np.argwhere((self.ymid + 0.5 * self.delc >= self.pg_ymin[iPG]) & \
                    (self.ymid - 0.5 * self.delc <= self.pg_ymax[iPG]))[-1])
                # np.where((self.xmid < right) & (self.xmid > left))    
            except IndexError as e:
                print(e,"Set rowidx_min and rowidx_max to None.")
                rowidx_min, rowidx_max = 0, 0            

            try:
                # Determine layer indices
                layidx_min = int(np.argwhere((self.zmid + 0.5 * self.delv >= self.pg_zmin[iPG]) & \
                    (self.zmid - 0.5 * self.delv <= self.pg_zmax[iPG]))[0]) # shallow layer
                layidx_max = int(np.argwhere((self.zmid + 0.5 * self.delv >= self.pg_zmin[iPG]) & \
                    (self.zmid - 0.5 * self.delv <= self.pg_zmax[iPG]))[-1])  # deepest layer
                # np.where((self.xmid < right) & (self.xmid > left))    
            except IndexError as e:
                print(e,"Set layidx_min and layidx_max to 0.")
                layidx_min, layidx_max = 0, 0            

            # Particles ids (use counter)
            self.pids[iPG] = []
            # Particle location list per particle group
            self.part_locs[iPG] = []
            # Relative location within grid cells (per particle group)
            self.localx[iPG] = []
            self.localy[iPG] = []
            self.localz[iPG] = []
            # Particle group nodes
            self.pg_nodes[iPG] = []

            for iRow in range(rowidx_min,rowidx_max+1):
                for iLay in range(layidx_min,layidx_max+1):
                    for iCol in range(colidx_min,colidx_max+1):
                    
                        # Particle row
                        p_row = iRow
                        # Particle layer
                        p_lay = iLay
                        # Particle column
                        p_col = iCol

                        # check for active cell
                        if self.ibound[p_lay,p_row,p_col] == 0:
                            continue
                        
                        # locations for reading output in modpath model
                        self.pg_nodes[iPG].append((p_lay,p_row,p_col))
                        for iPart_cell in range(nparticles_cell):
                            # Add particle locations (lay,row,col)
                            self.part_locs[iPG].append((p_lay,p_row,p_col))                    
                            # Relative location of the particles in the cells
                            self.localx[iPG].append((iPart_cell + 0.5)/(float(nparticles_cell)))
                            self.localy[iPG].append(localy)
                            self.localz[iPG].append(localz)
                            # particle count
                            self.pcount += 1 
                            self.pids[iPG].append(self.pcount)

                            # flowline_type: diffuse_source (or point_source)
                            self.flowline_type[self.pcount] = 'diffuse_source'  


            # Particle distribution package - particle allocation
            #modpath.mp7particledata.Part...
            self.pd[iPG] = flopy.modpath.ParticleData(partlocs = self.part_locs[iPG], structured=True,
                                                drape=drape, localx= self.localx[iPG], 
                                                localy= self.localy[iPG], localz= self.localz[iPG],
                                                timeoffset = timeoffset, 
                                                particleids = self.pids[iPG])


            # particle group object
            # modpath.mp7particlegroup.Part.......
            self.pg[iPG] = flopy.modpath.ParticleGroup(particlegroupname=iPG,
                                                                filename=self.pg_filenames[iPG],
                                                                releasedata=releasedata,
                                                                particledata=self.pd[iPG])
            ''' ParticleGroup class to create MODPATH 7 particle group data for
                location input style 1.  '''
            self.trackingdirection = trackingdirection

    def _create_point_particles(self, point_parameters = None,
                                                   localx = 0.5, 
                                                   localy = 0.5, 
                                                   localz = 0.5,
                                                   timeoffset=0.0, drape=0,
                                                   trackingdirection = 'forward',
                                                   releasedata=0.0):
        ''' Class to create the most basic particle data type (starting location
        input style 1). Input style 1 is the most general input style and provides
        the highest flexibility in customizing starting locations, see flopy docs.
        ###########################################################################
        Create points using dictionary 'point_parameters' with start position(s) 'x_start', 'y_start' and 'z_start'.
        '''

        
        if point_parameters is None:
            point_parameters = self.schematisation_dict.get('point_parameters')
        else:
            point_parameters = self.schematisation_dict.get(point_parameters)

        ## Particle group data ##
        if point_parameters is None:
            pgroups = []
        else:
            pgroups = list(point_parameters.keys())

        # xmin and xmax per pg
        if not hasattr(self, "x_start_particle"):
            self.x_start_particle = {}
        if not hasattr(self, "y_start_particle"):
            self.y_start_particle = {}
        if not hasattr(self, "z_start_particle"):
            self.z_start_particle = {}            

        # particle group filenames
        if not hasattr(self, "pg_filenames"):
            self.pg_filenames = {}
        for iPG in pgroups:
            self.pg_filenames[iPG] = iPG + ".sloc"

    	# flowline_type: point_source (or diffuse_source)
        if not hasattr(self, "flowline_type"):
            self.flowline_type = {}     
        if not hasattr(self, "point_discharge"):
            self.point_discharge = {}                

        # particle starting locations [(lay,row,col),(k,i,j),...]
        if not hasattr(self, "part_locs"):
            self.part_locs = {}  
        # Relative location within grid cells (per particle group)
        if not hasattr(self, "localx"):
            self.localx = {}  
        if not hasattr(self, "localy"):
            self.localy = {}  
        if not hasattr(self, "localz"):
            self.localz = {}  

        # Particle id [part group, particle id] - zero based integers
        if not hasattr(self, "pids"):
            print("Create a new particle id dataset dict.")
            self.pids = {}
        if not hasattr(self, "pg"):
            print("Create a new particle group dataset dict.")
            self.pg = {}
        # Particle group nodes
        if not hasattr(self, "pg_nodes"):
            print("Create a new particle group nodes dataset dict.")
            self.pg_nodes = {}

        # Particle data objects
        if not hasattr(self, "pd"):
            print("Create a new particle dataset dict.")
            self.pd = {}


        # Particle counter
        if not hasattr(self, "pcount"):
            print("Create a new particle counter.")
            self.pcount = -1
            
        for iPG in pgroups:

            if len(iPG) == 0:
                # Particle group not entering via point source
                continue

            # x_start
            if "x_start" in point_parameters.get(iPG):
                self.x_start_particle[iPG] = point_parameters.get(iPG).get("x_start")
            else:
                self.x_start_particle[iPG] = self.xmid[0]
            # y_start
            if "y_start" in point_parameters.get(iPG):
                self.y_start_particle[iPG] = point_parameters.get(iPG).get("y_start")
            else:
                self.y_start_particle[iPG] = self.ymid[0]
            # z_start
            if "z_start" in point_parameters.get(iPG):
                self.z_start_particle[iPG] = point_parameters.get(iPG).get("z_start")
            else:
                self.z_start_particle[iPG] = self.zmid[0]

            try:
                # Determine particle column idx
                p_col = np.argmin(abs(self.xmid - self.x_start_particle[iPG]))
            except IndexError as e:
                print(e,"Set p_col to 0")
                p_col = 0  

            try:   
                # Determine particle row idx  
                p_row = np.argmin(abs(self.ymid - self.y_start_particle[iPG]))
            except IndexError as e:
                print(e,"Set p_row to 0")
                p_row = 0    

            try:
                # Determine particle layer idx
                p_lay = np.argmin(abs(self.zmid - self.z_start_particle[iPG]))
            except IndexError as e:
                print(e,"Set p_lay to 0")
                p_lay = 0  

            # Particles ids (use counter)
            self.pids[iPG] = []
            # Particle location list per particle group
            self.part_locs[iPG] = []
            # Relative location within grid cells (per particle group)
            self.localx[iPG] = []
            self.localy[iPG] = []
            self.localz[iPG] = []
            # Particle group nodes
            self.pg_nodes[iPG] = []

            # Cell boundary locations
            left, right = self.xmid[p_col] - 0.5 * self.delr[p_col], self.xmid[p_col] + 0.5 * self.delr[p_col]
            back, front = self.ymid[p_row] - 0.5 * self.delc[p_row], self.ymid[p_row] + 0.5 * self.delc[p_row]
            bot, top = self.zmid[p_lay] - 0.5 * self.delv[p_lay], self.zmid[p_lay] + 0.5 * self.delv[p_lay]

            # Calculate localx, localy, localz using cell boundary locations
            localx = (self.x_start_particle[iPG] - left) / (right - left)
            localy = (self.y_start_particle[iPG] - back) / (front - back)
            localz = (self.z_start_particle[iPG] - bot) / (top - bot)

            if self.ibound[p_lay,p_row,p_col] != 0:
                # locations for reading output in modpath model
                self.pg_nodes[iPG].append((p_lay,p_row,p_col))
                # Add particle locations (lay,row,col)
                self.part_locs[iPG].append((p_lay,p_row,p_col))                    
                # Relative location of the particles in the cells
                self.localx[iPG].append(localx)
                self.localy[iPG].append(localy)
                self.localz[iPG].append(localz)
                # particle count
                self.pcount += 1 
                self.pids[iPG].append(self.pcount)

                # flowline_type: diffuse_source or point_source
                self.flowline_type[self.pcount] = 'point_source'
                # Add point_discharge
                self.point_discharge[self.pcount] = point_parameters.get(iPG).get("discharge") 


                # Particle distribution package - particle allocation
                #modpath.mp7particledata.Part...
                self.pd[iPG] = flopy.modpath.ParticleData(partlocs = self.part_locs[iPG], structured=True,
                                                    drape=drape, localx= self.localx[iPG], 
                                                    localy= self.localy[iPG], localz= self.localz[iPG],
                                                    timeoffset = timeoffset, 
                                                    particleids = self.pids[iPG])


                # particle group object
                # modpath.mp7particlegroup.Part.......
                self.pg[iPG] = flopy.modpath.ParticleGroup(particlegroupname=iPG,
                                                                    filename= self.pg_filenames[iPG],
                                                                    releasedata=releasedata,
                                                                    particledata=self.pd[iPG])
                ''' ParticleGroup class to create MODPATH 7 particle group data for
                    location input style 1.  '''
            self.trackingdirection = trackingdirection


    ####################################
    ### ModPath 7 input and analyses ###
    ####################################
    ####### Fill Modpath modules #######
    def create_MP7object(self, mp_exe = 'mpath7', mf_model = None,
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
        
    def create_mpbas(self):
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

    def _export_to_df(self, mppth):
        """ Makes 'df_flowline' and 'df_particle' for ModPath model simulation
  
            Parameters
            -------
            mppth: str
                File location to ModPath pathline output (*.mppth)

            Returns
            -------
            df_flowline: pandas.DataFrame
                Column 'flowline_id': Integer
                Column 'flowline_type': string
                    Described the type of contamination associated with the flowline,
                    either 'diffuse_source' or 'point_source'.
                Column 'flowline_discharge': Float
                    Discharge associated with the flowline, [m3/d].
                Column 'particle_release_day': Float
                Column 'input_concentration'
                Column 'endpoint_id': Integer
                    ID of Well (or drain) where the flowline ends.
                Column 'well_discharge': float
                Column 'substance': string
                Column 'removal_function': string

            df_particle: pandas.DataFrame
                Column 'flowline_id': int
                Column 'zone': string
                    Zone in the aquifer, ground surface 'surface', 'vadose_zone', 'shallow_aquifer' or 'target_aquifer'
                Column 'travel_time': float
                    Travel time in the respective aquifer zone given in column 'zone.
                Column 'xcoord': float
                Column 'ycoord': float
                Column 'zcoord': float
                Column 'redox': float
                    'suboxic', 'anoxic', deeply_anoxic'
                Column 'temperature': float
                    Of the respective aquifer zone.
                Column 'travel_distance': float
                Column 'porosity': float
                Column 'dissolved_organic_carbon': float
                Column 'pH': float
                Column 'fraction_organic_carbon': float
                Column 'solid_density': float
        """

        # MK:  this is something that I could imagine that you use it as an argument to the
        # function.
        #AH_todo, @MartinK, what is the advantage over using the attribute?
        # what_to_export = self.schematisation.what_to_export

        # pthline_data converted to rec.arrays (with fields appended)
        self.particle_data = {}
        # dataframes of particle data (dict)
        self.df_particle_data = {}
        # list the particle_data dataframes
        df_particle_list = []

        # Create rec.arrays for porosity, pH, T, etc. to append to particle_data
        parm_list = ["prsity_uncorr","solid_density","fraction_organic_carbon",
                    "redox","dissolved_organic_carbon", "pH","temperature","material"]

        # Create df_particle
        self.df_particle, self.df_particle_data = self.fill_df_particle(
                            particle_group = self.pg,
                            pg_nodes = self.pg_nodes,
                            parm_list = parm_list,
                            mppth = mppth)

        ''' TODO 211123: append phreatic pathline data '''

        # df_particle file name
        particle_fname = os.path.join(self.dstroot,self.schematisation_type + "_df_particle.csv")
        # Save df_particle
        self.df_particle.to_csv(particle_fname)
        
        # Create dataframe df_flowline
        self.df_flowline = self.fill_df_flowline(df_particle = self.df_particle,
                                                model_cbc = self.model_cbc)

        # df_flowline file name
        flowline_fname = os.path.join(self.dstroot,self.schematisation_type + "_df_flowline.csv")
        # Save df_flowline
        self.df_flowline.to_csv(flowline_fname)

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
                self.create_mfmodel(self.mf_namfile)
        else:
            if hasattr(self, "mf"):
                # modflow model "mf" already exists in object
                pass
            else:
                print("Load mf model using given namfile location:", mf_namfile)
                # mf_namefile
                self.mf_namfile = mf_namfile
                self.create_mfmodel(self.mf_namfile)

        
        # Load the MP7 module/object
        if mp_exe is None:
            self.mp_exe = "mpath7"
        else:
            self.mp_exe = mp_exe

        # Load mp7 object into Aximodel class
        self.create_MP7object(mp_exe = self.mp_exe, mf_model = self.mf)#,

        # Load modpath basic package
        self.create_mpbas()
        
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

    def read_hdsobj(self, fname = None, time = None):
        
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

    def read_binarycbc_flow(self, fname = None):
        ''' Read binary cell budget file (fname). 
            This is modflow output.

            return frf, flf, fff.'''
            
        cbcobj = bf.CellBudgetFile(fname)
        print(cbcobj.list_records())
        try:
            frf = cbcobj.get_data(text='FLOW RIGHT FACE')[0]
            flf = cbcobj.get_data(text='FLOW LOWER FACE')[0]
            fff = cbcobj.get_data(text='FLOW FRONT FACE')[0]
        finally:
            cbcobj.close()
        
        return frf, flf, fff 

    def get_node_ID(self,locs):
        ''' Obtain/return model node index (int) belonging to
            layer 'iLay', row 'iRow' and column 'iCol'.
            for point locations 'locs'. '''
        nodes = []
        for iLoc in locs:
            iLay,iRow,iCol = iLoc[0], iLoc[1], iLoc[2]
            nodes.append(iLay * self.nrow * self.ncol + iRow * self.ncol + iCol)
        return nodes

    def xyz_to_layrowcol(self,xyz_point, round_digits = 3):
        ''' obtain lay, row, col index of an xyz_point list [x,y,z]. '''

        # Round XYZ to 4 digits
        xyz_point = [round(iVal,round_digits) for iVal in xyz_point]

        # min and max arrays in X-direction
        xmin_arr = np.round(self.xmid-0.5*self.delr,round_digits)
        xmax_arr = np.round(self.xmid+0.5*self.delr,round_digits)
        # min and max arrays in Y-direction
        ymin_arr = np.round(self.ymid-0.5*self.delc,round_digits)
        ymax_arr = np.round(self.ymid+0.5*self.delc,round_digits)
        # min and max arrays in Z-direction
        zmin_arr = np.round(self.zmid-0.5*self.delv,round_digits)
        zmax_arr = np.round(self.zmid+0.5*self.delv,round_digits)

        # node layer, row, column
        node_lay = np.argwhere((xyz_point[2] > zmin_arr) & \
                    (xyz_point[2] <= zmax_arr))[0][0]
        node_row = np.argwhere((xyz_point[1] > ymin_arr) & \
                    (xyz_point[1] <= ymax_arr))[0][0]
        node_col = np.argwhere((xyz_point[0] > xmin_arr) & \
                    (xyz_point[0] <= xmax_arr))[0][0]

        # Node index (lay,row,col)
        node_idx = (node_lay,node_row,node_col)

        return node_idx

    def get_node_indices(self,xyz_nodes, particle_list: list or None = None):
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
            node_indices[iPart] = []
            for iNode in range(nr_nodes):
                # list lay, row, col of nodes and store them in 'node_indices' dict
                node_idx = self.xyz_to_layrowcol(xyz_point = xyz_nodes[iPart][iNode].tolist(), round_digits = 2)
                node_indices[iPart].append(node_idx)

        return node_indices

    def calc_flux_cell(self, frf,flf,fff, loc):
        ''' Calculate the total volume flux along the cell boundary, using
        frf, flf and fff.
        With: 
            frf: np.array
                flux right face (= positive to the right)
            flf: np.array
                flux lower face (= positive in downward direction)
            fff: np.array
                flux front face (= positive in 'southward' direction)
            loc: tuple
                cell location (lay,row,col)
        '''

        # fluxes along x-direction
        flux_east = max(0,frf[loc[0],loc[1],loc[2]])
        if loc[2]-1 >= 0:
            flux_west = max(0,-frf[loc[0],loc[1],loc[2]-1])
        else:
            flux_west = 0.
        
        # fluxes along y-direction
        flux_south = max(0,fff[loc[0],loc[1],loc[2]])
        if loc[1]-1 >= 0:
            flux_north = max(0,-fff[loc[0],loc[1]-1,loc[2]])
        else:
            flux_north = 0.
            
        # fluxes along z-direction
        flux_lower = max(0, flf[loc[0],loc[1],loc[2]])
        if loc[0]-1 >= 0:
            flux_upper = max(0,-flf[loc[0],loc[1],loc[2]])
        else:
            flux_upper = 0.
        
        # Total flux in all directions
        flux_total = flux_east + flux_west + flux_north + flux_south + flux_upper + flux_lower

        return flux_total

    def read_pathlinedata(self,fpth, nodes):
        ''' read pathlinedata from file fpth (extension: '*.mppth'), 
            given the particle release node index 'nodes' obtained from 
            a tuple or list of tuples (iLay,iRow,iCol).
            Cell node indices can be obtained using the method get_node_ID((iLay,iRow,iCol))
            
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
        # Round rec.arrays to 3 decimals
        for idx,iPart in enumerate(pth_data):
            pth_data[idx]["x"] = pth_data[idx]["x"].round(3)
            pth_data[idx]["y"] = pth_data[idx]["y"].round(3)
            pth_data[idx]["z"] = pth_data[idx]["z"].round(3)

            # # Test unique
            # test_unique = np.unique(pth_data[idx],return_index=True, return_inverse=True, axis = 0)

        time_nodes, xyz_nodes, txyz, dist, tdiff = {}, {}, {}, {}, {}
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
            time_nodes[idx] = txyz[idx][:,0]
            # xyz data
            xyz_nodes[idx] = txyz[idx][:,1:]

            # Determine number of remaining nodes
            n_nodes = min(xyz_nodes[idx].shape[0],time_nodes[idx].shape[0])
            if xyz_nodes[idx].shape[0] != time_nodes[idx].shape[0]:
                print(xyz_nodes[idx].shape[0],time_nodes[idx].shape[0])
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
                tdiff[idx][iNode-1] = (time_nodes[idx][iNode] - time_nodes[idx][iNode - 1])
                
            # Total distance covered per particle
            dist_tot[idx] = dist[idx].sum()
            # Total time covered per particle
            time_tot[idx] = time_nodes[idx][-1]
                        
        return time_nodes, xyz_nodes, dist, tdiff, dist_tot, time_tot, pth_data

    # Make df_particle
    def fill_df_particle(self, particle_group,
                        pg_nodes,
                        parm_list: list = list(),
                        mppth = None
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

        '''

        if mppth is None:
            # Use default
            mppth = os.path.join(self.workspace, self.modelname + '_mp.mppth')

        # if particle_data is None:
        # Temporary storage (dict) of particle_data
        particle_data = {}
        # dataframes of particle data (dict)
        df_particle_data = {}
        # list the particle_data dataframes
        df_particle_list = []

        for iPG in particle_group:  # use endpoint_id dict or list

            # Node indices to retrieve particle data
            nodes = self.get_node_ID(pg_nodes.get(iPG))

            # Nodes to retrieve
            for id_,iNode in enumerate(nodes):
                # print(id_,iNode)
                # try:


                # nodes_rel[ = flopy.utils.ra_slice(m.wel.stress_period_data[0], ['k', 'i', 'j'])
                #     nodes_well = prf.get_nodes(locs = wel_locs, nrow = nrow, ncol = ncol) # m.dis.get_node(wel_locs.tolist())

                # Read pathline data
                time_nodes, xyz_nodes, dist_data,  \
                time_diff, dist_tot,   \
                time_tot, pth_data =  \
                            self.read_pathlinedata(fpth = mppth,
                                                nodes = iNode)
                
                '''                             
                xyz_points, \    # XYZ data 
                dist_data,  \    # Distance array between nodes
                time_diff,  \    # Save flow duration (time_diff) of pathlines: array
                dist_tot,   \    # Total distance covered per particle
                time_tot,   \    # Total time covered per particle
                pth_data =  \    # Raw pathline data
                '''

                # Check, if axisymmetric or 2D --> y should be self.ymid[0]
                if (self.model_type == "axisymmetric") | (self.model_type == "2D"):
                
                    for iKey in xyz_nodes.keys():
                        xyz_nodes[iKey] = np.array([(round(idx_[0],3),
                                                    round(self.ymid[0],3),
                                                    round(idx_[2],3)) for idx_ in xyz_nodes[iKey]])
                    
                # col, lay, row index
                node_indices = self.get_node_indices(xyz_nodes = xyz_nodes)
                # Particle indices
                part_idx = [iPart for iPart in xyz_nodes]
                # Test array travel times (days)
                tot_time_arr = {iPart: time_diff[iPart].sum() for iPart in part_idx}

                for iPart in part_idx:
                    # Loop through converted rec.arrays of pth_data using part_idx 
                    # Fill recarray using pathline_data
                    particle_data[f"{iPG}-{iNode}-{iPart}"] = copy.deepcopy(pth_data[iPart][:]) # [:-1]
                    # Export rec.arrays as pd dataframe
                    df_particle_data[f"{iPG}-{iNode}-{iPart}"] = pd.DataFrame.from_records(data = particle_data[f"{iPG}-{iNode}-{iPart}"],
                                                                                index = "particleid",
                                                                                exclude = ["k"]).iloc[:,:]
                    # Drop duplicates                                                            
                    df_particle_data[f"{iPG}-{iNode}-{iPart}"] = df_particle_data[f"{iPG}-{iNode}-{iPart}"].drop_duplicates(
                        subset=["x","y","z","time"], keep = 'first') 

                    
                    for iParm in parm_list:

                        # Material property array
                        material_property_arr = getattr(self,iParm)
                        # dtype of array
                        mat_dtype = material_property_arr.dtype
                        # if mat_dtype == 'object':
                        #     dtype_ = '|S20'
                        # else:
                        dtype_ = mat_dtype
                            
                        # Numpy array values
                        # Use parm values for the first row in both 2D and axisymmetric models
                        if self.model_type in ["axisymmetric","2D"]:
                            parm_values = np.array([material_property_arr[idx[0],0,idx[2]] for idx in node_indices[iPart]])
                        else: # else: Use pathline data columns
                            parm_values = np.array([material_property_arr[idx[0],idx[1],idx[2]] for idx in node_indices[iPart]])


                        # Append recarray to particle_data (dict of dicts of np.recarray)
                        df_particle_data[f"{iPG}-{iNode}-{iPart}"].loc[:,iParm] = parm_values

                    # Change index name of df_particle                                                           
                    df_particle_data[f"{iPG}-{iNode}-{iPart}"].index.name = "flowline_id"
                    # Pseudonyms for df_particle column names
                    colnames_df_particle = {"x": "xcoord","y":"ycoord","z":"zcoord","time":"total_travel_time","prsity_uncorr":"porosity",
                                "solid_density": "solid_density", "fraction_organic_carbon": "fraction_organic_carbon", "redox": "redox", 
                                "dissolved_organic_carbon":	"dissolved_organic_carbon", "pH": "pH",	"temperature": "temperature", "material": "material"}
                    df_particle_data[f"{iPG}-{iNode}-{iPart}"].rename(columns = colnames_df_particle, 
                                                                                    inplace = True, errors = "raise")
                    df_particle_data[f"{iPG}-{iNode}-{iPart}"] = df_particle_data[f"{iPG}-{iNode}-{iPart}"].drop_duplicates(subset=["xcoord","ycoord","zcoord","total_travel_time"], keep = 'first')                                                                

                    # Append dataframes
                    df_particle_list.append(df_particle_data[f"{iPG}-{iNode}-{iPart}"])

        # Concatenate df_particle dataframes ('combined')
        if len(df_particle_list) > 0:
            df_particle = pd.concat(df_particle_list, axis = 0, ignore_index = False) #.drop_duplicates()
        else:
            colnames_df_particle = {"x": "xcoord","y":"ycoord","z":"zcoord","time":"total_travel_time","prsity_uncorr":"porosity",
                                "solid_density": "solid_density", "fraction_organic_carbon": "fraction_organic_carbon", "redox": "redox", 
                                "dissolved_organic_carbon":	"dissolved_organic_carbon", "pH": "pH",	"temperature": "temperature", "material": "material"}
            df_particle = pd.DataFrame(columns = colnames_df_particle)

        # y-coordinate equals 0.5 * self.delc[0] in axisymmetric or 2D model
        if self.model_type in ["axisymmetric","2D"]:
            df_particle.loc[:,"ycoord"] = 0.5 * self.delc[0]

        # remove duplicate records from df_particle
        df_particle = df_particle.drop_duplicates(
            subset=["xcoord","ycoord","zcoord","total_travel_time"], keep = 'first')

        return df_particle, df_particle_data
        

    # Fill df_flowline
    def fill_df_flowline(self, df_particle, model_cbc):
        '''
        Fill df_flowline 

        Parameters
        ----------
        Returns
        -------
        df_flowline: pandas.DataFrame
            Column 'flowline_id': Integer
            Column 'flowline_type': string
                Described the type of contamination associated with the flowline,
                either 'diffuse_source' or 'point_source'.
            Column 'flowline_discharge': Float
                Discharge associated with the flowline, [m3/d].
            Column 'particle_release_day': Float
            Column 'input_concentration'
            Column 'endpoint_id': Integer
                ID of Well (or drain) where the flowline ends.
            Column 'well_discharge': float
            Column 'substance': string
            Column 'removal_function': string
        '''

        # Default value for model_cbc
        if model_cbc is None:
            model_cbc = self.model_cbc

        # flowline ID as index
        flowline_id = list(df_particle.index.unique())

        # Column names in df_flowline
        colnames_df_flowline = ["flowline_type","flowline_discharge","particle_release_day","endpoint_id",
                        "well_discharge", "substance", "removal_function"]

        df_flowline = pd.DataFrame(index = flowline_id, columns = colnames_df_flowline)
        # Change df_flowline index name
        df_flowline.index.name = "flowline_id"

        ## Obtain flowline discharge ##
        # flux at starting location (grid cell) divided by total flowlines starting there

        # flow right face (frf) and flow lower face (flf) (third option: flow front face)
        frf, flf, fff = self.read_binarycbc_flow(model_cbc)
        # flux of pathlines
        flux_pathline = {}
        # startpoints and endpoints of flowlines
        node_start, node_end = {}, {}
        # startpoint and endpoint counter
        count_startpoints, count_endpoints = {}, {}
        # endpoint ids
        endpoint_id = {}

        # Removal function
        df_flowline.loc[:,"removal_function"] = self.schematisation.removal_function
        
        for fid in flowline_id:
            # X,Y,Z,time data per particle flowline 
            xyzt_data = df_particle.loc[df_particle.index == fid,["xcoord","ycoord","zcoord","total_travel_time"]].sort_values(by = "total_travel_time")
            
            # Startpoint and endpoint (XYZ-data) of particle flowlines
            startpoint = xyzt_data.loc[xyzt_data.total_travel_time == xyzt_data.total_travel_time.min(),
                                        ["xcoord","ycoord","zcoord"]].values.tolist()[0]
            endpoint = xyzt_data.loc[xyzt_data.total_travel_time == xyzt_data.total_travel_time.max(),
                                        ["xcoord","ycoord","zcoord"]].values.tolist()[0]
            # Nodes of start and endpoints
            node_start[fid] = self.xyz_to_layrowcol(xyz_point = startpoint, round_digits = 5)
            node_end[fid] = self.xyz_to_layrowcol(xyz_point = endpoint, round_digits = 5)
            # Count start points in cell
            count_startpoints[node_start[fid]] = count_startpoints.get(node_start[fid],0) + 1
            count_endpoints[node_end[fid]] = count_endpoints.get(node_end[fid],0) + 1

            # flowline_type
            df_flowline.loc[fid,"flowline_type"] = self.flowline_type[fid]




        for fid in flowline_id:
            if self.trackingdirection == "forward":

                if df_flowline.loc[fid,"flowline_type"] in ["diffuse_source",]:
                    # starting point is used to calculate flux of pathline (flux_pathline)
                    flux_pathline[fid] = round(self.calc_flux_cell(frf,flf,fff, loc = node_start[fid]) / \
                                        count_startpoints[node_start[fid]],4)
                    # flux_pathline[fid] = round(abs(frf[node_start[fid][0],node_start[fid][1],node_start[fid][2]]) + \
                    #                             abs(flf[node_start[fid][0],node_start[fid][1],node_start[fid][2]]) + \
                    #                             abs(fff[node_start[fid][0],node_start[fid][1],node_start[fid][2]]) / \
                    #                     count_startpoints[node_start[fid]],4)
                elif df_flowline.loc[fid,"flowline_type"] in ["point_source",]:
                    flux_pathline[fid] = self.point_discharge[fid]

                # flux_pathline[fid] = round(np.sqrt((frf[node_start[fid][0],node_start[fid][1],node_start[fid][2]]**2 + \
                #                             flf[node_start[fid][0],node_start[fid][1],node_start[fid][2]]**2 + \
                #                             fff[node_start[fid][0],node_start[fid][1],node_start[fid][2]]**2)) / \
                #                     count_startpoints[node_start[fid]],4)
            elif self.trackingdirection == "backward":

                if df_flowline.loc[fid,"flowline_type"] in ["diffuse_source"]:
                    # endpoint is used to calculate flux of pathline
                    flux_pathline[fid] = round(self.calc_flux_cell(frf,flf,fff, loc = node_end[fid]) / \
                                            count_endpoints[node_end[fid]],4)
                    # flux_pathline[fid] = round(abs(frf[node_end[fid][0],node_end[fid][1],node_end[fid][2]]) + \
                    #                             abs(flf[node_end[fid][0],node_end[fid][1],node_end[fid][2]]) + \
                    #                             abs(fff[node_end[fid][0],node_end[fid][1],node_end[fid][2]]) / \
                    #                     count_startpoints[node_end[fid]],4)
                elif df_flowline.loc[fid,"flowline_type"] in ["point_source",]:
                    flux_pathline[fid] = self.point_discharge[fid]
                # flux_pathline[fid] = round(np.sqrt((frf[node_end[fid][0],node_end[fid][1],node_end[fid][2]]**2 + \
                #                             flf[node_end[fid][0],node_end[fid][1],node_end[fid][2]]**2 + \
                #                             fff[node_end[fid][0],node_end[fid][1],node_end[fid][2]]**2)) / \
                #                     count_startpoints[node_end[fid]],4)
            # endpoint id
            endpoint_id[fid] = self.material[node_end[fid][0],node_end[fid][1],node_end[fid][2]]

        # fill flowline_df
        for fid in flowline_id:
            df_flowline.loc[fid,"flowline_discharge"] = flux_pathline[fid]
            df_flowline.loc[fid,"endpoint_id"] = endpoint_id[fid]

        for end_id in df_flowline["endpoint_id"].unique():
            # well (=endpoint) discharge (using cbc-file)       
            well_discharge = round(abs(frf[(self.material == end_id) & (self.ibound != 0)]).sum() + \
                                    abs(flf[(self.material == end_id) & (self.ibound != 0)]).sum() + \
                                    abs(fff[(self.material == end_id) & (self.ibound != 0)]).sum(), 4)
            
            # well_discharge = df_flowline.loc[df_flowline.endpoint_id == end_id,"flowline_discharge"].astype('float').values.sum()
            df_flowline.loc[df_flowline.endpoint_id == end_id,"well_discharge"] = well_discharge

        return df_flowline
     


    def plot_pathtimes(self,df_particle, 
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
            marker_size=0.05  # 's' in functie scatter
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

    def calc_advective_microbial_removal(self,df_particle,df_flowline, 
                                        endpoint_id, trackingdirection = "forward",
                                        mu1 = 0.023, grainsize = 0.00025, alpha0 = 1.E-5, pH0 = 6.8, const_BM = 1.38e-23,
                                        temp_water = 11., rho_water = 999.703, species_diam = 2.33e-8,
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
                                                    rho_water = rho_water, species_diam = species_diam)

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
        if type(perlen) != dict:
            self.perlen = {0: perlen}    
        else:
            self.perlen = perlen
        # Nr of time periods per stress period (int)
        if type(nstp) != dict:
            self.nstp = {0: nstp} 
        else:
            self.nstp = nstp 

        # Nr of stress periods  (int)
        self.nper = nper

        # Steady state model run (True/False)
        if type(steady) != dict:
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

        # Remove vadose_zone from geo_parameters and add to separate dict "vadose_parameters"
        self.schematisation_dict = self.extract_vadose_zone_parameters(schematisation = self.schematisation_dict,
                                                                        remove_keys = True)
        ## Following code only relevant if there are different input requirements for each schematisation type
        if self.schematisation_type in ["phreatic","semiconfined"]:
            self.model_type = "axisymmetric"
            print("Run phreatic model")
            self._check_schematisation_input()

        ## NOTE to programmer: for additional schematisation_types add to functionality here (elif)...

        elif self.schematisation_type in ["BAR"]:
            self.model_type = "2D"
            print("Run 2D model")
            self._check_schematisation_input()

        #### 29-11-'21: generalize the code based on required modflow_packages
        self.create_modflow_input()
        if self.run_mfmodel:

            # Run modflow model
            self.mfmodelrun()
            
        # Modpath simulation
        if self.run_mpmodel:

            # Create radial distance array with particle locations
            # self._create_radial_distance_array() # Analytische fluxverdeling
            self._create_diffuse_particles(recharge_parameters = 'concentration_boundary_parameters',
                                                            nparticles_cell = 1,
                                                            localy=0.5, localz=1.,
                                                            timeoffset=0.0, drape=0,
                                                            trackingdirection = self.trackingdirection,
                                                            releasedata=0.0)  

            self._create_point_particles(point_parameters = "point_parameters",
                                                   localx = 0.5, 
                                                   localy = 0.5, 
                                                   localz = 0.5,
                                                   timeoffset = 0.0, drape = 0,
                                                   trackingdirection = self.trackingdirection,  ## 'forward'
                                                   releasedata = 0.0)

            # Default flux interfaces
            defaultiface = {'RECHARGE': 6, 'ET': 6}
        
            ## Model mpbas input ##
            self.mpbas_input(prsity = self.porosity, defaultiface = defaultiface)

            # Run modpath model
            self.MP7modelrun()

            # self.success_mp = True
            print("modelrun of type", self.schematisation_type, "completed.")

            # Pathline output file
            self.mppth = os.path.join(self.workspace, self.modelname + '_mp.mppth')

            # Export output data to particle_df and flowline_df
            self._export_to_df(mppth = self.mppth)

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

