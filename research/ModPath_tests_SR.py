# -*- coding: utf-8 -*-
"""
Created on Mon May 31 10:20:39 2021

@author: rosst
"""
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

from scipy.special.orthogonal import chebyc
try:
    from tqdm import tqdm  # tqdm gives a progress bar for the simultation
except Exception as e:
    print(e)
    pass
# import pyarrow.parquet as pq
import math
import re # regular expressions
from scipy.special import kn as besselk

path = os.getcwd()  # path of working directory

# Import schematisation functions

from pathlib import Path
try:
    from project_path import module_path #the dot says looik in the current folder, this project_path.py file must be in the folder here
except ModuleNotFoundError:
    from project_path import module_path

from greta.draft_transport_function import *
from greta.ModPath_functions import ModPathWell

from testing.test_transatomic import *
# get directory of this file
path = Path(__file__).parent #os.getcwd() #path of working directory

#%%
# output Alex "phreatic_dict_nogravel.txt" --> saved in "testing dir"
research_dir = os.path.join(path,"..","research")
with open(os.path.join(research_dir,"phreatic_dict_nogravel.txt"),"r") as file_:
    dict_content = file_.read()
    phreatic_scheme = eval(dict_content)

check_schematisation = False # Check input dict (T/F)
if check_schematisation:
    for iKey,iVal in phreatic_scheme.items():
        print(iKey,iVal,"\n")


#%%

# # check of creating delr in make_radial_discretisation
# ncols_filterscreen = None
# ncols_gravelpack = None
# ncols_near_well = None
# ncols_far_well = None
# diameter_filterscreen = 0.1
# diameter_gravelpack = 0.75
# thickness_aquifer = 20.
# model_radius = 500.


def calc_modeldim_phrea(geo_parameters: dict):
    ''' Calculate the model dimensions using geo_parameters:
        - model_top (model top in m)
        - model_bot (model bottom in m)
        - model_thickness (vertical extent (thickness) of the model in m)
        - model_radius (horizontal extent (distance from well to outer boundary)
    '''
    model_top = None
    model_bot = None
    model_thickness = None
    model_radius = None
    for iDict in geo_parameters.keys():
        if "top" in geo_parameters[iDict]:
            if model_top is None:
                model_top = geo_parameters[iDict]["top"]
            elif model_top < geo_parameters[iDict]["top"]:
                # Determine maximum top value from dicts
                model_top = geo_parameters[iDict]["top"]
        if "bot" in geo_parameters[iDict]:
            if model_bot is None:
                model_bot = geo_parameters[iDict]["bot"]
            elif model_bot > geo_parameters[iDict]["bot"]:
                # Determine minimum bottom value from dicts
                model_bot = geo_parameters[iDict]["bot"]
        if "rmax" in geo_parameters[iDict]:
            if model_radius is None:
                model_radius = geo_parameters[iDict]["rmax"]
            elif model_radius < geo_parameters[iDict]["rmax"]:
                # Determine maximum model radius from dicts
                model_radius = geo_parameters[iDict]["rmax"]

        # Calculate model thickness from top and bottom values
        model_thickness = model_top - model_bot
    return model_top, model_bot, model_thickness, model_radius

def _assign_cellboundaries(schematisation: dict, dict_keys: list or None = None,
                            bound_min: str = "xmin", bound_max: str = "xmax",
                            n_refinement: str = "ncols", ascending: bool = True,
                            res_max: int or float or None = None):
    
    # Keep track of grid boundaries
    bound_list = []

    if dict_keys is None:
        dict_keys = [iDict for iDict in schematisation.keys()]

    # Loop through schematisation keys (dict_keys)
    for iDict in dict_keys:
        # Loop through subkeys of schematisation dictionary
        for iDict_sub in schematisation[iDict]:   

            try:
                # minimum bound
                val_min = schematisation[iDict][iDict_sub][bound_min]
            except KeyError as e:
                print(e,"continue")
                continue

            try:
                # maximum bound
                val_max = schematisation[iDict][iDict_sub][bound_max]
            except KeyError as e:
                print(e,"continue")
                continue

            try:
                # number of refinements
                n_ref = schematisation[iDict][iDict_sub][n_refinement]
            except KeyError:
                if res_max is None:
                    n_ref = 1
                else:
                    # Limit the cell resolution using 'res_max' (if not None)
                    n_ref = max(1,math.ceil((bound_max-bound_min)/res_max))
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
        center_points[0] = cell_sizes[0] * 0.5
        for idx in range(1, len_arr):
            center_points[idx] = (center_points[(idx - 1)] + ((cell_sizes[(idx)]) + (cell_sizes[(idx - 1)])) * 0.5)  
    else: # zmid arrays are decreasing with increasing index number
        center_points[0] = cell_sizes[0] * 0.5
        for idx in range(1, len_arr):
            center_points[idx] = (center_points[(idx - 1)] - ((cell_sizes[(idx)]) + (cell_sizes[(idx - 1)])) * 0.5)  

    return len_arr, cell_sizes, center_points, bound_list


def make_discretisation(schematisation: dict, dict_keys = None,
                        axisym = True):
    ''' Generate distance between columns for axisymmetric model.
    Sets it to self.delr
    - schematisation is of type dict with (sub)dictionaries with keys 'dict_keys'.
    - res_hor_max: maximum horizontal resolution if horizontal resolution 'res_hor'
      is not given in the subdictionary.

    The function returns numpy arrays of: 
    - delr: column widths of the model columns [np.array]
      rounded to three decimals [mm scale].
    - xmid: x-coordinates (middle) of the model columns [np.array]
    - refinement_bounds: boundaries where horizontal resolution values may change.
    '''
    # mesh refinement boundaries
    refinement_bounds = {"X": [],
                         "Y": [],
                         "Z": []}
    # Grid boundaries along X,Y,Z axes
    grid_bounds = {"X": [],
                   "Y": [],
                   "Z": []}
    # Assign delr and xmid
    ncol, delr, xmid,_ = _assign_cellboundaries(schematisation, dict_keys,
                               bound_min = "xmin", bound_max = "xmax",
                               n_refinement = "ncols", ascending = True)
    # Assign delv and zmid    
    nlay, delv, zmid, lay_bounds = _assign_cellboundaries(schematisation, dict_keys,
                               bound_min = "bot", bound_max = "top",
                               n_refinement = "nlayers", ascending = False)                                   
    # Model top
    top = max(lay_bounds)
    # Model bottoms
    bot = top - delv.cumsum() 

    # Assign delc and ymid
    if axisym:
        nrow = 2
        delc = np.ones((nrow),dtype = 'float')
        ymid = np.array([0.5,1.5])
    else:
        nrow,delc,ymid,_ = _assign_cellboundaries(schematisation, dict_keys,
                               bound_min = "ymin", bound_max = "ymax",
                               n_refinement = "nrows", ascending = True)
    
    
    return nlay,nrow,ncol, delv,delc,delr, zmid,ymid,xmid,top,bot



#%%


# # Method 1: horizontal discretisation based on number of columns
# # phreatic_scheme['geo_parameters']['mesh_refinement1'] =   {"rmin": 0.,
# #                                          "rmax": 0.5,
# #                                          "ncols": 20}
# # phreatic_scheme['geo_parameters']['mesh_refinement2'] =   {"rmin": 0.5,
# #                                          "rmax": model_thickness,
# #                                          "ncols": 20}
# # phreatic_scheme['geo_parameters']['mesh_refinement3'] =   {"rmin": model_thickness,
# #                                          "rmax": model_radius,
# #                                          "ncols": 40}

# # Refinement boundaries and column boundaries
# # Horizontal discretisation dictionary keys
# dict_keys = ["geo_parameters","recharge_parameters","ibound_parameters",
#                       "well_parameters"]

                      
# delr, xmid, refinement_bounds = make_radial_discretisation(schematisation = phreatic_scheme,
#                                                             dict_keys = dict_keys,
#                                                             res_hor_max = 10.)
# print(delr, xmid, refinement_bounds)
# print(len(delr), len(xmid), len(refinement_bounds))

# #%%
# # Method 2: horizontal and vertical discretisation based on grid resolution
# ''' change mesh refinement dicts SR '''
# phreatic_scheme['geo_parameters']['mesh_refinement1'] =   {"rmin": 0.,
#                                          "rmax": 0.5,
#                                          "res_hor": 0.025}
# phreatic_scheme['geo_parameters']['mesh_refinement2'] =   {"rmin": 0.5,
#                                          "rmax": model_thickness,
#                                          "res_hor": 0.25}
# phreatic_scheme['geo_parameters']['mesh_refinement3'] =   {"rmin": model_thickness,
#                                          "rmax": model_radius,
#                                          "res_hor": 5.}

# # Refinement boundaries and column boundaries
# # Horizontal discretisation dictionary keys
# # dict_keys = ["geo_parameters","recharge_parameters","ibound_parameters",
# #                       "well_parameters"]
# # delr, xmid, refinement_bounds = make_radial_discretisation(schematisation = phreatic_scheme,
# #                                                             dict_keys = dict_keys,
# #                                                             res_hor_max = None)
# # print(delr, xmid, refinement_bounds)
# # print(len(delr), len(xmid), len(refinement_bounds))

#%%
''' Add changes to dictionary scheme (for programming purposes). '''

model_top, model_bot,model_thickness, model_radius = calc_modeldim_phrea(phreatic_scheme["geo_parameters"])
''' change mesh refinement dicts SR '''
phreatic_scheme['geo_parameters']['mesh_refinement1'] =   {"rmin": 0.,
                                         "rmax": 0.5,
                                         "res_hor": 0.025,
                                         "ncols": 20}
phreatic_scheme['geo_parameters']['mesh_refinement2'] =   {"rmin": 0.5,
                                         "rmax": model_thickness,
                                         "res_hor": 0.25,
                                         "ncols": 100}
phreatic_scheme['geo_parameters']['mesh_refinement3'] =   {"rmin": model_thickness,
                                         "rmax": model_radius,
                                         "res_hor": 5.,
                                         "ncols": 100}
                                         
dict_keys = ["geo_parameters","recharge_parameters","ibound_parameters",
                      "well_parameters"]

# Change well discharge of well1 to -7665.6
phreatic_scheme["well_parameters"]['well1']['Q'] = -7665.6
# phreatic_scheme["well_parameters"]['well1']['res_vert'] = 0.5
phreatic_scheme["well_parameters"]['well1']['xmin'] = phreatic_scheme["well_parameters"]['well1']['rmin']
phreatic_scheme["well_parameters"]['well1']['xmax'] = phreatic_scheme["well_parameters"]['well1']['rmax']
# Add test well_parameters with leak from 9.9 to 10.0 m
phreatic_scheme["well_parameters"]['well_leak'] = {'Q': -1.,
  'top': 10.0,
  'bot': 9.9,
  'rmin': 0.0,
  'rmax': 0.1,
  'nlayers': 1}

# Add ibound parameters
phreatic_scheme["ibound_parameters"]["ibound_type"] = -1
# # Create model discretisation using schematisation dict
# nlay,nrow,ncol,delv,delc,delr,zmid,ymid,xmid,top,bot = make_discretisation(schematisation = phreatic_scheme,
#                                                             dict_keys = dict_keys)
# print(nlay,nrow,ncol,delv,delc,delr,zmid,ymid,xmid,top,bot)

# # Add delv, zmid, top and bottoms
# delv, zmid, top, bot, refinement_bounds = make_vertical_discretisation(schematisation = phreatic_scheme,
#                                                     dict_keys = dict_keys,
#                                                     res_vert_max = None)

# Gaat het goed met freatische winning? --> grid + toewijzen concentraties


#%%
''' Inititalize ModPath class'''
modpath_phrea = ModPathWell(phreatic_scheme)
# modpath_phrea.schematisation

# Refinement boundaries and column boundaries
# Horizontal discretisation dictionary keys
dict_keys = ["geo_parameters","recharge_parameters","ibound_parameters",
                      "well_parameters"]
# Adds to object: delr, ncol, xmid
# modpath_phrea.make_discretisation(dict_keys = dict_keys)
# Adds to object: delv, nlay, top, bot, zmid
# modpath_phrea.make_vertical_discretisation(dict_keys = dict_keys)

# Print all attributes in object
# print(modpath_phrea.__dict__)
# modpath_phrea.phreatic()
modpath_phrea.run_model()
#### HIER GEBLEVEN 14-6-2021 ####

# Check attributes in ModPath object
check_attr_list = ["nlay","nrow","ncol","delv","delc","delr","zmid","ymid","xmid","top","bot"]
for iAttr in check_attr_list:
    print(iAttr,getattr(modpath_phrea,iAttr))

modpath_phrea.schematisation["ibound_parameters"]
# # Create bas_parameters
# modpath_phrea.assign_bas_parameters()
# modpath_phrea.bas_parameters

#%%
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
# Horizontal boundaries determined by rmin, rmax
# 
# Refinement1 --> near well 


#%%
# Horizontal spacing based on rmin, rmax, ncol

# output Alex "phreatic_dict_nogravel.txt" --> saved in "testing dir"
path = 'd:\\Greta\\greta\\greta'
testing_dir = os.path.join(path,"..","testing")
with open(os.path.join(testing_dir,"phreatic_dict_nogravel.txt"),"r") as file_:
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
# Unique boundaries
delr_bounds = sorted(list(set(delr_bounds)))
# List of distances between horizontal boundaries
delr_list = np.diff(delr_bounds)
# Number of columns based on pre-specified distances
ncol_list = []
for i_delr in delr_list:
    # cumulative radius
    radius_min = sum(delr_list[:i_delr])
    radius_max = sum(delr_list[:i_delr+1])
    if radius_min > 1.

for iBound in delr_bounds:
    
    ncol, rmin, rmax = iBound[0:3]
    dr = (rmax-rmin) / ncol
    print(rmin, rmax, dr)

    

# Proposal
# Horizontal boudnaries determined by rmin, rmax
# 
# Refinement1 --> near well 

# Refinement2 --> well --> aquifer thickness

#%%
# Vertical spacing based on top, bot, ncol

# output Alex "phreatic_dict_nogravel.txt" --> saved in "testing dir"
path = 'd:\\Greta\\greta\\greta'
testing_dir = os.path.join(path,"..","testing")
with open(os.path.join(testing_dir,"phreatic_dict_nogravel.txt"),"r") as file_:
    dict_content = file_.read()
    phreatic_scheme = eval(dict_content)

for iKey,iVal in phreatic_scheme.items():
    print(iKey,iVal,"\n")

delv_bounds = []
for iDict_main in phreatic_scheme:
    
    if iDict_main in ["geo_parameters","ibound_parameters",
                      "well_parameters"]:
        
        for iDict_sub in phreatic_scheme[iDict_main]:
            if (iDict_main == "geo_parameters") & (iDict_sub == "vadose"):
                continue
            print(iDict_main,iDict_sub)
            try:
                delv_bounds.append(phreatic_scheme[iDict_main][iDict_sub]["top"])
            except Exception as e: print(e)
            try:
                delv_bounds.append(phreatic_scheme[iDict_main][iDict_sub]["bot"])
            except Exception as e: print(e)

# Unique boundaries
delv_bounds = sorted(list(set(delv_bounds)))
print(delv_bounds)

nlay_list = [max(1,math.ceil((delv_bounds[idx]-delv_bounds[idx-1]) / 0.5)) for idx in range(1,len(delv_bounds))]
ncol_tot = sum(nlay_list)
delv = np.zeros((ncol_tot), dtype = 'float') + 0.5

