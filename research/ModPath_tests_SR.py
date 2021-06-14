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

def make_radial_discretisation(schematisation: dict, dict_keys = None,
                                res_hor_max = 10.):
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
    delr = np.round(np.diff(col_bounds),3)
    # Number of model columns
    ncol = len(delr)
    # Create xmid array from delr array
    xmid = np.empty((ncol), dtype= 'float')
    xmid[0] = delr[0] * 0.5
    for iCol in range(1, ncol):
        xmid[iCol] = (xmid[(iCol - 1)] + ((delr[(iCol)]) + (delr[(iCol - 1)])) * 0.5)  

    return delr, xmid, refinement_bounds


#%%

model_top, model_bot,model_thickness, model_radius = calc_modeldim_phrea(phreatic_scheme["geo_parameters"])
''' change mesh refinement dicts SR '''
phreatic_scheme['geo_parameters']['mesh_refinement1'] =   {"rmin": 0.,
                                         "rmax": 0.5,
                                         "res_hor": 0.025}
phreatic_scheme['geo_parameters']['mesh_refinement2'] =   {"rmin": 0.5,
                                         "rmax": model_thickness,
                                         "res_hor": 0.25}
phreatic_scheme['geo_parameters']['mesh_refinement3'] =   {"rmin": model_thickness,
                                         "rmax": model_radius,
                                         "res_hor": 5.}

# Refinement boundaries and column boundaries
# Horizontal discretisation dictionary keys
dict_keys = ["geo_parameters","recharge_parameters","ibound_parameters",
                      "well_parameters"]
delr, xmid, refinement_bounds = make_radial_discretisation(schematisation = phreatic_scheme,
                                                            dict_keys = dict_keys)
print(delr, xmid, refinement_bounds)
print(len(delr), len(xmid), len(refinement_bounds))

#%%
modpath_phrea = ModPathWell(phreatic_scheme)
modpath_phrea.schematisation
# Create bas_parameters
modpath_phrea.assign_bas_parameters()
modpath_phrea.bas_parameters

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


# check of creating delr in make_radial_discretisation
ncols_filterscreen = None
ncols_gravelpack = None
ncols_near_well = None
ncols_far_well = None
diameter_filterscreen = 0.1
diameter_gravelpack = 0.75
thickness_aquifer = 20.
model_radius = 500.

# Calculated parameters (resulting in actual model column widths)
delr_filterscreen = (0.5 * diameter_filterscreen)
delr_gravelpack = 0.5 * (diameter_gravelpack - diameter_filterscreen) 
delr_near_well = (thickness_aquifer - 0.5 * diameter_gravelpack)
delr_far_well = (model_radius - thickness_aquifer)

if ncols_gravelpack is None:
    ncols_filterscreen = 1
if ncols_gravelpack is None:
    ncols_gravelpack = max(1,math.ceil(delr_gravelpack/0.025))
if ncols_near_well is None:
    ncols_near_well = max(1,math.ceil(delr_near_well/0.25))
if ncols_far_well is None:
    ncols_far_well = max(1,math.ceil(delr_far_well/2.5))

# Use an iterator 'list' of delr-values to create delr attribute (array)
delr_list = [delr_filterscreen,delr_gravelpack, 
            delr_near_well, delr_far_well]
# List of number of columns
ncol_list = [ncols_filterscreen, ncols_gravelpack,
             ncols_near_well,ncols_far_well] 
# Total number of model culumns
ncols = sum(ncol_list)
# Create "delr" attribute from column widths (m)                                                                    num = ncol_list[idx] + 1)[1:]
delr = np.zeros((ncols), dtype = 'float')
for idx in range(len(delr_list)):
    if idx == 0: 
        delr[0:ncol_list[idx]] = np.diff(np.linspace(0.,delr_list[idx], 
                                                       num = ncol_list[idx] + 1))
    
    else: # linear spacing
        delr[sum(ncol_list[:idx]):sum(ncol_list[:idx+1])] = np.diff(np.linspace(sum(delr_list[:idx]),
                                                                    sum(delr_list[:idx+1]),
                                                                    num = ncol_list[idx] + 1))
    # else: # log sspacing
    #     delr[sum(ncol_list[:idx]):sum(ncol_list[:idx+1])] = np.diff(np.logspace(np.log10(sum(delr_list[:idx])),
    #                                                                 np.log10(sum(delr_list[:idx+1])),
    #                                                                 num = ncol_list[idx] + 1))


        
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

