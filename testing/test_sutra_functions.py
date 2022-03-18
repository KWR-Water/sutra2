

import pytest
import numpy as np
import pandas as pd
import os
import sys
from pathlib import Path

from zmq import zmq_version_info

import greta.sutra_functions as SF


#%%
def test_mbo_removal_function(organism_name = "MS2"):

    # Calculate advective microbial removal
    modpath_removal = ST.SubstanceTransport(modpath_phrea,
                                            organism = organism_name)
 
    # Calculate advective microbial removal
    # Final concentration per endpoint_id
    C_final = {}
    for endpoint_id in modpath_phrea.schematisation_dict.get("endpoint_id"):
        df_particle, df_flowline, C_final[endpoint_id] = modpath_removal.calc_advective_microbial_removal(
                                            modpath_phrea.df_particle, modpath_phrea.df_flowline, 
                                            endpoint_id = endpoint_id,
                                            trackingdirection = modpath_phrea.trackingdirection,
                                            mu1 = 0.023, grainsize = 0.00025, alpha0 = 1.E-5, reference_pH = 6.8, const_BM = 1.38e-23,
                                            temp_water = 11., rho_water = 999.703, organism_diam = 2.33e-8,
                                            conc_start = 1., conc_gw = 0.)

