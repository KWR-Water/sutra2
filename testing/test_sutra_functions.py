

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
    mbo_removal = SF.MicrobialRemoval(organism = organism_name)
 
    # Calculate advective microbial removal
    C_final = mbo_removal.calc_advective_microbial_removal(grainsize = 0.00025,
                                            temp_water = 11., rho_water = 999.703, 
                                            conc_start = 1., conc_gw = 0.,
                                            redox = 'anoxic',
                                            distance_traveled = 100., time_diff = 10.,
                                            # organism_diam = 2.33e-8,
                                            # mu1 = 0.023,
                                            # alpha0 = 1.E-5,
                                            # reference_pH = 6.8
                                            )

    assert C_final > 0.

