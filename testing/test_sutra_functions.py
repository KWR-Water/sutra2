

from dis import dis
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
 
    # Calculate final concentration after advective microbial removal
    C_final,_,_ = mbo_removal.calc_advective_microbial_removal(grainsize = 0.00025,
                                            temp_water = 11., rho_water = 999.703, 
                                            conc_start = 1., conc_gw = 0.,
                                            redox = 'anoxic',
                                            distance_traveled = 100., traveltime = 10.,
                                            # organism_diam = 2.33e-8,
                                            # mu1 = 0.023,
                                            # alpha0 = 1.E-5,
                                            # reference_pH = 6.8
                                            )

    assert C_final > 0.

#%%
def test_mbo_removal_function_check_default(organism_name = "MS2",
                                            redox = 'anoxic',
                                            alpha0 = 1.e-5,
                                            reference_pH = 6.8,
                                            mu1 = 0.023,
                                            organism_diam = 2.33e-8,
                                            por_eff = 0.33,
                                            grainsize = 0.00025,
                                            pH_water = 7.5,
                                            temp_water = 11.,
                                            rho_water = 999.703,
                                            conc_start = 1.,
                                            conc_gw = 0.,
                                            distance_traveled = 1.,
                                            traveltime = 100.):

    ## Default test
    # Calculate advective microbial removal
    mbo_removal_default = SF.MicrobialRemoval(organism = organism_name)
    # Calculate advective microbial removal
    C_final_default, lambda_default, k_att_default = mbo_removal_default.calc_advective_microbial_removal()

    '''
    ## Default parameters: ##
    redox = 'anoxic',
    alpha0 = 1.e-5,
    reference_pH = 6.8,
    mu1 = 0.023,
    organism_diam = 2.33e-8,
    por_eff = 0.33,
    grainsize = 0.00025,
    pH_water = 7.5,
    temp_water = 11.,
    rho_water = 999.703,
    conc_start = 1.,
    conc_gw = 0.,
    distance_traveled = 1.,
    traveltime = 100.
    '''

    # Calculate advective microbial removal
    mbo_removal = SF.MicrobialRemoval(organism = organism_name)
    # Calculate advective microbial removal
    C_final_test, lambda_test, k_att_test = mbo_removal.calc_advective_microbial_removal(grainsize = grainsize,
                                            temp_water = temp_water, rho_water = rho_water,
                                            pH = pH_water, por_eff = por_eff, 
                                            conc_start = 1., conc_gw = 0.,
                                            redox = 'anoxic',
                                            distance_traveled = distance_traveled, 
                                            traveltime = traveltime,
                                            organism_diam = organism_diam,
                                            mu1 = mu1,
                                            alpha0 = alpha0,
                                            reference_pH = reference_pH
                                            )

    assert round(k_att_default,4) == round(k_att_test,4) 
    assert round(lambda_default,4) == round(lambda_test,4) 

    assert round(C_final_default,4) == round(C_final_test,4)



def test_manual_input_mbo_removal(organism_name = "MS2"):

    # test parameters
    organism_name = organism_name
    por_eff = 0.33
    grainsize = 0.00025
    pH_water = 7.5
    temp_water = 10.
    rho_water = 999.703

    # Removal parms
    # alpha  'sticky coefficient'
    alpha0 = 0.001 # [-]
    reference_pH = 7.5
    # --> if pH == reference_pH, then coll_eff == alpha0
    # coll_eff = 0.001

    # time dependent inactivation coefficient mu1 [day-1]
    mu1 = 0.149
    # org. diameter [m]
    organism_diam = 2.33e-8

    distance_traveled = 1.
    traveltime = 100.
    porewater_velocity = distance_traveled / traveltime

    # Calculate advective microbial removal
    mbo_removal = SF.MicrobialRemoval(organism = organism_name)
    # Calculate advective microbial removal
    C_final, lambda_, k_att = mbo_removal.calc_advective_microbial_removal(grainsize = grainsize,
                                            temp_water = temp_water, rho_water = rho_water,
                                            pH = pH_water, por_eff = por_eff, 
                                            conc_start = 1., conc_gw = 0.,
                                            redox = 'anoxic',
                                            distance_traveled = distance_traveled, 
                                            traveltime = traveltime,
                                            organism_diam = organism_diam,
                                            mu1 = mu1,
                                            alpha0 = alpha0,
                                            reference_pH = reference_pH
                                            )

    assert round(k_att,4) == round(0.7993188853572424,4) 
    assert round(lambda_,4) == round(0.7993188853572424 + mu1,4) 

    assert round(C_final,3) == round(6.531818379725895e-42,3)
    