
# MWK verwijder alle ongebruikte imports in alle files.
# gebruik tools als pylint en pylance hiervoor dan wordt het weergegeven in de GUI

from dis import dis
import pytest
import numpy as np
import pandas as pd
import os
import sys # MWK: unused import. install pylance (vscode extensie) en pylint (python pakket), dan worden dit sorrt dingen weergegeven
from pathlib import Path
from pandas.testing import assert_frame_equal # MWK unused import
import warnings #MWK unused import

from zmq import zmq_version_info # MWK unused import

import sutra2.removal_functions as rf

# get directory of this file
path = Path(__file__).parent


#%%
def test_scenarios_mbo_removal_function(organism_name = "MS2"):
    ''' Verify manual input for a species that is not yet available in the 'Organism'
        class (containing default removal parameters for some researched species).
    '''

    # Location of test file for scenarios
    # MWK:
    # - Deze file is niet in de repo. Wellicht moet je even .gitignore aanpassen zodat xlsx niet meer worden genegeerd
    #   dan deze file toevoegen aan de repo en dan xlsx files weer negeren in .gitignore.
    # - als het enigszins mogelijk is, hier een CSV van maken. CSV werkt samen met git. dit geldt voor alle xlsx in de
    #   hele repo (als het kan vervangen door CSV)
    # - naam beter test_berekening... ipv Testbereke...
    # - probeer alle input test data in een subfolder (e.g. input_data) te stoppen, dan blijft het overzichtelijker.
    scenarios_fpath = os.path.join(path,"Testberekeningen_sutra_mbo_removal_220321.xlsx")
    sheet_name = "Scenarios"
    # Read scenario excel file
    df_test = pd.read_excel(scenarios_fpath, sheet_name = sheet_name, skiprows = 1)

    # df_output
    columns_output = ["k_att","lambda","steady_state_concentration"]
    df_output = pd.DataFrame(index = df_test.index, columns = columns_output)

    for fid in df_test.index:
        # MWK: I learned something: at is faster than loc :D
        organism_name = organism_name
        redox = df_test.at[fid,'redox']
        alpha0 = df_test.at[fid,'alpha0']
        pH0 = df_test.at[fid,'pH0']
        mu1 = df_test.at[fid,'mu1']
        organism_diam = df_test.at[fid,'organism_diam']
        por_eff = df_test.at[fid,'porosity']
        grainsize = df_test.at[fid,'grainsize']
        pH_water = df_test.at[fid,'pH']
        temp_water = df_test.at[fid,'temperature']
        rho_water = df_test.at[fid,'rho_water']
        conc_start = 1.  # normally in df_flowline; use relative concentration as output
        conc_gw = 0.     # normally in df_flowline; use relative concentration as output
        distance_traveled = df_test.at[fid,'relative_distance']
        traveltime = df_test.at[fid,'total_travel_time']

        # Calculate advective microbial removal
        mbo_removal = rf.MicrobialRemoval(organism = organism_name)
        # Calculate final concentration after advective microbial removal
        C_final = mbo_removal.calc_advective_microbial_removal(grainsize = grainsize,
                                                temp_water = temp_water, rho_water = rho_water,
                                                pH = pH_water, por_eff = por_eff,
                                                conc_start = conc_start, conc_gw = conc_gw,
                                                redox = redox,
                                                distance_traveled = distance_traveled,
                                                traveltime = traveltime,
                                                organism_diam = organism_diam,
                                                mu1 = mu1,
                                                alpha0 = alpha0,
                                                pH0 = pH0)

        # k_att, calculated
        df_output.loc[fid,"k_att"] = mbo_removal.k_att
        # lambda, calculated
        df_output.loc[fid,"lambda"] = mbo_removal.lamda
        # (relative) concentration, calculated
        df_output.loc[fid,"steady_state_concentration"] = C_final

    # Calculate the difference between test dataframe and generated dataframe
    # NWK: dit is wel slim om te doen, je hebt ook een
    # https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.testing.assert_frame_equal.html
    # maar dat levert ook wel eens gezeur op. deze methode is lekker straigthforward
    diff_perc = np.abs((df_output.loc[:,columns_output].values - \
                        df_test.loc[:,columns_output].values) / \
                        df_test.loc[:,columns_output].values) * 100.

    assert not np.any(diff_perc > 0.5)
    #print("dataframe values differ to much: " + str(round(diff_perc.max(),2)) + " %")

#%%
def test_mbo_removal_function_check_default(organism_name = "carotovorum",
                                            redox = 'anoxic', # unused argument
                                            alpha0 = 0.577,
                                            pH0 = 7.5,
                                            mu1 = 0.1279,
                                            organism_diam = 1.803e-6,
                                            por_eff = 0.33,
                                            grainsize = 0.00025,
                                            pH_water = 7.5,
                                            temp_water = 11.,
                                            rho_water = 999.703,
                                            conc_start = 1., # unused argumetn
                                            conc_gw = 0.,# unused argumetn
                                            distance_traveled = 1.,
                                            traveltime = 100.):
    # MWK op zich een goede test, maar als je de default waardes aanpast moet je het hier ook aanpassen. Niet een
    # probleem maar wel een huishoudelijk dingetje wat in de gaten gehouden moet worden
    ''' Verify whether the default removal parameters is loaded successfully and gives
        the same result as manual input for the 'default parameters'.

    ## Default parameters: ##
    organism_name = "carotovorum"
    redox = 'anoxic',
    alpha0 = 0.577,
    pH0 = 7.5,
    mu1 = 0.1279,
    organism_diam = 1.803e-6,
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

    ## Default test
    # Calculate advective microbial removal
    mbo_removal_default = rf.MicrobialRemoval(organism = organism_name)
    # MWK Het is niet nodig om de hele advectie te bepalen als je alleen maar de default waardes wil checken
    # dat kun je gewoon doen door MicrobialRemoval instanties mbo_removal_default en mbo_removal_test en dan
    # de attributen te vergelijken waar die default waardes worden uitgevoerd
    # Calculate final concentration after advective microbial removal
    C_final_default= mbo_removal_default.calc_advective_microbial_removal()

    # Lambda (default): inactivation
    lambda_default = mbo_removal_default.lamda

    # Calculate advective microbial removal
    mbo_removal_test = rf.MicrobialRemoval(organism = organism_name)

    # Calculate final concentration after advective microbial removal
    C_final_test = mbo_removal_test.calc_advective_microbial_removal(grainsize = grainsize,
                                            temp_water = temp_water, rho_water = rho_water,
                                            pH = pH_water, por_eff = por_eff,
                                            conc_start = 1., conc_gw = 0.,
                                            redox = 'anoxic',
                                            distance_traveled = distance_traveled,
                                            traveltime = traveltime,
                                            organism_diam = organism_diam,
                                            mu1 = mu1,
                                            alpha0 = alpha0,
                                            pH0 = pH0
                                            )

    # Lambda: inactivation
    lambda_test = mbo_removal_test.lamda

    assert round(lambda_default,4) == round(lambda_test,4)
    assert round(C_final_default,4) == round(C_final_test,4)



# MWK, waarom heeft test functie een argument? Dat werkt toch aleen met fixtures?
def test_manual_input_mbo_removal(organism_name = "MS2"):
    # MWK wat is mbo ? microbiology?
    # volgens mij gaat het hier niet om manual input maar de test of het model
    # klopt bij een opgelegde input (waarvan je de output kent). dan kun je beter
    # de typering van deze input (anoxic MS2 oid, dit weet jij beter dan ik) meegeven in de naam van de functie
    # en zeggen dat het verificatie van model output is. dus functie naam wordt dan
    # iets als `test_verificatie_anoxic_MS2_case`
    # je moet dan de docstring overeenkomstig aanpassen.
    '''
    Verify manual input to function 'calc_advective_microbial_removal'
    against an earlier result.
    '''
    # test parameters
    por_eff = 0.33
    grainsize = 0.00025
    pH_water = 7.5
    temp_water = 10.
    rho_water = 999.703

    # Removal parms
    # alpha  'sticky coefficient'
    alpha0 = 0.001 # [-]
    pH0 = 7.5
    # MWK ik snap dit commentaar niet. of weghalen of uitleggen waarom het hier staat
    # --> if pH == pH0, then alpha == alpha0
    # coll_eff = 0.001

    # MWK als je zo;n omschrijving nodig hebt voor je variable, dan is de naam waarschijnlijk niet duidelijk genoeg
    # al is het hier wel lastig omdat het wel een hele lange naam wordt time_dep_inactivation_coeff moet wel kunnen
    # time dependent inactivation coefficient mu1 [day-1]
    mu1 = 0.149
    # org. diameter [m] # dit wordt eigenlijk al weergegeven in de variable naam
    organism_diam = 2.33e-8 # MWK beter alleen # [m] asl inline commentraar hier

    distance_traveled = 1.
    traveltime = 100.
    porewater_velocity = distance_traveled / traveltime # MWK unused variable

    # Calculate advective microbial removal
    mbo_removal = rf.MicrobialRemoval(organism = organism_name)
    # Calculate advective microbial removal
    C_final = mbo_removal.calc_advective_microbial_removal(grainsize = grainsize,
                                            temp_water = temp_water, rho_water = rho_water,
                                            pH = pH_water, por_eff = por_eff,
                                            conc_start = 1., conc_gw = 0.,
                                            redox = 'anoxic',
                                            distance_traveled = distance_traveled,
                                            traveltime = traveltime,
                                            organism_diam = organism_diam,
                                            mu1 = mu1,
                                            alpha0 = alpha0,
                                            pH0 = pH0
                                            )
    # Lambda: inactivation
    lamda = mbo_removal.lamda

    assert round(lamda,4) == round(0.7993188853572424 + mu1,4)

    assert round(C_final,3) == round(6.531818379725895e-42,3)


# MWK alle commentraar van hiervoven geldt hier ook
def test_manual_input_AW_mbo_removal(organism_name = "solani"):
    '''
    Verify manual input to function 'calc_advective_microbial_removal'
    against an earlier result.
    '''
    # test parameters
    organism_name = organism_name
    por_eff = 0.35
    grainsize = 0.00025
    pH_water = 7
    temp_water = 11.
    rho_water = 999.703

    # Removal parms
    # alpha  'sticky coefficient'
    alpha0 = 0.00037 # [-]
    pH0 = 7.5
    # --> if pH == pH0, then alpha == alpha0
    # coll_eff = 0.001

    # time dependent inactivation coefficient mu1 [day-1]
    mu1 = 0.1151
    # org. diameter [m]
    organism_diam = 2.73e-6

    distance_traveled = 5.
    traveltime = 2477.013
    porewater_velocity = distance_traveled / traveltime

    # Calculate advective microbial removal
    mbo_removal = rf.MicrobialRemoval(organism = organism_name)
    # Calculate advective microbial removal
    C_final = mbo_removal.calc_advective_microbial_removal(grainsize = grainsize,
                                            temp_water = temp_water, rho_water = rho_water,
                                            pH = pH_water, por_eff = por_eff,
                                            conc_start = 1., conc_gw = 0.,
                                            redox = 'anoxic',
                                            distance_traveled = distance_traveled,
                                            traveltime = traveltime,
                                            organism_diam = organism_diam,
                                            mu1 = mu1,
                                            alpha0 = alpha0,
                                            pH0 = pH0
                                            )
    # Lambda: inactivation
    lamda = mbo_removal.lamda

    assert round(lamda,4) == round(0.019023489784646908 + mu1,4)

    assert round(C_final,3) == round(5.2028714869195504e-145,3)