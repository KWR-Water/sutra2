========================================================================================================
Tutorial
========================================================================================================

Overview
========

SSTR (SubSurface TRansport) is a model to calculate the behavior of Organic
MicroPollutants (OMPs) and pathogens for 4 standard types of Public Supply Well
Fields (PSWFs), which cover the most frequently occurring and most vulnerable
groundwater resources for drinking water supply in the Netherlands (and Flanders).
One of the aims of this approach is to forecast the behavior of new OMPs in
groundwater. Groundwater is often overlooked in current environmental risk
assessment methods, which are a priori or a posteriori applied when new organic
chemicals appear on the market

The 4 standard PSWF types consist of a phreatic, a semiconfined, a Basin Artificial
Recharge (BAR) and River Bank Filtration (RBF) well field, each predefined with
representative, standard hydrogeological, hydrological and hydrogeochemical
characteristics.

This python version is based on the Lite+ version of the OMP transport model 'TRANSATOMIC'
(acronym: TRANS Aquifer Transport Of MicroContaminants, developed by P.Stuyfzand)
in which concentration changes are calculated with analytical solutions set in Excel spreadsheet.

The model has been expanded to include Modflow solutions, in addition to the analytical
solutions and to include microbial organisms ('mbo') in addition to OMP.

Steps
------
Operating the analytical module typically involves 5 steps:

#. Define the hydrogeochemical system using the HydroChemicalSchematisation class. 
#. Run the AnalyticalWell class to calculate the travel time distribution in the different aquifer zones
#. Run the Substance class to retrieve the substance (removal) parameters
#. Run the SubstanceTransport class to calculate the removal and concentration in each zone and in the well
#. Plot the desired functions

Now, let’s try some examples. First we import the necessary python packages

.. ipython:: python

    import pandas as pd
    from pathlib import Path
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import os
    from pandas import read_csv
    from pandas import read_excel
    import math
    from scipy.special import kn as besselk
    from pathlib import Path
    import sutra2
    import sutra2.Analytical_Well as AW
    import sutra2.ModPath_Well as mpw
    import sutra2.Transport_Removal as TR
    # path = Path(__file__).parent # path of working directory


Step 1: Define the HydroChemicalSchematisation
==============================================
The first step is to define the hydrogeochemistry of the system using the HydroChemicalSchematisation class.
In this class you specify the:

    * Computational method ('analytical' or 'modpath').
    * The schematisation type ('phreatic', 'semiconfined') 
    * schematisation types 'riverbankfiltration', 'basinfiltration' yet to be supported
    * The removal function ('omp' or 'mbo').
    * Input the relevant parameters for the porous media, the hydrochemistry, hydrology and the contamination of interest

.. ('riverbankfiltration', 'basinfiltration' coming soon).

.. The class parameters can be roughly grouped into the following categories;

.. * System.
.. * Settings.
.. * Porous Medium
.. * Hydrochemistry
.. * Hydrology
.. * Contaminant
.. * Diffuse contamination
.. * Point Contamination
.. * Model size

Units of input are:
.. * Discharge : m3/d
.. * Time: days
.. * Length: meters
.. * Concentration: ug/L
.. * Temperature: degree C
.. * Depth: meters above sea level (m ASL)
.. * Density: kg/L
.. * DOC/TOC: mg/L

Lets start with a simple example defining a HydroChemicalSchematisation object for a phreatic aquifer:

.. ipython:: python

    phreatic_schematisation = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                                        computation_method = 'modpath',
                                                        well_discharge=-7500, #m3/day
                                                        recharge_rate=0.0008, #m/day
                                                        thickness_vadose_zone_at_boundary=5, #m
                                                        thickness_shallow_aquifer=10,  #m
                                                        thickness_target_aquifer=40, #m
                                                        hor_permeability_target_aquifer=35, #m/day
                                                        redox_vadose_zone='anoxic',
                                                        redox_shallow_aquifer='anoxic',
                                                        redox_target_aquifer='deeply_anoxic',
                                                        pH_target_aquifer=7.,
                                                        temp_water=11.,
                                                        name='benzene',
                                                        diffuse_input_concentration = 100, #ug/L
                                                        )

The parameters from the HydroChemicalSchematisation class are added as attributes to
the class and can be accessed for example:

.. ipython:: python

    phreatic_schematisation.schematisation_type
    phreatic_schematisation.well_discharge
    phreatic_schematisation.porosity_shallow_aquifer

If not defined, default values are used for the rest of the parameters. To view all parameters in the schematisation:

.. ipython:: python

    phreatic_schematisation.__dict__

Then, we create a ModpathWell object for the HydroChemicalSchematisation object that we just made.
The ModpathWell object requires a dictionary of the subsurface schematisation and a set of boundary conditions
the numerical model has to abide by in calculating flow velocity and direction of flow.

.. ipython:: python

    phreatic_schematisation.make_dictionary()

To view the created dictionary use the following snippet of code.

.. ipython:: python

    schematisation_dict = {'simulation_parameters' : phreatic_schematisation.simulation_parameters,
        'endpoint_id': phreatic_schematisation.endpoint_id,
        'mesh_refinement': phreatic_schematisation.mesh_refinement,
        'geo_parameters' : phreatic_schematisation.geo_parameters,
        'ibound_parameters' : phreatic_schematisation.ibound_parameters,
        'recharge_parameters' : phreatic_schematisation.recharge_parameters,
        'well_parameters' : phreatic_schematisation.well_parameters,
        'point_parameters' : phreatic_schematisation.point_parameters,
        'concentration_boundary_parameters' : phreatic_schematisation.concentration_boundary_parameters,
    }
    schematisation_dict

The schematisation dict contains the following data:
.. * simulation_parameters: simulation data such as schematisation_type and computation_method
.. * endpoint_id: object location to compute final concentration for after removal like 'well1'
.. * mesh_refinement: optional additional grid refinement parameters
.. * geo_parameters: chemical/material data for creating geological layers [porosity,hydraulic conductivity,foc,DOC, pH, etc,]
.. * ibound_parameters: boundary conditions for flow
.. * recharge_parameters: groundwater recharge [unit: m] in a specified region
.. * well_parameters: collection of well locations and discharge to simulate.
.. * point_parameters: (starting) point source contamination(s) to calculate removal for
.. * concentration_boundary_parameters: diffuse contamination(s) to calculate removal for

.. ipython:: python

    print(os.getcwd())
    import pathlib
 
    path_cur = pathlib.Path('.')
    
    full_path = path_cur.absolute()
    
    my_path = full_path.as_posix()
    
    print(my_path)
    
    type(my_path)

Step 2: Run the ModpathWell class
=====================================
Next we create an ModpathWell object for the HydroChemicalSchematisation object we just made.
The data files will be stored in location workspace using a given modelname.

.. ipython:: python

    package_folder = Path(sutra2.__file__).parent
    
    mf_exe = os.path.join(package_folder,"mf2005.exe")
    mp_exe = os.path.join(package_folder,"mpath7.exe")

    p = Path(package_folder).glob('**/*')
    files = [x for x in p if x.is_file()]

    package_folder
    print(package_folder)
    mf_exe
    print(mf_exe)
    mp_exe
    print(mp_exe)
    print(files)

    mf_exe_pos = package_folder / "mf2005.exe" 
    print(os.path.exists(mf_exe_pos))
    print(mf_exe_pos)
    mf_exe_pos

    env_var = os.environ
    import pprint
    # Print the list of user's
    # environment variables
    print("User's Environment variable:")

    # pwd = os.environ.get("PWD")
    # p = Path(pwd).glob('**/*')
    # files = [x for x in p if x.is_file()]
    # print(files)

    rtd_venv = r"/home/docs/checkouts/readthedocs.org/user_builds/sutra2/envs/latest"
    mf_exe_rtd = os.path.join(rtd_venv, "sutra2","mf2005.exe")
    mp_exe_rtd = os.path.join(rtd_venv, "sutra2","mpath7.exe")

    mf_exe_git = r"https://github.com/KWR-Water/sutra2/blob/main/sutra2/mf2005.exe"
    mp_exe_git = r'https://github.com/KWR-Water/sutra2/blob/main/sutra2/mpath7.exe'

    print(os.path.exists(mf_exe_rtd))
    print(os.path.exists(mf_exe))
    print(os.path.exists(mp_exe))
    print(os.path.exists(mf_exe_git))

    package_folder_rtd = r"/home/docs/checkouts/readthedocs.org/user_builds/sutra2/checkouts/latest/sutra2"
    print(os.path.exists(package_folder_rtd))

    ## p = Path(package_folder_rtd).glob('**/*')
    ## files = [x for x in p if x.is_file()]
    ## print(files)


    ## mf_exe = r"/home/docs/checkouts/readthedocs.org/user_builds/sutra2/checkouts/latest/sutra2/mf2005.exe"
    ## mp_exe = r"/home/docs/checkouts/readthedocs.org/user_builds/sutra2/checkouts/latest/sutra2/mpath7.exe"

    # mf_exe = "../mf2005.exe"
    # mp_exe = "../mpath7.exe"
    mf_exe = "mf2005"  # r"d:\Sutra2_tool\sutra2\sutra2\mf2005.exe"
    mp_exe = "mpath7"  #r"d:\Sutra2_tool\sutra2\sutra2\mpath7.exe"

.. ipython:: python
    
    modpath_phrea = mpw.ModPathWell(phreatic_schematisation,
                                workspace = "phreatic_test",
                                modelname = "phreatic",
                                mf_exe = mf_exe, #"mf2005.exe",
                                mp_exe = mp_exe, #"mpath7.exe"
                                )



.. .. mf_exe = "..//mf2005.exe",
.. .. mp_exe = "..//mpath7.exe")

.. Now we run the Modpath model, which numerically calculates the flow in the subsurface using the 
.. 'schematisation' dictionary stored in the HydroChemicalSchematisation object. By default the model will
.. calculate both the hydraulic head distribution (using modflow: 'run_mfmodel' = True) and
.. the particle pathlines [X,Y,Z,T-data] (using modpath: 'run_mpmodel' = True) with which OMP removal
.. or microbial organism ('mbo') removal is later calculated.

.. ipython:: python

    modpath_phrea.run_model(run_mfmodel = False,
                        run_mpmodel = False)
    print(modpath_phrea.dstroot)
    modpath_phrea.dstroot

.. The traveltime distribution can be plotted as cross-section using either a linear or logarithmic distribution,
.. with lognorm = True: logarithmic distribution, using for example a 'viridis_r' (viridis reversed) color map.

.. .. ipython:: python

..     # time limits
..     tmin, tmax = 0.1, 10000.
..     # xcoord bounds
..     xmin, xmax = 0., 100.
..     # Create travel time plots (lognormal)
..     modpath_phrea.plot_pathtimes(df_particle = modpath_phrea.df_particle, 
..             vmin = tmin,vmax = tmax,
..             fpathfig = None, figtext = None,x_text = 0,
..             y_text = 0, lognorm = True, xmin = xmin, xmax = xmax,
..             line_dist = 1, dpi = 192, trackingdirection = "forward",
..             cmap = 'viridis_r')

.. .. fpath_plot = os.path.join(modpath_phrea.dstroot,"log_travel_times_test.png")
.. .. image: fpath_plot


.. From the ModpathWell class two other important outputs are:

.. * df_particle - Pandas dataframe with data about the different flowlines per particle node (vadose/shallow/target)
.. * df_flowline - Pandas dataframe with data about the flowlines per flowline (eg. total travel time per flowline)

.. Step 3: Collect removal parameters
.. ===========================================

.. Step 3a: View the Substance class (Optional)
.. ============================================
.. You can retrieve the default removal parameters used to calculate the removal of organic micropollutants [OMP] 
.. in the SubstanceTransport class. The data are stored in a dictionary

.. .. ipython:: python
    
..     test_substance = TR.Substance(substance_name='benzene')
..     test_substance.substance_dict

.. Step 3b: View the Organism class (Optional)
.. ===========================================
.. You can retrieve the default removal parameters used to calculate the removal of microbial organisms [mbo] 
.. in the SubstanceTransport class. The data are stored in a dictionary

.. .. ipython:: python
    
..     test_organism = TR.MicrobialOrganism(organism_name='MS2')
..     test_organism.organism_dict

.. Step 4: Run the SubstanceTransport class
.. ========================================
.. To calculate the removal and the steady-state concentration in each zone (analytical solution) or per particle node (modpath), create a concentration
.. object by running the SubstanceTransport class with the phreatic_well object and specifying
.. the OMP or microbial organism (mbo) of interest. 
.. The type of removal is defined using the option 'removal_function: 'omp' or 'mbo'
.. All required parameters for removal are stored as 'removal_parameters'.

.. In this example we use solani, which is a plant pathogen. First we create the object and view the organism properties:

.. .. ipython:: python

..     # Define removal parameters of microbial organism
..     organism_solani = TR.MicrobialOrganism(organism_name='solani')
..     # Connect to Transport class
..     phreatic_concentration = TR.Transport(modpath_phrea, pollutant = organism_solani)
..     phreatic_concentration.removal_parameters 

.. Optional: You may specify a different value for the removal_parameters, for example
.. a different inactivation rate 'mu1' or collission related removal 'alpha' and optional reference pH for 
.. calculating collision efficiency (pH0) for the anoxic redox zone while keeping other values as default.
.. This can be input in the SubstanceTransport object and this will be used in the calculation for 
.. the removal for the mbo.

.. .. ipython:: python

..     organism_solani_anox = TR.MicrobialOrganism(organism_name = 'solani',
..                                         alpha0_suboxic=None,
..                                         alpha0_anoxic=1.e-4,
..                                         alpha0_deeply_anoxic=None,
..                                         pH0_suboxic=None,
..                                         pH0_anoxic=7.5,
..                                         pH0_deeply_anoxic=None,
..                                         mu1_suboxic=None,
..                                         mu1_anoxic=0.01,
..                                         mu1_deeply_anoxic=None,)

..     phreatic_concentration = TR.Transport(modpath_phrea,pollutant = organism_solani_anox)

.. Step 4a: Calculate the removal of a (non-default) microbial organism ('mbo')
.. =============================================================================
.. In this example we calculate the removal of 'MS2' from a diffuse source, given 
.. that the modpath_model has completed successfully.

.. First we add removal parameters and create the 
.. SubstanceTransport object.

.. .. ipython:: python

..     # microbial removal properties of microbial organism
..     organism_name = 'MS2'
..     # reference_collision_efficiency [-]
..     alpha0 = {"suboxic": 1.e-3, "anoxic": 1.e-5, "deeply_anoxic": 1.e-5}
..     # reference pH for calculating collision efficiency [-]
..     pH0 = {"suboxic": 6.6, "anoxic": 6.8, "deeply_anoxic": 6.8}
..     # diameter of pathogen/species [m]
..     organism_diam =  2.33e-8
..     # inactivation coefficient [1/day]
..     mu1 = {"suboxic": 0.149,"anoxic": 0.023,"deeply_anoxic": 0.023}

..     # removal parameters for MS2 (manual input MicrobialOrganism)
..     organism_ms2 = TR.MicrobialOrganism(organism_name = organism_name,
..                                         alpha0_suboxic = alpha0["suboxic"],
..                                         alpha0_anoxic = alpha0["anoxic"],
..                                         alpha0_deeply_anoxic = alpha0["deeply_anoxic"],
..                                         pH0_suboxic = pH0["suboxic"],
..                                         pH0_anoxic = pH0["anoxic"],
..                                         pH0_deeply_anoxic = pH0["deeply_anoxic"],
..                                         mu1_suboxic = mu1["suboxic"],
..                                         mu1_anoxic = mu1["anoxic"],
..                                         mu1_deeply_anoxic = mu1["deeply_anoxic"],
..                                         organism_diam = organism_diam)

..     # Calculate advective microbial removal
..     modpath_removal = TR.Transport(modpath_phrea,
..                             pollutant = organism_ms2)
..     # Removal parameters organism
..     modpath_removal.removal_parameters

.. Then we calculate the final concentration after advective microbial removal of microbial organisms for a given endpoint_id
.. using the function 'calc_advective_microbial_removal'. This function calls a separate function 'calc_lambda'
.. which calculates the rate with which mbo's are removed per node along each given pathline. As input we use the
.. dataframes df_particle and df_flowline, which have been created by the ModpathWell class. These pandas dataframes
.. will be updated with calculated removal parameters and final_concentration per node. 
.. Also, we can plot the log removal along pathlines in a cross-section (optional)



.. .. ipython:: python
..     :okwarning:

..     C_final = {}
..     for endpoint_id in modpath_phrea.schematisation_dict.get("endpoint_id"):
..         df_particle, df_flowline, C_final[endpoint_id] = modpath_removal.calc_advective_microbial_removal(
..                                             modpath_phrea.df_particle, modpath_phrea.df_flowline, 
..                                             endpoint_id = endpoint_id,
..                                             conc_start = 1., conc_gw = 0.)
..         # relative conc limits
..         cmin, cmax = 1.e-11, 1.
..         # xcoord bounds
..         xmin, xmax = 0., 50.
..         # Create travel time plots (lognormal)
..         modpath_removal.plot_logremoval(df_particle=df_particle,
..                 df_flowline=df_flowline,
..                 vmin = cmin,vmax = cmax,
..                 fpathfig = None,
..                 y_text = 0, lognorm = True, xmin = xmin, xmax = xmax,
..                 trackingdirection = "forward",
..                 cmap = 'viridis_r')

.. Step 4b: Calculate the OMP removal
.. ========================================
.. Alternatively, you can calculate the removal of organic micropollutants (OMP). As example,
.. we take the default removal parameters for the substances 'AMPA'.
.. Note: For OMP you will have to specify values relevant for substances (e.g. half-life, pKa, log_Koc).
.. Any/all default values will be stored and used in the calculation of the removal. 
.. Note that by default the class expects the removal of microbial organisms copied from removal_function 
.. entered in modpath_phrea. We have to explicitly enter the removal_function below for removal op substances.
.. removal_function == 'omp'

.. .. ipython:: python

..     # substance (AMPA)
..     substance_name = 'AMPA'

..     # Load default removal parameters of AMPA
..     substance_ampa_default = TR.Substance(substance_name = substance_name,
..                                         partition_coefficient_water_organic_carbon=None,
..                                         dissociation_constant=None,
..                                         molar_mass = None,
..                                         halflife_suboxic=None,
..                                         halflife_anoxic=None,
..                                         halflife_deeply_anoxic=None
..                                         )
..     # Calculate removal of organic micro-pollutants (removal_function = 'omp')
..     modpath_removal = TR.Transport(well = modpath_phrea,
..                                     pollutant = substance_ampa_default
..                                     )

.. View the updated removal_parameters dictionary from the SubstanceTransport object

.. .. ipython:: python

..     modpath_removal.removal_parameters

.. We compute the removal by running the 'compute_omp_removal' function:
.. modpath_removal.compute_omp_removal()

.. .. ipython:: python
..     :okwarning:
    
..     modpath_removal.compute_omp_removal()


.. Once the removal has been calculated, you can view the steady-state concentration
.. and breakthrough time per zone for the OMP in the df_particle:

.. .. ipython:: python

..     phreatic_concentration.df_particle.loc[:,['zone', 'steady_state_concentration', 'travel_time']].head(4)

.. View the steady-state concentration of the flowline or the steady-state
.. contribution of the flowline to the concentration in the well

.. .. ipython:: python

..     phreatic_concentration.df_flowline.loc[:,['breakthrough_concentration', 'total_breakthrough_travel_time']].head(5)

.. .. Maak 'modpath' varianten voor de afbraak. Plots via jup nb

.. Plot the breakthrough curve at the well over time:

.. .. ipython:: python

..     benzene_plot = phreatic_concentration.plot_concentration(ylim=[0,10 ])

.. .. image:: https://github.com/KWR-Water/sutra2/blob/main/docs/_images/benzene_plot.png?raw=true
..   :width: 600
..   :alt: benzene_plot.png

.. You can also compute the removal for a different OMP of interest:

.. * OMP-X: a ficticous OMP with no degradation or sorption
.. * AMPA
.. * benzo(a)pyrene

.. To do so you can use the original schematisation, but specify a different OMP when you create
.. the SubstanceTransport object.

.. .. ipython:: python

..     # removal parameters OMP-X (default)
..     substance_ompx = TR.Substance(substance_name = "OMP-X")

..     phreatic_concentration = TR.Transport(modpath_phrea, pollutant = substance_ompx)
..     phreatic_concentration.compute_omp_removal()
..     omp_x_plot = phreatic_concentration.plot_concentration(ylim=[0,100 ])

.. .. image:: https://github.com/KWR-Water/sutra2/blob/main/docs/_images/omp_x_plot.png?raw=true
..   :width: 600
..   :alt: omp_x_plot.png

.. .. ipython:: python

..     # removal parameters benzo(a)pyrene (default)
..     substance_benzpy = TR.Substance(substance_name = "benzo(a)pyrene")
..     phreatic_concentration = TR.Transport(modpath_phrea, pollutant = substance_benzpy)
..     phreatic_concentration.compute_omp_removal()
..     benzo_plot = phreatic_concentration.plot_concentration(ylim=[0,1])

.. .. image:: https://github.com/KWR-Water/sutra2/blob/main/docs/_images/benzo_plot.png?raw=true
..   :width: 600
..   :alt: benzo_plot.png

.. .. ipython:: python

..     # removal parameters AMPA (default)
..     substance_ampa = TR.Substance(substance_name = "AMPA")
..     phreatic_concentration = TR.Transport(modpath_phrea, pollutant = substance_ampa)
..     phreatic_concentration.compute_omp_removal()
..     ampa_plot = phreatic_concentration.plot_concentration( ylim=[0,1])

.. .. image:: https://github.com/KWR-Water/sutra2/blob/main/docs/_images/ampa_plot.png?raw=true
..   :width: 600
..   :alt: ampa_plot.png

.. .. Other examples in the Bas_tutorial.py file are:

.. .. * diffuse/point source example for phreatic 
.. .. * semiconfined example



