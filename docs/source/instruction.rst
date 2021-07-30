========================================================================================================
Analytical solution 
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
solutions and to include pathogens in addition to OMP.

Steps
-----

Operating the analytical module typically involves 5 steps:

#. Define the hydrogeochemical system using the HydroChemicalSchematisation class. In the HydroChemicalSchematisation class choose the:

    * Computational method ('analytical' or 'modpath').
    * The schematisation type ('phreatic', 'semiconfined', 'riverbankfiltration', 'basinfiltration').
    * The removal function ('omp', 'pathogen).
    * Input the relevant geohydrological data

#. Run the AnalyticalWell class to calculate the travel time distribution in the different aquifer zones
#. Run the Substance class to retrieve the substance (removal) paramters
#. Run the SubstanceTransport class to calculate the removal and concentration in each zone and in the well
#. Plot the desired functions 

Now, let’s try some examples.

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

    from greta.Analytical_Well import *
    from greta.Substance_Transport import *


Step 1: Define the HydroChemicalSchematisation
===============================================
The first step is to define the hydrogeochemistry of the system using the HydroChemicalSchematisation class.
In this class you specify which solution you want (Analytical or Modflow), the 
schematisation from the choice of 4 PSWFs, and the relevant parameters for the porous 
media, the hydrochemistry, hydrology and the contamination of interest (OMP or 
pathogen). 

The class paramters can be roughly grouped into the following categories;

* System
* Settings
* Porous Medium
* Hydrochemistry
* Hydrology
* Contaminant
* Diffuse contamination
* Point Contamination
* Model size

Units of input ate 

* Time: days
* Length: meters
* Concentration: ug/L
* Temperature: degree C
* Depth: meters above sea level (m ASL)
* Density: kg/L
* DOC/TOC: mg/L

1A System parameters
--------------------------------------
Choose the schematisaiton type of the model from one of the four well types;

*'phreatic', 
*'semiconfined', 
*'riverbankfiltration', 
*'basinfiltration'

1B Settings
--------------------------------------

* computation_method
* removal_function
* temp_correction_Koc
* temp_correction_halflife
* biodegradation_sorbed_phase
* compute_thickness_vadose_zone

1C Porous Medium
--------------------------------------

1D Hydrochemistry
--------------------------------------

1E Hydrology
--------------------------------------

1F Contaminant
--------------------------------------
* Diffuse contamination
* Point Contamination

1G Modflow
--------------------------------------
Additional paramters about the model domain are input here

In this example we calculate the analytical solution for a phreatic well, with a diffuse 
contamination over the whole model domain.

.. ipython:: python
    
    from greta.Analytical_Well import HydroChemicalSchematisation
    phreatic_schematisation = HydroChemicalSchematisation(schematisation_type='phreatic',
                                        well_discharge=7500, #m3/day
                                        vertical_anistropy_shallow_aquifer=0.0006,
                                        porosity_vadose_zone=0.38,
                                        porosity_shallow_aquifer=0.35,
                                        porosity_target_aquifer=0.35,
                                        recharge_rate=0.00082, #m/day
                                        moisture_content_vadose_zone=0.15,
                                        ground_surface=22.,
                                        thickness_vadose_zone_at_boundary=5.,
                                        thickness_shallow_aquifer=10.,
                                        thickness_target_aquifer=40.,
                                        hor_permeability_target_aquifer=35.,
                                        thickness_full_capillary_fringe=0.4,
                                        redox_vadose_zone='suboxic',
                                        redox_shallow_aquifer='anoxic',
                                        redox_target_aquifer='deeply_anoxic',
                                        pH_vadose_zone=5.,
                                        pH_shallow_aquifer=6.,
                                        pH_target_aquifer=7.,
                                        dissolved_organic_carbon_vadose_zone=10., #mg/L
                                        dissolved_organic_carbon_shallow_aquifer=4., 
                                        dissolved_organic_carbon_target_aquifer=2.,
                                        fraction_organic_carbon_vadose_zone=0.001,
                                        fraction_organic_carbon_shallow_aquifer=0.0005,
                                        fraction_organic_carbon_target_aquifer=0.0005, 
                                        temperature=11.,
                                        solid_density_vadose_zone=2.650, 
                                        solid_density_shallow_aquifer=2.650, 
                                        solid_density_target_aquifer=2.650, 
                                        diameter_borehole=0.75,
                                        diffuse_input_concentration=600, #ug/L
                                        )

The paramters from the HydroChemicalSchematisation class are added as attributes to 
the class and can be accessed for example: 

.. .. ipython:: python
..     phreatic_schematisation.schematisation_type
..     phreatic_schematisation.well_discharge
..     phreatic_schematisation.porosity_shallow_aquifer

Step 2: Run the AnalyticalWell class 
=====================================
In the AnalyticalWell class the analytical solution for the chosen PSWF is run and 
the travel time is calculated for each of the zones. 

.. .. ipython:: python
..     phreatic_well = AnalyticalWell(phreatic_schematisation) # pass the phreatic_well object to initailize the well object
..     phreatic_well.phreatic() #calculate the travel time distribution for the phreatic analytical solution

You can plot the travel time distribution of the AnalyticalWell function here, as
well as the cumulative fraction of abstracted water

.. .. ipython:: python
..     phreatic_well.plot_travel_time_versus_radial_distance(xlim=[0, 2000], ylim=[1e3, 1e6])
..     phreatic_well.plot_travel_time_versus_cumulative_abstracted_water(xlim=[0, 1], ylim=[1e3, 1e6])

.. .. include:: 
..     travel_time_versus_cumulative_fraction_abstracted_water_phreatic.png
..     travel_time_versus_radial_distance_phreatic.png

From the AnalyticalWell class two important outputs are:

* df_particle - Pandas dataframe with the travel time per zone and 
* df_flowline


Step 3: View the Substance class (Optional)
===========================================
You can retrieve the default substance parameters used to calculate the removal in the 
SubstanceTransport class. 

.. .. ipython:: python
..     test_substance = Substance(substance_name='benzene')

You may specify a different value for the substance parameters, for example
a different half-life for the anoxic redox zone. This can be input in the HydroChemicalSchematisation
and this will be used in the calculation for the removal for the OMP.

.. .. ipython:: python
..     phreatic_schematisation = HydroChemicalSchematisation(schematisation_type='phreatic',
..                                 ....
..                                 halflife_anoxic= 650, )

Step 4: Run the SubstanceTransport class 
========================================
To calculate the removal and the steady-state concentration in each zone, create a concentration 
object by running the SubstanceTransport class with the phreatic_well object and specifying
the OMP (or pathogen) of interest. 

In this example we use benzene. First we create the object then compute the removal:

.. .. ipython:: python
..     phreatic_concentration = SubstanceTransport(phreatic_well, substance = 'benzene')
..     phreatic_concentration.compute_omp_removal()

If you have specified a values for the substance (e.g. half-life, pKa, log_Koc),
the default value is overriden and used in the calculation of the removal. You can
view the updated substance dictionary from the concentration object:

.. .. ipython:: python
..     phreatic_concentration.substance_dict

Once the removal has been calculated, you can view the steady-state concentration
and breakthrough time for the OMP in the df_particle:
.. .. ipython:: python
..     phreatic_concentration.df_particle['steady_state_concentration]
..     phreatic_concentration.df_particle['total_breakthrough_time']

View the steady-state concentration of the flowline or the steady-state 
contribution of the flowline to the concentration in the well
.. .. ipython:: python
..     phreatic_concentration.df_particle['breakthrough_concentration]
..     phreatic_concentration.df_particle['concentration_in_well']

Plot the breakthrough curve at the well over time
.. .. ipython:: python
..     phreatic_concentration.plot_concentration(xlim=[0, 500], ylim=[0,0.1 ])
.. .. include:: 
..     travel_time_versus_cumulative_fraction_abstracted_water_phreatic.png
..     travel_time_versus_radial_distance_phreatic.png

Other possibilities

* semiconfined example
* point source example



