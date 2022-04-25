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
-----

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

    from greta.Analytical_Well import *
    from greta.Substance_Transport import *



Step 1: Define the HydroChemicalSchematisation
==============================================
The first step is to define the hydrogeochemistry of the system using the HydroChemicalSchematisation class.
In this class you specify the:

    * Computational method ('analytical' or 'modpath').
    * The schematisation type ('phreatic', 'semiconfined').
    .. ('riverbankfiltration', 'basinfiltration' coming soon).
    * The removal function ('omp' or 'mbo').

    * Input the relevant parameters for the porous media, the hydrochemistry, hydrology and the contamination of interest

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
* Dischage : m3/d
* Time: days
* Length: meters
* Concentration: ug/L
* Temperature: degree C
* Depth: meters above sea level (m ASL)
* Density: kg/L
* DOC/TOC: mg/L

Lets start with a simple example defining a HydroChemicalSchematisation object for a phreatic aquifer:

.. ipython:: python

    phreatic_schematisation = HydroChemicalSchematisation(schematisation_type='phreatic',
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
                                                        temperature=11.,
                                                        substance='benzene',
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

Next we create a ModpathWell object for the HydroChemicalSchematisation object that we just made.

Step 2: Run the AnalyticalWell class
=====================================
Next we create an AnalyticalWell object for the HydroChemicalSchematisation object we just made.

.. ipython:: python

    phreatic_well = AnalyticalWell(phreatic_schematisation)

Then we calculate the travel time for each of the zones unsaturated, shallow aquifer and target aquifer zones
by running the .phreatic() function for the well object. 

.. ipython:: python

    phreatic_well.phreatic()

The total travel time can be plotted as a function of radial distance from the well, or as a function
of the cumulative fraction of abstracted water: 

.. ipython:: python

    radial_plot = phreatic_well.plot_travel_time_versus_radial_distance(xlim=[0, 2000], ylim=[1e3, 1e6])
    cumulative_plot = phreatic_well.plot_travel_time_versus_cumulative_abstracted_water(xlim=[0, 1], ylim=[1e3, 1e6])

.. image:: travel_time_versus_radial_distance_phreatic.png

.. image:: travel_time_versus_cumulative_abs_water_phreatic.png

From the AnalyticalWell class two other important outputs are:

* df_particle - Pandas dataframe with data about the different flowlines per zone (unsaturated/shallwo/target)
* df_flowline - Pandas dataframe with data about the flowlines per flowline (eg. total travel time per flowline)

Step 3: View the Substance class (Optional)
===========================================
You can retrieve the default substance parameters used to calculate the removal in the
SubstanceTransport class. The data are stored in a dictionary

.. ipython:: python
    
    test_substance = Substance(substance_name='benzene')
    test_substance.substance_dict


Step 4: Run the SubstanceTransport class
========================================
To calculate the removal and the steady-state concentration in each zone, create a concentration
object by running the SubstanceTransport class with the phreatic_well object and specifying
the OMP (or pathogen) of interest.

In this example we use benzene. First we create the object and view the substance properties:

.. ipython:: python

    phreatic_concentration = SubstanceTransport(phreatic_well, substance = 'benzene')
    phreatic_concentration.substance_dict

Optional: You may specify a different value for the substance parameters, for example
a different half-life for the anoxic redox zone. This can be input in the HydroChemicalSchematisation
and this will be used in the calculation for the removal for the OMP. The AnalyticalWell and 
phreatic() functions must be rerun:

.. ipython:: python

    phreatic_schematisation = HydroChemicalSchematisation(schematisation_type='phreatic',
                                                            well_discharge=-7500, #m3/day
                                                            recharge_rate=0.0008, #m/day
                                                            thickness_vadose_zone_at_boundary=5,
                                                            thickness_shallow_aquifer=10,
                                                            thickness_target_aquifer=40,
                                                            hor_permeability_target_aquifer=35,
                                                            redox_vadose_zone='anoxic',
                                                            redox_shallow_aquifer='anoxic',
                                                            redox_target_aquifer='deeply_anoxic',
                                                            pH_target_aquifer=7.,
                                                            temperature=11.,
                                                            substance='benzene',
                                                            diffuse_input_concentration = 100, #ug/L
                                                            partition_coefficient_water_organic_carbon=2,
                                                            dissociation_constant=1,
                                                            halflife_suboxic=12, 
                                                            halflife_anoxic=420, 
                                                            halflife_deeply_anoxic=6000,
                                                            )
    phreatic_well = AnalyticalWell(phreatic_schematisation)
    phreatic_well.phreatic() 
    phreatic_concentration = SubstanceTransport(phreatic_well, substance = 'benzene')
    
If you have specified a values for the substance (e.g. half-life, pKa, log_Koc),
the default value is overriden and used in the calculation of the removal. You can
view the updated substance dictionary from the concentration object:

.. ipython:: python

    phreatic_concentration.substance_dict

Then we compute the removal by running the 'compute_omp_removal' function:
phreatic_concentration.compute_omp_removal()

.. ipython:: python
    
    phreatic_concentration.compute_omp_removal()


Once the removal has been calculated, you can view the steady-state concentration
and breakthrough time per zone for the OMP in the df_particle:

.. ipython:: python

    phreatic_concentration.df_particle[['flowline_id', 'zone', 'steady_state_concentration', 'breakthrough_travel_time']].head(4)

View the steady-state concentration of the flowline or the steady-state
contribution of the flowline to the concentration in the well

.. ipython:: python

    phreatic_concentration.df_flowline[['flowline_id', 'breakthrough_concentration', 'total_breakthrough_travel_time']].head(5)

Plot the breakthrough curve at the well over time:

.. ipython:: python

    benzene_plot = phreatic_concentration.plot_concentration(ylim=[0,10 ])

.. image:: benzene_plot.png

You can also compute the removal for a different OMP of interest:

* OMP-X: a ficticous OMP with no degradation or sorption
* AMPA
* benzo(a)pyrene

To do so you can use the original schematisation, but specify a different OMP when you create
the SubstanceTransport object.

.. ipython:: python

    phreatic_well = AnalyticalWell(phreatic_schematisation)
    phreatic_well.phreatic() 
    phreatic_concentration = SubstanceTransport(phreatic_well, substance = 'OMP-X')
    phreatic_concentration.compute_omp_removal()
    omp_x_plot = phreatic_concentration.plot_concentration(ylim=[0,100 ])

.. image:: omp_x_plot.png


.. ipython:: python

    phreatic_well = AnalyticalWell(phreatic_schematisation)
    phreatic_well.phreatic() 
    phreatic_concentration = SubstanceTransport(phreatic_well, substance = 'benzo(a)pyrene')
    phreatic_concentration.compute_omp_removal()
    benzo_plot = phreatic_concentration.plot_concentration(ylim=[0,1])

.. image:: benzo_plot.png

.. ipython:: python

    phreatic_well = AnalyticalWell(phreatic_schematisation)
    phreatic_well.phreatic() 
    phreatic_concentration = SubstanceTransport(phreatic_well, substance = 'AMPA')
    phreatic_concentration.compute_omp_removal()
    ampa_plot = phreatic_concentration.plot_concentration( ylim=[0,1])

.. image:: ampa_plot.png


Other examples in the Bas_tutorial.py file are:

* diffuse/point source example for phreatic 
* semiconfined example



