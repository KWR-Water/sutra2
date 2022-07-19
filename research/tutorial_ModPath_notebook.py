# <markdowncell>
# ## Modpath tutorial notebook


# <codecell>
# set working directory to main directory of package
import codecs
import imp

from matplotlib.pyplot import ylabel
from set_cwd_to_project_root import project_root
from pathlib import Path

from pandas import RangeIndex, read_csv
import pandas as pd
from pandas._testing import assert_frame_equal
import numpy as np
import datetime as dt
import os
import ast  # abstract syntax trees
import sys
import copy
# path = os.getcwd()  # path of working directory
from pathlib import Path

from zmq import zmq_version_info

# Plotting modules
import matplotlib.pyplot as plt
import matplotlib 
import matplotlib.colors as colors

# sutra2 modules    
import sutra2.Analytical_Well as AW
import sutra2.ModPath_Well as mpw
import sutra2.Transport_Removal as TR

# get directory of this file
path = Path(__file__).parent
%load_ext autoreload
%autoreload 2

# Create research files dir
researchfiles_dir = os.path.join(path,"test_files")
if not os.path.exists(researchfiles_dir):
    os.makedirs(researchfiles_dir)

# modflow executable location
mf_exe = os.path.join(path, "mf2005.exe")
# modpath executable location
mp_exe = os.path.join(path, "mpath7.exe")

# <markdowncell>
# Test document for tutorial_modpath.rst
# Create plots and other data to show in html-format (Read The Docs)


# <codecell>

''' calculate the removal of a default substance using ModPathWell class.'''
# Lets start with a simple example defining a HydroChemicalSchematisation object for a phreatic aquifer:
phreatic_schematisation = AW.HydroChemicalSchematisation(schematisation_type='phreatic',
                                                    computation_method = 'modpath',
                                                    removal_function = 'omp',
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

# <markdowncell>
# Then, we create a ModpathWell object for the HydroChemicalSchematisation object that we just made.
# The ModpathWell object requires a dictionary of the subsurface schematisation and a set of boundary conditions
# the numerical model has to abide by in calculating flow velocity and direction of flow.

# <codecell>
phreatic_schematisation.make_dictionary()

# <markdowncell>
# Step 2: Run the ModpathWell class
# =====================================
# Next we create an ModpathWell object for the HydroChemicalSchematisation object we just made.
# The data files will be stored in location workspace using a given modelname.

#<codecell>
modpath_phrea = mpw.ModPathWell(phreatic_schematisation,
                            workspace = os.path.join(path,"test7_omp_removal"),
                            modelname = "phreatic")

# <markdowncell>
# Now we run the Modpath model, which numerically calculates the flow in the subsurface using the 
# 'schematisation' dictionary stored in the HydroChemicalSchematisation object. By default the model will
# calculate both the hydraulic head distribution (using modflow: 'run_mfmodel' = True) and
# the particle pathlines [X,Y,Z,T-data] (using modpath: 'run_mpmodel' = True) with which OMP removal
# or microbial organism ('mbo') removal is later calculated.

# <codecell>
modpath_phrea.run_model(run_mfmodel = True,
                    run_mpmodel = True)

# <markdowncell>
# Step 3: Collect removal parameters
# ===========================================

# Step 3a: View the Substance class (Optional)
# ============================================
# You can retrieve the default removal parameters used to calculate the removal of organic micropollutants [OMP] 
# in the SubstanceTransport class. The data are stored in a dictionary

# <codecell>
test_substance = TR.Substance(substance_name='benzene')
test_substance.substance_dict

#<markdowncell>
# Step 4: Run the SubstanceTransport class
# ========================================
# To calculate the removal and the steady-state concentration in each zone (analytical solution) or per particle node (modpath), create a concentration
# object by running the SubstanceTransport class with the phreatic_well object and specifying
# the OMP or microbial organism (mbo) of interest. 
# The type of removal is defined using the option 'removal_function: 'omp' or 'mbo'
# All required parameters for removal are stored as 'removal_parameters'.

# Step 4b: Calculate the OMP removal
# ========================================
# As example, we take the default removal parameters for the substances 'AMPA'.
# Note: For OMP you will have to specify values relevant for substances (e.g. half-life, pKa, log_Koc).
# Any/all default values will be stored and used in the calculation of the removal. 
# Note that by default the class expects the removal of microbial organisms copied from removal_function 
# entered in modpath_phrea. We have to explicitly enter the removal_function below for removal op substances.
# removal_function == 'omp'

# <codecell>
# substance (AMPA)
substance_name = 'AMPA'
# Calculate removal of organic micro-pollutants (removal_function = 'omp')
modpath_removal = TR.Transport(well = modpath_phrea,
                        substance = substance_name,
                        partition_coefficient_water_organic_carbon=None,
                        dissociation_constant=None,
                        halflife_suboxic=None,
                        halflife_anoxic=None,
                        halflife_deeply_anoxic=None,
                        removal_function = 'omp',
                        )

# <markdowncell>
# View the updated removal_parameters dictionary from the SubstanceTransport object

# <codecell>
modpath_removal.removal_parameters

# <markdowncell>
# We compute the removal by running the 'compute_omp_removal' function:
# modpath_removal.compute_omp_removal()

# <codecell>
modpath_removal.compute_omp_removal()

# <markdowncell>
# Once the removal has been calculated, you can view the steady-state concentration
# and breakthrough time per zone for the OMP in the df_particle:

# <codecell>
modpath_removal.df_particle.loc[:,['zone', 'steady_state_concentration', 'travel_time']].head(4)

# <markdowncell>
# View the steady-state concentration of the flowline or the steady-state
# contribution of the flowline to the concentration in the well

# <codecell>
modpath_removal.df_flowline.loc[:,['breakthrough_concentration', 'total_breakthrough_travel_time']].head(5)

# <markdowncell>
# Maak 'modpath' varianten voor de afbraak. Plots via jup nb
# 
# Plot the breakthrough curve at the well over time

# <codecell>
benzene_plot = phreatic_concentration.plot_concentration(ylim=[0,10 ])

# .. image:: benzene_plot.png

# <markdowncell>
# You can also compute the removal for a different OMP of interest:
# 
# * OMP-X: a ficticous OMP with no degradation or sorption
# * AMPA
# * benzo(a)pyrene
# 
# To do so you can use the original schematisation, but specify a different OMP when you create
# the SubstanceTransport object.
#
#.. ipython:: python

    phreatic_concentration = SubstanceTransport(phreatic_well, substance = 'OMP-X')
    phreatic_concentration.compute_omp_removal()
    omp_x_plot = phreatic_concentration.plot_concentration(ylim=[0,100 ])

.. image:: omp_x_plot.png

.. ipython:: python

    phreatic_concentration = SubstanceTransport(phreatic_well, substance = 'benzo(a)pyrene')
    phreatic_concentration.compute_omp_removal()
    benzo_plot = phreatic_concentration.plot_concentration(ylim=[0,1])

.. image:: benzo_plot.png

.. ipython:: python

    phreatic_concentration = SubstanceTransport(phreatic_well, substance = 'AMPA')
    phreatic_concentration.compute_omp_removal()
    ampa_plot = phreatic_concentration.plot_concentration( ylim=[0,1])

.. image:: ampa_plot.png

.. Other examples in the Bas_tutorial.py file are:

.. * diffuse/point source example for phreatic 
.. * semiconfined example

#%%
# if __name__ == "__main__":
#     test_modpath_run_phreatic_withgravelpack()