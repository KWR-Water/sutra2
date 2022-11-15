#%% ----------------------------------------------------------------------------
# A. Hockin, March 2021
# KWR BO 402045-247
# ZZS verwijdering bodempassage
# AquaPriori - Transport Model
# With Martin Korevaar, Martin vd Schans, Steven Ros, Bas vd Grift
#
# Based on Stuyfzand, P. J. (2020). Predicting organic micropollutant behavior
#               for 4 public supply well field types, with TRANSATOMIC Lite+
#               (Vol. 2). Nieuwegein, Netherlands.
# ------------------------------------------------------------------------------

#### Notes ####

# things which must be checked indicated in comments with AH
# specific questions flagged for;
# @MartinvdS // @steven //@martinK

#To Do marked as ah_todo

####

#%% ----------------------------------------------------------------------------
# INITIALISATION OF PYTHON e.g. packages, etc.
# ------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from pandas import read_csv
from pandas import read_excel
import math
from scipy.special import kn as besselk
import datetime as dt
from datetime import timedelta
from datetime import datetime
import warnings

path = os.getcwd()  # path of working directory


# @Steven_todo --> create_schematisation_dictionary (dict als uitvoer)
# Maak er een functie van wat dict oplevert ()

# @Steven_todo: HydroChemicalSchematisation verplaatsen naar andere python-file.
# Daarna deze file aanroepen (Analytical_Well.py); idem voor imports in ModPath_Well.py
# from sutra2.HydroChemicalSchematisation import *


# @ Steven_todo: functie --> run everything
# EEN WRAPPER! --> private functie aanmaken "_aquapriori_wrapper" Analytical_Well of ModPath_well als argument link met OMP removal
# Overige argumenten als verwijzing naar klasse

# @Steven_todo: functie voor ModPath --> Alles voor Analytical_Well klasse eruit halen
# @steven_todo @mvds alleen nog mar de properties in een dict op slaan, niets meer in de attributeten van deze
# klassse. hiervoor moet er dus wel veel worden gewijzigd (al is het eigenlijk alleen hernoemen)


def create_schematisation_dictionary():

    """ This function returns a dictionary used by the ModPath simulation
        holding all the parameters that define the (surrounding of the) well such that travel
        times and removal rates can be determined.

    Parameters
    ----------
    schematisation_type: string
        Choose the schematisation of the aquifer, choice of 4 public supply well fields
        'phreatic', 'semiconfined', 'riverbankfiltration', 'basinfiltration'
    well_name: string
        Name of the well in the simulaiton. Default is 'well1'
    ground_surface: float
        Meters above sea level (ASL)
    thickness_vadose_zone_at_boundary, thickness_shallow_aquifer, thickness_target_aquife: float
        Thickness of each of the aquifer zones of interest (vadose, shallow and target), [m].
    thickness_full_capillary_fringe: float
        Thickness of the capillary fringe (subsurface layer in which groundwater seeps up from a water table by capillary action to fill pores), [m].
    porosity_vadose_zone, porosity_shallow_aquifer, porosity_target_aquifer, porosity_gravelpack, porosity_clayseal: float
        Porosity of each of the zones of interest, dimensionless [volume fraction].
    moisture_content_vadose_zone: float
        Mean multi-annual moisture content of unsaturated zone, dimensionless [volume fraction].
    solid_density_vadose_zone, solid_density_shallow_aquifer, solid_density_target_aquifer: float
        Solid density, in kg/L, of each of the aquifer zones of interest.
    fraction_organic_carbon_vadose_zone, fraction_organic_carbon_shallow_aquifer, fraction_organic_carbon_target_aquifer: float
        Mass fraction of the organic carbon in the aquifer zone of interest.
    grainsize_vadose_zone, grainsize_shallow_aquifer, grainsize_target_aquifer: float
        Grain diameter (d50) of the sediment, [m].
    redox_vadose_zone, redox_shallow_aquifer, redox_target_aquifer: string
            Definition of redox zones follows the scheme presented by Stuyfzand (1993, 2012b), which is in simplified form:
            'suboxic': NO3 ≥ 1 mg/L, Fe = Mn < 0.1 mg/L;
            'anoxic': O2 and NO3 < 1 mg/L, no SO4-reduction;
            'deeply_anoxic': = O2 and NO3 < 1 mg/L, with SO4-reduction and/or methanogenesis.
    dissolved_organic_carbon_vadose_zone, dissolved_organic_carbon_shallow_aquifer, dissolved_organic_carbon_target_aquifer: float
        Concentration of the dissolved organic carbon in each of the aquifer zones, [mg/L].
    dissolved_organic_carbon_infiltration_water: float
        Concentrtaion of the dissolved organic carbon in the infiltrating water, [mg/L].
    total_organic_carbon_infiltration_water: float
        Concentrtaion of the total organic carbon in the infiltrating water, [mg/L].
    pH_vadose_zone, pH_shallow_aquifer, pH_target_aquifer: float
        pH in each aquifer zone, dimensionless.
    temp_water, temperature_vadose_zone, temperature_shallow_aquifer, temperature_target_aquifer: float
        Temperature, in degrees Celcius, for each aquifer zone. If only 'temperature' is given, the temperature in each zone is set to 'temp_water'.
    recharge_rate: float
        Mean multi-annual groundwater recharge rate, [m/d].
    well_discharge: float
        Mean annual well discharge, [m3/d].

    diameter_borehole: float
        Diameter of the borehole of the well, [m].
    top_filterscreen, bottom_filterscreen: float
        Top and bottom of the filterscreen of the well [mASL].
    diameter_filterscreen, inner_diameter_filterscreen: float
        Diameter and inner diameter of the filterscreen of the well, [m].
    top_gravelpack,  bottom_gravelpack: float
        Top and bottom of the gravelpack of the well [mASL].
    diameter_gravelpack, inner_diameter_gravelpack: float
        Diameter and inner diameter of the gravelscreen of the well, [m].
    grainsize_gravelpack: float
        Grain diameter (d50) of the gravelpack, [m].
    top_clayseal, bottom_clayseal: float
        Top and bottom of the clayseal of the well [mASL].
    diameter_clayseal, inner_diameter_clayseal: float
        Diameter and inner diameter of the clayseal of the well, [m].
    grainsize_clayseal: float
        Grain diameter (d50) of the clay seal, [m].

    hor_permeability_shallow_aquifer, hor_permeability_target_aquifer, hor_permebility_gravelpack, hor_permeability_clayseal: float
        Horizonatal permeability of each zone of the aquifer and the clayseal, [m/d].
    hor_permeability_clayseal: float
        Horizonatal permeability of the clayseal, [m/d].
    vertical_anisotropy_shallow_aquifer, vertical_anisotropy_target_aquifer, vertical_anisotropy_gravelpack, vertical_anisotropy_clayseal: float
        Vertical anisotropy of each zone of the aquifer and the clayseal, [m/d].
        @MartinvdS  ratio of horizontal to vertical hydraulic conductivity?
        #AH MartinvdS error in this?
    name: string
        Name of the OMP of MBO to use in the model.
    diffuse_input_concentration: float
        Concentration of the diffuse source of the OMP in the groundwater recharge, [ug/L].
    point_input_concentration: float
        Concentration of the point source contamintation [ug/L].
    distance_point_contamination_from_well, depth_point_contamination, discharge_point_contamination: float
        Distance [m] from the well, depth [mASL] and discharge [m3/d] of the point source contamination.

    source_area_radius: float
        Radius of the source area, [m].

    number_of_spreading_distance: float
        the higher the number of spreading distance between well and boundary, the smaller any boundary effects will be. (Ignored in this version)
    model_radius, model_width: float
        Radius and width of the model, [m].
    ncols_filterscreen, ncols_gravelpack, ncols_near_well, ncols_far_well: int
        Number of model columns in the filterscreen, gravelpack, near and far from the well.
    delr_filterscreen, delr_gravelpack, delr_near_well, delr_far_well: float
        delr is the column spacing in the row direction for each of the filerscreen, gravelpak, near and far from the well.
    delc: list of floats
        delc is the row spacing in the column direction.
    nlayer_shallow_aquifer, nlayer_target_aquifer:
        The number of layers in the model grid for the shallow and target aquifers.
    relative_position_starting_points_radial, relative_position_starting_points_in_basin, relative_position_starting_points_outside_basin: list of floats -
        not yet included functionality for relative positions starting points flowlines

    BAR specific parameters: #AH_todo, will define later
        basin_length: float m
        basin_width: float m
        basin_xmin: float m
        basin_xmax: float m
        bottom_basin: float mASL
        bottom_sludge: float mASL
        hor_permeability_sludge: float m/d
        vertical_anisotropy_sludge: float -
        head_basin: float mASL
        basin_infiltration_rate: float m3/d
        travel_time_h20_shallow_aquifer: float d
        minimum_travel_time_h20_shallow_aquifer: float d
        travel_time_h20_deeper_aquifer: float d
        minimum_travel_time_h20_target_aquifer: float d

    """
    # initialize non-defaults to None to start.
    def __init__(self,
                schematisation_type,

                computation_method='analytical',
                removal_function = 'omp',
                name = None,  # substance or species name
                what_to_export = 'all', # MK: cryptic, AH_todo w/MartinK - better name for this
                temp_correction_Koc=True,
                temp_correction_halflife=True,
                biodegradation_sorbed_phase=True,
                compute_thickness_vadose_zone=True,
                well_name='well1',

                ground_surface=0.0,
                thickness_vadose_zone_at_boundary=1,
                # head_boundary=None, @MartinvdS is this used for anything? In the list of params but not used elsewhere...

                thickness_shallow_aquifer=10.0,
                thickness_target_aquifer=10.0,
                thickness_full_capillary_fringe=0.0,
                porosity_vadose_zone=0.35,
                porosity_shallow_aquifer=0.35,
                porosity_target_aquifer=0.35,
                porosity_gravelpack=0.35,
                porosity_clayseal=0.35,
                moisture_content_vadose_zone=0.2,
                solid_density_vadose_zone=2.65,
                solid_density_shallow_aquifer=2.65,
                solid_density_target_aquifer=2.65,
                fraction_organic_carbon_vadose_zone=0.001,
                fraction_organic_carbon_shallow_aquifer=0.0005,
                fraction_organic_carbon_target_aquifer=0.0005,
                grainsize_vadose_zone= 0.00025,
                grainsize_shallow_aquifer=0.00025,
                grainsize_target_aquifer=0.00025,

                redox_vadose_zone='suboxic',
                redox_shallow_aquifer='anoxic',
                redox_target_aquifer='deeply_anoxic',
                dissolved_organic_carbon_vadose_zone=0.0,
                dissolved_organic_carbon_shallow_aquifer=0.0,
                dissolved_organic_carbon_target_aquifer=0.0,
                dissolved_organic_carbon_infiltration_water=0.0,
                total_organic_carbon_infiltration_water=0.0,
                pH_vadose_zone=7.0,
                pH_shallow_aquifer=7.0,
                pH_target_aquifer=7.0,
                temp_water=10.0,
                temperature_vadose_zone=None,
                temperature_shallow_aquifer=None,
                temperature_target_aquifer=None,

                recharge_rate=0.001,
                well_discharge=-1000.0,

                basin_length=None, # BAR parameters \/
                basin_width=None,
                basin_xmin=None,
                basin_xmax=None,
                bottom_basin=None,
                bottom_sludge=None,
                hor_permeability_sludge=None,
                vertical_anisotropy_sludge=None,
                head_basin=None,
                basin_infiltration_rate=None,
                travel_time_h20_shallow_aquifer=None,
                minimum_travel_time_h20_shallow_aquifer=None,
                travel_time_h20_deeper_aquifer=None,
                minimum_travel_time_h20_target_aquifer=None, # BAR parameters /\

                diameter_borehole=0.75,
                top_filterscreen=None,
                bottom_filterscreen=None,
                diameter_filterscreen=0.2,
                inner_diameter_filterscreen=None,
                top_gravelpack=None,
                bottom_gravelpack=None,
                diameter_gravelpack=None,
                inner_diameter_gravelpack=None,
                grainsize_gravelpack=0.001,
                top_clayseal=None,
                bottom_clayseal=None,
                diameter_clayseal=None,
                inner_diameter_clayseal=None,
                grainsize_clayseal=0.000001,
                hor_permeability_shallow_aquifer=1.0,
                hor_permeability_target_aquifer=10.0,
                hor_permeability_gravelpack=None,
                hor_permeability_clayseal=None,
                vertical_anisotropy_shallow_aquifer=1.0,
                vertical_anisotropy_target_aquifer=1.0,
                vertical_anisotropy_gravelpack=None,
                vertical_anisotropy_clayseal=None,

                # replaced by name: depends on removal_function if it is an microbial organism ['mbo'] or substance ['omp']
                # substance=None,

                diffuse_input_concentration=1,

                point_input_concentration=None,
                distance_point_contamination_from_well=np.array([1]),
                depth_point_contamination=None,
                discharge_point_contamination=None,

                source_area_radius=None,
                number_of_spreading_distance=None,
                model_radius=None,
                model_width=1, #AH @MartinvdS for the ycoord? cell width, see def _export_to_df(self, ):
                relative_position_starting_points_radial=None,
                relative_position_starting_points_in_basin=None,
                relative_position_starting_points_outside_basin=None,

                particle_release_day=None,

                #modpath params
                 ncols_near_well = 20,
                 ncols_far_well = 30,
                 nlayers_shallow_aquifer = 10,
                 nlayers_target_aquifer = 10,
                 ):

        ''' Assigns the parameters to be attributes of the class and sets the default values'''

        # System
        self.schematisation_type = schematisation_type

        # Settings
        self.computation_method = computation_method
        self.removal_function = removal_function
        self.what_to_export = what_to_export
        self.well_name = well_name

        def check_parameter_choice(parameters_choice, options ):
            '''Check if parameter_choice is one of thelisted options,
            if not raise error'''
            for var in parameters_choice:
                value = getattr(self, var)
                if value not in options:
                    raise ValueError(f'Invalid {var}. Expected one of: {options}')

        check_parameter_choice(parameters_choice = ['schematisation_type'], options =['phreatic', 'semiconfined',])
        check_parameter_choice(parameters_choice = ['computation_method'], options =['analytical', 'modpath',])
        check_parameter_choice(parameters_choice = ['removal_function'], options =['omp','mbo',])
        check_parameter_choice(parameters_choice = ['what_to_export'], options =['all','omp', 'mbo',])

        self.temp_correction_Koc = temp_correction_Koc
        self.temp_correction_halflife = temp_correction_halflife
        self.biodegradation_sorbed_phase = biodegradation_sorbed_phase
        self.compute_thickness_vadose_zone = compute_thickness_vadose_zone

        # Porous Medium
        self.ground_surface = ground_surface #as meters above sea level (m ASL)
        self.thickness_vadose_zone_at_boundary = abs(thickness_vadose_zone_at_boundary)
        self.bottom_vadose_zone_at_boundary = ground_surface - abs(thickness_vadose_zone_at_boundary)
        # self.head_boundary = head_boundary #AH @MartinvdS is this used anywhere? I don't think so...

        # check that thicknesses for shallow and target aquifer > 0
        if thickness_shallow_aquifer <=0:
            raise ValueError('Error, thickness_shallow_aquifer should be > 0.')
        if thickness_target_aquifer <=0:
            raise ValueError('Error, thickness_target_aquifer should be > 0.')

        self.thickness_shallow_aquifer = abs(thickness_shallow_aquifer)
        self.bottom_shallow_aquifer = ground_surface - abs(thickness_vadose_zone_at_boundary) - abs(thickness_shallow_aquifer)
        self.thickness_target_aquifer = abs(thickness_target_aquifer)
        self.bottom_target_aquifer = self.bottom_shallow_aquifer - abs(thickness_target_aquifer)
        self.thickness_full_capillary_fringe = abs(thickness_full_capillary_fringe)
        self.porosity_vadose_zone = porosity_vadose_zone
        self.porosity_shallow_aquifer = porosity_shallow_aquifer
        self.porosity_target_aquifer = porosity_target_aquifer
        self.porosity_gravelpack = porosity_gravelpack
        self.porosity_clayseal = porosity_clayseal
        self.moisture_content_vadose_zone = moisture_content_vadose_zone
        self.solid_density_vadose_zone = solid_density_vadose_zone
        self.solid_density_shallow_aquifer = solid_density_shallow_aquifer
        self.solid_density_target_aquifer = solid_density_target_aquifer
        self.fraction_organic_carbon_vadose_zone = fraction_organic_carbon_vadose_zone
        self.fraction_organic_carbon_shallow_aquifer = fraction_organic_carbon_shallow_aquifer
        self.fraction_organic_carbon_target_aquifer = fraction_organic_carbon_target_aquifer
        self.grainsize_vadose_zone = grainsize_vadose_zone
        self.grainsize_shallow_aquifer = grainsize_shallow_aquifer
        self.grainsize_target_aquifer = grainsize_target_aquifer

        # Hydrochemistry
        self.redox_vadose_zone = redox_vadose_zone
        self.redox_shallow_aquifer = redox_shallow_aquifer
        self.redox_target_aquifer = redox_target_aquifer

        #  Check redox zone options, if not one listed, raise error
        check_parameter_choice(parameters_choice = ['redox_vadose_zone'], options =['suboxic', 'anoxic', 'deeply_anoxic'])
        check_parameter_choice(parameters_choice = ['redox_shallow_aquifer'], options =['suboxic', 'anoxic', 'deeply_anoxic'])
        check_parameter_choice(parameters_choice = ['redox_target_aquifer'], options =['suboxic', 'anoxic', 'deeply_anoxic'])

        self.dissolved_organic_carbon_vadose_zone = dissolved_organic_carbon_vadose_zone
        self.dissolved_organic_carbon_shallow_aquifer = dissolved_organic_carbon_shallow_aquifer
        self.dissolved_organic_carbon_target_aquifer = dissolved_organic_carbon_target_aquifer
        self.dissolved_organic_carbon_infiltration_water = dissolved_organic_carbon_infiltration_water
        self.total_organic_carbon_infiltration_water = total_organic_carbon_infiltration_water
        self.pH_vadose_zone = pH_vadose_zone
        self.pH_shallow_aquifer = pH_shallow_aquifer
        self.pH_target_aquifer = pH_target_aquifer
        self.temp_water = temp_water
        self.temperature_vadose_zone = temperature_vadose_zone
        self.temperature_shallow_aquifer = temperature_shallow_aquifer
        self.temperature_target_aquifer = temperature_target_aquifer

        # Hydrology
        self.recharge_rate = recharge_rate
        self.well_discharge = well_discharge
        self.basin_length = basin_length
        self.basin_width = basin_width
        self.basin_xmin = basin_xmin
        self.basin_xmax = basin_xmax
        self.bottom_basin = bottom_basin
        self.bottom_sludge = bottom_sludge
        self.hor_permeability_sludge = hor_permeability_sludge
        self.vertical_anisotropy_sludge = vertical_anisotropy_sludge
        self.head_basin = head_basin
        self.basin_infiltration_rate = basin_infiltration_rate
        self.travel_time_h20_shallow_aquifer = travel_time_h20_shallow_aquifer
        self.minimum_travel_time_h20_shallow_aquifer = minimum_travel_time_h20_shallow_aquifer
        self.travel_time_h20_deeper_aquifer = travel_time_h20_deeper_aquifer
        self.minimum_travel_time_h20_target_aquifer = minimum_travel_time_h20_target_aquifer
        self.diameter_borehole = diameter_borehole
        self.top_filterscreen = top_filterscreen
        self.bottom_filterscreen = bottom_filterscreen
        self.diameter_filterscreen = diameter_filterscreen
        self.inner_diameter_filterscreen = inner_diameter_filterscreen
        self.top_gravelpack = top_gravelpack
        self.bottom_gravelpack = bottom_gravelpack
        self.diameter_gravelpack = diameter_gravelpack
        self.inner_diameter_gravelpack = inner_diameter_gravelpack
        self.grainsize_gravelpack = grainsize_gravelpack
        self.top_clayseal = top_clayseal
        self.bottom_clayseal = bottom_clayseal
        self.diameter_clayseal = diameter_clayseal
        self.inner_diameter_clayseal = inner_diameter_clayseal
        self.grainsize_clayseal = grainsize_clayseal

        # do we need this?
        self.hor_permeability_shallow_aquifer = hor_permeability_shallow_aquifer
        self.hor_permeability_target_aquifer = hor_permeability_target_aquifer
        self.hor_permeability_gravelpack = hor_permeability_gravelpack
        self.hor_permeability_clayseal = hor_permeability_clayseal

        # @MartinvdS vertical anisotropy not used?
        self.vertical_anisotropy_shallow_aquifer = vertical_anisotropy_shallow_aquifer
        self.vertical_anisotropy_target_aquifer = vertical_anisotropy_target_aquifer
        self.vertical_anisotropy_gravelpack = vertical_anisotropy_gravelpack
        self.vertical_anisotropy_clayseal = vertical_anisotropy_clayseal

        #Calculated
        self.vertical_permeability_target_aquifer = self.hor_permeability_target_aquifer / self.vertical_anisotropy_target_aquifer
        self.vertical_permeability_shallow_aquifer = self.hor_permeability_shallow_aquifer / self.vertical_anisotropy_shallow_aquifer

        self.vertical_resistance_shallow_aquifer = thickness_shallow_aquifer / self.vertical_permeability_shallow_aquifer

        self.KD = hor_permeability_target_aquifer * thickness_target_aquifer

        self.groundwater_level =self.ground_surface-self.thickness_vadose_zone_at_boundary

        self.organism = None
        self.substance = None
        if self.removal_function in ['omp',]:
            # Substance name
            self.substance = name
        elif self.removal_function in ['mbo',]:
            # Organism name
            self.organism = name


        # Diffuse contamination override if point contamination specified
        self.diffuse_input_concentration = diffuse_input_concentration
        self.point_input_concentration = point_input_concentration

        def check_date_format(check_date):
            for var in check_date:
                value = getattr(self, var)
                if not type(value) is dt.datetime:
                    raise TypeError(f"Error invalid date input, please enter a new {var} using the format dt.datetime.strptime('YYYY-MM-DD', '%Y-%m-%d')")

        self.start_date_well = start_date_well

        check_date_format(check_date=['start_date_well'])

        # Contamination
        if start_date_contamination is None:
            self.start_date_contamination = self.start_date_well
        else:
            self.start_date_contamination = start_date_contamination
            check_date_format(check_date=['start_date_contamination'])

        if compute_contamination_for_date is None:
            self.compute_contamination_for_date = self.start_date_well + timedelta(days=365.24*50)
        else:
            self.compute_contamination_for_date = compute_contamination_for_date
            check_date_format(check_date=['compute_contamination_for_date'])

        if end_date_contamination is None:
            self.end_date_contamination = end_date_contamination
        else:
            self.end_date_contamination = end_date_contamination
            check_date_format(check_date=['end_date_contamination'])

            if self.end_date_contamination < self.start_date_contamination:
                raise ValueError('Error, "end_date_contamination" is before "start_date_contamination". Please enter an new "end_date_contamination" or "start_date_contamination" ')


        ''' Check logical things here for contamination'''
        if self.compute_contamination_for_date < self.start_date_contamination:
            raise ValueError('Error, "compute_contamination_for_date" is before "start_date_contamination". Please enter an new "compute_contamination_for_date" or "start_date_contamination" ')
        if self.compute_contamination_for_date < self.start_date_well:
            raise ValueError('Error, "compute_contamination_for_date" is before "start_date_well". Please enter an new "compute_contamination_for_date" or "start_date_well" ')
        #AH_todo @MartinvdS -> if end_date_contamination < start_date_well what to do?



        if depth_point_contamination is None:
            self.depth_point_contamination = self.ground_surface
        elif depth_point_contamination > self.ground_surface:
            self.depth_point_contamination = self.ground_surface
        else:
            self.depth_point_contamination = depth_point_contamination

        # Point Contamination
        self.point_input_concentration = point_input_concentration
        self.distance_point_contamination_from_well = np.array([distance_point_contamination_from_well])
        self.discharge_point_contamination = discharge_point_contamination

        # Model size
        self.source_area_radius = source_area_radius
        self.number_of_spreading_distance = number_of_spreading_distance
        self.model_radius = model_radius
        self.model_width = model_width
        self.model_radius_computed = None #from comment MK, initialize to None
        self.relative_position_starting_points_radial = relative_position_starting_points_radial
        self.relative_position_starting_points_in_basin = relative_position_starting_points_in_basin
        self.relative_position_starting_points_outside_basin = relative_position_starting_points_outside_basin

        # Temperature
        if self.temperature_vadose_zone is None:
            self.temperature_vadose_zone = self.temp_water
        if self.temperature_shallow_aquifer is None:
            self.temperature_shallow_aquifer = self.temp_water
        if self.temperature_target_aquifer is None:
            self.temperature_target_aquifer = self.temp_water

        #Modpath params
        if particle_release_day is None:
            self.particle_release_day = particle_release_day
        else:
            self.particle_release_day = particle_release_day

        self.ncols_near_well = ncols_near_well
        self.ncols_far_well = ncols_far_well
        if nlayers_shallow_aquifer is None:
            self.nlayers_shallow_aquifer = int(self.thickness_shallow_aquifer)
        else:
            self.nlayers_shallow_aquifer = int(max(1,nlayers_shallow_aquifer))

        if nlayers_target_aquifer is None:
            self.nlayers_target_aquifer = int(self.thickness_target_aquifer)
        else:
            self.nlayers_target_aquifer = int(max(1,nlayers_target_aquifer))

        top_shallow_aquifer = self.bottom_vadose_zone_at_boundary
        top_target_aquifer = self.bottom_shallow_aquifer

        if top_filterscreen	is None:
            self.top_filterscreen = top_target_aquifer
        if bottom_filterscreen is None:
            self.bottom_filterscreen = self.bottom_target_aquifer
        if top_gravelpack is None:
            self.top_gravelpack = top_target_aquifer
        if bottom_gravelpack is None:
            self.bottom_gravelpack = self.bottom_target_aquifer
        if top_clayseal is None:
            self.top_clayseal = self.ground_surface
        if bottom_clayseal is None:
            self.bottom_clayseal = top_target_aquifer

        # DIAMETERS
        if diameter_gravelpack is None:
            self.diameter_gravelpack = self.diameter_borehole
        if inner_diameter_gravelpack is None:
            self.inner_diameter_gravelpack = diameter_filterscreen
        if inner_diameter_filterscreen is None:
            self.inner_diameter_filterscreen = self.diameter_filterscreen
        if diameter_clayseal is None:
            self.diameter_clayseal = diameter_borehole
        if inner_diameter_clayseal is None:
            self.inner_diameter_clayseal = self.diameter_filterscreen

        #other default params
        if hor_permeability_gravelpack is None:
            self.hor_permeability_gravelpack = 1000
        if hor_permeability_clayseal is None:
            self.hor_permeability_clayseal = 0.001
        if vertical_anisotropy_gravelpack is None:
            self.vertical_anisotropy_gravelpack = 1.0
        if vertical_anisotropy_clayseal is None:
            self.vertical_anisotropy_clayseal = 1.0

        if model_radius is None:
            if self.schematisation_type == 'phreatic':
                self.model_radius = math.sqrt(abs(self.well_discharge /
                    (math.pi * self.recharge_rate))) #AH SA*recharge = well_discharge
            elif self.schematisation_type == 'semiconfined':
                self.model_radius = math.sqrt(self.vertical_resistance_shallow_aquifer * self.KD) * 3 # spreading_distance*3


    def make_modflow_schematisation(self,):

        """ Returns dicitonaries of the different parameters for MODFLOW schematisation. """
        if self.schematisation_type == 'phreatic':
            # Additional meter added to the model radius for the fixed head boundary
            # only for the phreatic case, not for the semiconfined case
            # @MvdS: moet hier nog een meter bij?
            self.model_radius_computed = self.model_radius

            # only outer_boundary for phreatic
            ibound_parameters = {
                'inner_boundary_top_shallow_aquifer': {
                    'top': self.ground_surface, #self.bottom_vadose_zone_at_boundary,
                    'bot': self.bottom_vadose_zone_at_boundary,
                    'xmin': 0,
                    'xmax': self.diameter_filterscreen/2., #self.diameter_filterscreen/2.,
                    'ibound':0,
                    },
                'inner_boundary_shallow_aquifer': {
                    'top': self.bottom_vadose_zone_at_boundary,
                    'bot': self.bottom_shallow_aquifer,
                    'xmin': 0,
                    'xmax': self.diameter_filterscreen/2.,
                    'ibound':0,
                    },

                self.well_name: { #AH nov. 11 changed to from 'inner_boundary_target_aquifer' to well_name (default 'well1')
                    'head': self.bottom_shallow_aquifer,
                    'top': self.bottom_shallow_aquifer,
                    'bot': self.bottom_target_aquifer,
                    'xmin': 0,
                    'xmax': self.diameter_filterscreen/2.,
                    'ibound':-1,
                    }
                }

        elif self.schematisation_type == 'semiconfined':
            # ibound at the model radius (no additional meter added)
            self.model_radius_computed = self.model_radius

            # only top_boundary for semiconfined
            ibound_parameters = {
                'top_boundary_well': {
                    'head': self.bottom_vadose_zone_at_boundary,
                    'top': self.bottom_vadose_zone_at_boundary, # OLD: 10 cm ficticous thickness to allow head boundary
                    'bot': self.bottom_vadose_zone_at_boundary -0.1,  # NEW: ibound =-1 up to 10 cm below layer 1 top
                    'xmin': 0.,
                    'xmax': self.diameter_borehole/2.,
                    'ibound': 0
                    },
                'top_boundary1': {
                    'head': self.bottom_vadose_zone_at_boundary,
                    'top': self.bottom_vadose_zone_at_boundary, # OLD: 10 cm ficticous thickness to allow head boundary
                    'bot': self.bottom_vadose_zone_at_boundary -0.1,  # NEW: ibound =-1 up to 10 cm below layer 1 top
                    'xmin': self.diameter_borehole/2.,
                    'xmax': self.model_radius_computed,
                    'ibound': -1,
                    },
                    'inner_boundary_shallow_aquifer': {
                    'top': self.bottom_vadose_zone_at_boundary,
                    'bot':self.bottom_shallow_aquifer,
                    'xmin':0,
                    'xmax':self.diameter_filterscreen/2, #check this
                    'ibound':0,
                    }
                }

        # make dictionaries of each grouping of parameter
        simulation_parameters = {}
        for key in ['schematisation_type',
                    'computation_method',
                    'temp_correction_Koc',
                    'temp_correction_halflife',
                    'biodegradation_sorbed_phase',
                    'compute_thickness_vadose_zone',
                    'start_date_well',
                    'start_date_contamination',
                    'compute_contamination_for_date']:
            simulation_parameters[key] = getattr(self, key)



        endpoint_id = {
            self.well_name: { #AH nov. 11 changed to 'well_name'
                'top': self.top_filterscreen,
                'bot': self.bottom_filterscreen,
                'xmin': 0.0,
                'xmax': self.diameter_filterscreen/2,
                },
            }

        #mesh refinement for the Modflow grid
        mesh_refinement = {
            'mesh_refinement1': {
                'xmin': self.diameter_borehole/2,
                'xmax': self.thickness_target_aquifer,
                'ncols': self.ncols_near_well, #indicates the number of columns close to the well
                },

            'mesh_refinement2': {
                'xmin': self.thickness_target_aquifer, #@Martin from email... correct? self.diameter_gravelpack,
                'xmax': self.model_radius,
                'ncols': self.ncols_far_well
                },}

        # Aquifer parameters dcitionary
        geo_parameters  = {
            'vadose_zone': {
                'vadose': True,
                'top': self.ground_surface,
                'bot': self.bottom_vadose_zone_at_boundary,
                'xmin': self.diameter_filterscreen/2., # self.diameter_gravelpack/2,
                'xmax': self.model_radius_computed,
                'porosity': self.porosity_vadose_zone,
                'moisture_content': self.moisture_content_vadose_zone,
                'solid_density': self.solid_density_vadose_zone,
                'fraction_organic_carbon': self.fraction_organic_carbon_vadose_zone,
                'redox': self.redox_vadose_zone,
                'dissolved_organic_carbon': self.dissolved_organic_carbon_vadose_zone,
                'pH': self.pH_vadose_zone,
                'temp_water': self.temperature_vadose_zone,
                'grainsize': self.grainsize_vadose_zone,
                'hk': self.hor_permeability_shallow_aquifer,
                'vani': self.vertical_anisotropy_shallow_aquifer
                },
            'shallow_aquifer': {
                'top': self.bottom_vadose_zone_at_boundary,
                'bot': self.bottom_shallow_aquifer,
                'xmin': self.diameter_filterscreen/2., # self.diameter_gravelpack/2,
                'xmax': self.model_radius_computed,
                'porosity': self.porosity_shallow_aquifer,
                'solid_density': self.solid_density_shallow_aquifer,
                'fraction_organic_carbon': self.fraction_organic_carbon_shallow_aquifer,
                'redox': self.redox_shallow_aquifer,
                'dissolved_organic_carbon': self.dissolved_organic_carbon_shallow_aquifer,
                'pH': self.pH_shallow_aquifer,
                'temp_water': self.temperature_shallow_aquifer,
                'grainsize': self.grainsize_shallow_aquifer,
                'hk': self.hor_permeability_shallow_aquifer,
                'vani': self.vertical_anisotropy_shallow_aquifer,
                'nlayers': self.nlayers_shallow_aquifer,
                },
            'target_aquifer': {
                'top': self.bottom_shallow_aquifer,
                'bot': self.bottom_target_aquifer,
                'xmin': self.diameter_filterscreen/2., # self.diameter_gravelpack/2,
                'xmax': self.model_radius_computed,
                'porosity': self.porosity_target_aquifer,
                'solid_density': self.solid_density_target_aquifer,
                'fraction_organic_carbon': self.fraction_organic_carbon_target_aquifer,
                'redox': self.redox_target_aquifer,
                'dissolved_organic_carbon': self.dissolved_organic_carbon_target_aquifer,
                'pH': self.pH_target_aquifer,
                'temp_water': self.temperature_target_aquifer,
                'grainsize': self.grainsize_target_aquifer,
                'hk': self.hor_permeability_target_aquifer,
                'vani': self.vertical_anisotropy_target_aquifer,
                'nlayers': self.nlayers_target_aquifer, #AH maybe this will change to vertical_resolution
                },
            'gravelpack1': {
                'top': self.top_gravelpack,
                'bot': self.bottom_gravelpack,
                'xmin': self.inner_diameter_gravelpack/2.,
                'xmax': self.diameter_gravelpack/2.,
                'porosity': self.porosity_gravelpack,
                'hk': self.hor_permeability_gravelpack,
                'vani': self.vertical_anisotropy_gravelpack,
                'grainsize': self.grainsize_gravelpack,
                },
            'clayseal1':{
                'top': self.top_clayseal,
                'bot': self.bottom_clayseal,
                'xmin': self.inner_diameter_clayseal/2, #@MartinvdS correct?
                'xmax': self.diameter_clayseal/2,
                'porosity': self.porosity_clayseal,
                'hk': self.hor_permeability_clayseal,
                'vani': self.vertical_anisotropy_clayseal,
                'grainsize': self.grainsize_clayseal,
                },


            }

        if self.schematisation_type == 'semiconfined':
            well_parameters = {
                self.well_name: {
                    'well_discharge': self.well_discharge, 
                    'top': self.top_filterscreen,
                    'bot': self.bottom_filterscreen,
                    'xmin': 0.0,
                    'xmax': self.diameter_filterscreen/2,
                    },
                }
        # elif will come here for the BAR, RBF cases AH_todo
        else:
            # no well_paramters dictionary needed for the phreatic case
            well_parameters = {
            }

        if self.schematisation_type == 'semiconfined':
            recharge_parameters = {
            }
        else:
            recharge_parameters = {
                'source1': {
                    # source1 -> recharge & diffuse sources
                    'substance_name': self.substance,
                    'organism_name': self.organism,
                    'recharge': self.recharge_rate,
                    'xmin': self.diameter_borehole/2,
                    'xmax': self.model_radius_computed,
                    'dissolved_organic_carbon': self.dissolved_organic_carbon_infiltration_water,
                    'TOC': self.total_organic_carbon_infiltration_water,
                    'input_concentration': self.diffuse_input_concentration,
                    },
                # 'source2' :{}> surface water (BAR & RBF) #TODO @MartinvdS come back to this when we start this module
            }

        # concentration_boundary assigns location (and substance/organism) of diffuse contamination
        concentration_boundary_parameters = {
            'source1': { # source1 -> recharge & diffuse sources
                'substance_name': self.substance,
                'organism_name': self.organism,
                'recharge': self.recharge_rate,
                'xmin': self.diameter_borehole/2 + 0.1,
                'xmax': self.model_radius,
                'zmin': self.bottom_vadose_zone_at_boundary - 0.05,
                'zmax': self.bottom_vadose_zone_at_boundary - 0.05,
                'dissolved_organic_carbon': self.dissolved_organic_carbon_infiltration_water,
                'TOC': self.total_organic_carbon_infiltration_water,
                'input_concentration': self.diffuse_input_concentration,
                },
        # 'source2' :{}> surface water (BAR & RBF) #@MartinvdS come back to this when we start this module
        }

        # Create point diciontary if point source concentration specified,
        # otherwise pass empty dictionary
        if self.point_input_concentration is None:
            point_parameters = {}
        else:
            point_parameters = {
                'point1': {
                    'substance_name': self.substance,
                    'organism_name': self.organism,
                    'input_concentration': self.point_input_concentration,
                    'x_start': float(self.distance_point_contamination_from_well[0]), #AH Nov. 5 -> I think this might casue problems later in the loops for calculations later but deal with it when it comes
                    'z_start': self.depth_point_contamination,
                    'discharge': self.discharge_point_contamination,
                    },
                # 'point2': {} #AH_todo if there is more than one point source,
                # input a dictionary directly, @martinvdS we will discuss this
                }

        #AH eventially to be computed by QSAR"
        #MK: this does not belong to the schematization class. The same
        # schematizaiton can used with different substances. So substance information
        # and analyticalwell instance are input for the Concentration class
        # And indeed, they can be computed by the QSAR, but that is out of the scope
        # of this package, this package should be able to accept a dictionary like
        # you define below.
        # AH_todo @MartinK -> ok, what should be different then?

        ### substance_parameters replaced by removal_parameters in SubstanceTransport class (Substance_Transport.py)
        # substance_parameters = {
        #         'substance_name': self.substance,
        #         'log_Koc': self.partition_coefficient_water_organic_carbon,
        #         'pKa': self.dissociation_constant,
        #         'omp_half_life': {
        #             'suboxic': self.halflife_suboxic,
        #             'anoxic': self.halflife_anoxic,
        #             'deeply_anoxic': self.halflife_deeply_anoxic,
        #             },
        #     }

        # MK:should this be removed? #AH_todo @MartinvdS -> do we keep this dict if it is blank?
        bas_parameters = {
            }

        # Dictionaries returned as attributes of the object
        # simulation parameters (and maybe more) are now stored exactly the same way in two places (as dict and as attribute,
        # where the dict is a new attribute) this is not good: confusing and source for errors. self.simulation_parameters
        # is only used in research and test functions, so not the functionality of the code, so it should be removed
        # altogether
        # @Steven_todo: check: worden attributen tm regel 896 gebruikt voor de functionaliteit in de modelcode zel (dus niet alleen voor research/tests)?
        # Als alleen voor research/tests: dan attribuut verwijderen. Code in research / testbestanden aanpassen dien overeenkomstig.
        self.simulation_parameters = simulation_parameters
        self.endpoint_id = endpoint_id
        self.mesh_refinement = mesh_refinement
        self.geo_parameters = geo_parameters
        self.ibound_parameters = ibound_parameters
        self.recharge_parameters = recharge_parameters
        self.concentration_boundary_parameters = concentration_boundary_parameters
        self.well_parameters = well_parameters
        self.point_parameters = point_parameters
        # self.substance_parameters = substance_parameters
        self.bas_parameters = bas_parameters

    # # Here are functions to calculate the travel time through vadose zone, shared functions for
    # # Analytical and Modflow models
    # def _create_radial_distance_array(self):

    #     ''' Create array of radial distances from the well to a maximum value, radial distance recharge, which
    #     is the distance from the well needed to recharge the well to meet the pumping demand. '''
    #     # ah_todo change this to single array of 0.001 to 100
    #     # right now we have this set up to directly compare with P. Stuyfzand's results
    #     # in the 'final' version of the model can be a more fine mesh
    #     # but need to think about how we change the testing
    #     fraction_flux = np.array([0.00001, 0.0001, 0.001, 0.005])
    #     fraction_flux = np.append(fraction_flux, np.arange(0.01, 1, 0.01))
    #     self.fraction_flux = np.append(fraction_flux, [0.995, 0.9999])

    #     radial_distance = self.radial_distance_recharge * \
    #         np.sqrt(self.fraction_flux)

    #     if self.schematisation_type == 'semiconfined':

    #         radial_distance = np.append(radial_distance,
    #                                 [(radial_distance[-1] + ((self.spreading_distance * 3) - radial_distance[-1]) / 3),
    #                                 (radial_distance[-1] + 2 *
    #                                 ((self.spreading_distance * 3) - radial_distance[-1]) / 3),
    #                                     (self.spreading_distance * 3)])

    #     self.radial_distance = radial_distance

    # def _calculate_hydraulic_head_phreatic(self, distance):
    #     ''' Calcualtes the hydraulic head distribution for the phreatic schematisation case

    #     Parameters
    #     ----------
    #     distance: array
    #         Array of distance(s) [m] from the well, for diffuse sources given as
    #         the array radial_distance, for point sources given as the distance_point_contamination_from_well.

    #     Returns
    #     -------
    #     head: array
    #         Hydraulic head for each point in the given distance array, [mASL].
    #     '''
    #     ##@Alex: check if this function works properly if a negative well_discharge is given)
    #     head = (self.groundwater_level + (self.well_discharge
    #                 / (2 * math.pi * self.KD))
    #                 * np.log(self.radial_distance_recharge / distance))
    #     return head


    # #ah_todo possibly split this into making the thickness, head and travel time (3 functions)? Advantages?
    # def _calculate_unsaturated_zone_travel_time_phreatic (self,
    #                                                     distance,
    #                                                     depth_point_contamination=None):
    #     '''
    #     Calculates the travel time in the unsaturated zone for the phreatic case.

    #     Parameters
    #     ----------
    #     distance: array
    #         Array of distance(s) [m] from the well.
    #         For diffuse sources 'distance' is the 'radial_distance' array.
    #         For point sources 'distance' is 'distance_point_contamination_from_well'.
    #     depth_point_contamination: float
    #         Depth [mASL] of the point source contamination, if only diffuse contamination None is passed.
    #         MK: depth_point_contamination is an attribute right? dont need for argument here. Moreover, the
    #             if depth_point_contamination is None:
    #         case seem to be already caught in the __init__ of Schematisation class
    #         AH_todo @MartinK -> yes it is an attribute, but this needs to be passed in the function to
    #         distinguish between the diffuse and point source calculations
    #         #AH_todo review the logic of how the diffusion/point sources are calculated with @martinK

    #     Returns
    #     -------
    #     travel_time_unsaturated: array
    #         Travel time in the unsaturated zone for each point in the given distance array, [days].
    #     thickness_vadose_zone_drawdown: array
    #         Drawdown of the hydraulic head in the vadose zone for each point in the given distance array, [m].
    #     head: array
    #         Hydraulic head for each point in the given distance array, [mASL].
    #     '''

    #     head = self._calculate_hydraulic_head_phreatic(distance=distance)

    #     if depth_point_contamination is None:
    #         thickness_vadose_zone_drawdown = (self.groundwater_level
    #                                             + self.thickness_vadose_zone_at_boundary) - head

    #         travel_time_unsaturated = (((thickness_vadose_zone_drawdown
    #                                 - self.thickness_full_capillary_fringe)
    #                                 * self.moisture_content_vadose_zone
    #                                 + self.thickness_full_capillary_fringe
    #                                 * self.porosity_vadose_zone)
    #                             / self.recharge_rate)

    #     elif depth_point_contamination >= head:
    #         thickness_vadose_zone_drawdown = depth_point_contamination - head
    #         if thickness_vadose_zone_drawdown < 0:
    #             travel_time_unsaturated =  np.array([0])
    #         else:
    #             travel_time_unsaturated = (((thickness_vadose_zone_drawdown
    #                                 - self.thickness_full_capillary_fringe)
    #                                 * self.moisture_content_vadose_zone
    #                                 + self.thickness_full_capillary_fringe
    #                                 * self.porosity_vadose_zone)
    #                             / self.recharge_rate)
    #     else:
    #         travel_time_unsaturated = np.array([0])
    #         thickness_vadose_zone_drawdown = 0 #AH_todo possibly replace this with the travel distance, not thickness_vadose because this is a stand in for the travel distance?

    #     # raise warning if the thickness of the drawdown at the well is such that
    #     # it reaches the target aquifer

    #     #@MartinK -> how to raise this warning properly in the web interface? #AH_todo
    #     self.drawdown_at_well = self.ground_surface - thickness_vadose_zone_drawdown
    #     if self.drawdown_at_well[0] < self.bottom_target_aquifer:
    #         raise ValueError('The drawdown at the well is lower than the bottom of the target aquifer. Please select a different schematisation.')

    #     elif self.drawdown_at_well[0] < self.bottom_shallow_aquifer:
    #         warnings.warn('The drawdown at the well is lower than the bottom of the shallow aquifer')

    #     else:
    #         pass

    #     return travel_time_unsaturated, thickness_vadose_zone_drawdown, head


    # def _calculate_travel_time_unsaturated_zone(self,
    #                                             distance=None,
    #                                             depth_point_contamination=None):

    #     ''' Calculates the travel time in the unsaturated zone for the phreatic and semiconfined cases.
    #     The travel time is returned as an attribute of the object.

    #     MK: interesting that you put this method in this class and not in the Analytical Well class. From
    #     a data point of view it is good that you do it here because all the data you need is in this class. However,
    #     from a workflow point of view it makes less sense. Because the schematization is defined  in this class,
    #     not its derived quantities. I would move it to AnalyticalWell class. /UNLESS/ it is used by Modflow as well.
    #     # AH_todo, @MartinK -> yes, this is needed by Modflow schematisation, so leave here
    #     # @martinvdS, confirm that Modflow schematisation still needs the unsaturated zone travel times calculated in this class

    #     Parameters
    #     ----------
    #     distance: array, optional
    #         Array of distance(s) [m] from the well.
    #         For diffuse sources 'distance' is the 'radial_distance' array.
    #         For point sources 'distance' is 'distance_point_contamination_from_well'.
    #         MK: what if None? see: https://numpydoc.readthedocs.io/en/latest/format.html#parameters
    #         -> #ah_todo @MartinK - ok I assume then I just add optional to the description?
    #     depth_point_contamination: float, optional
    #         Depth [mASL] of the point source contamination, if only diffuse contamination None is passed.
    #         MK: what if None? -> #ah_todo @MartinK - ok I assume then I just add optional to the description?

    #     Returns
    #     -------
    #     travel_time_unsaturated: array
    #         Travel time in the unsaturated zone for each point in the given distance array returned as
    #         attrubute of the function, [days].

    #     '''

    #     self.spreading_distance = math.sqrt(self.vertical_resistance_shallow_aquifer * self.KD)

    #     # AH do not change to model_radius, since the radial distance for recharge is based on the phreatic value for BOTH cases
    #     self.radial_distance_recharge =  math.sqrt(abs(self.well_discharge
    #                                                 / (math.pi * self.recharge_rate )))

    #     # Diffuse source or regular travel time calculation use radial distance array
    #     # AH_todo alter here if we use multiple point sources
    #     # what does None mean? or not None?
    #     if distance is None:
    #         self._create_radial_distance_array()
    #         distance= self.radial_distance
    #     else:
    #         distance = distance

    #     '''Equation A.11 in TRANSATOMIC report '''

    #     if self.schematisation_type =='phreatic':
    #         # MK: I think a better structure would have been:
    #         # _calculate_unsaturated_zone_travel_time_phreatic as main method that calls another method with the
    #         # logic it has in common with schematization_type=semiconfined
    #         #AH I don't see the advantage of this, @martinK
    #         travel_time_unsaturated, self.thickness_vadose_zone_drawdown, self.head = self._calculate_unsaturated_zone_travel_time_phreatic(distance= distance,depth_point_contamination=depth_point_contamination)

    #     elif self.schematisation_type == 'semiconfined':
    #         # MK: couldn't you have put depth_point_contamination at self.ground_surface in the __init__? then you
    #         # don't need this if statement here.
    #         #AH_todo review the logic of how the diffusion/point sources are calculated with @martinK
    #         if depth_point_contamination is None:
    #             travel_distance = self.ground_surface - self.groundwater_level - self.thickness_full_capillary_fringe
    #         else:
    #             #if point contamination at depth, assign ground surface to depth
    #             travel_distance =  depth_point_contamination - self.groundwater_level - self.thickness_full_capillary_fringe

    #         travel_time_unsaturated =(((travel_distance)
    #                                     * self.moisture_content_vadose_zone
    #                                     + self.thickness_full_capillary_fringe
    #                                     * self.porosity_vadose_zone)
    #                                 / self.recharge_rate)

    #         if travel_distance < 0:
    #             travel_time_unsaturated = 0

    #         # travel time in semiconfined is one value, make it array by repeating the value
    #         # MK: do you have a test for this? ...
    #         # AH_todo no not yet
    #         if not isinstance(distance, float):
    #             travel_time_unsaturated = [travel_time_unsaturated] * (len(distance))

    #     self.travel_time_unsaturated = travel_time_unsaturated