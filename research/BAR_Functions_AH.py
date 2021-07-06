
# ------------------------------------------------------------------------------
# Basin Aquifer Recharge (BAR) Functions
# ------------------------------------------------------------------------------


def calculate_travel_time_recharge_basin(length,
                                width, 
                                depth,
                                annual_basin_discharge):
    '''
    Calculated the median travel time of infiltration water in a
    recharge basin, according to Huisman & Olsthoorn (1983)

    Parameters
    ----------
    width = width of recharge basin [m];
    length = total length of recharge basin and drainage gallery [m];
    depth = mean water-filled depth of recharge basin [m];
    annual_basin_discharge = total volume of infiltration water discharged
        into basin [m3/yr]

    x = longitudinal distance from inlet, situated on one of the ends of an
        elongate basin [m];
    calculated for median distance for x = 0.632 * length #AH why 0.632?

    Returns
    -------
    travel_time: [years]
    '''

    travel_time_recharge_basin = (length * width * depth
                                  / (annual_basin_discharge * 10 ** 6
                                     / 365.25) * np.log(1 / (1 - 0.632)))
    # AH why 0.632 = median?

    return travel_time_recharge_basin


def calculate_travel_time_upper_aquifer( 
                    horizontal_distance_basin_gallery,   
                    hydraulic_conductivity_reference_value,
                    temperature,
                    porosity_recharge_basin,
                    groundwater_level_above_saturated_zone,
                     ):
    '''
    The travel time in the upper aquifer is based on design
    characteristics of the BAR system, which is rigorously
    schematized as indicated in Fig.A3.1. The schematisation
    and analytical solutions are based on Huisman & Olsthoorn
    (1983) and Stuyfzand & van der Schans (2018), with a
    MAR system composed of parallel spreading basins and
    drainage galleries in between, with phreatic aquifer
    characteristics and typical recharge system requirements.

    The design follows from the equations given below, assuming that:
    (i) flow is predominantly horizontal (Dupuit assumption;
    also excluding unsaturated zone beneath basin), requiring
    that both the spreading basin and drainage gallery fully
    penetrates the saturated thickness of the aquifer;
    (ii) the aquifer is homogeneous and isotropic,
    (iii) rainfall and evaporation can be ignored;
    (iv) all infiltration water is recovered (QIN = QOUT).

    Parameters
    ----------
    horizontal_distance_basin_gallery = horizontal distance between basin bank and drainage gallery [m];
    hydraulic_conductivity_reference_value = m/d, at specific deg. C
    temperature = deg C
    porosity_recharge_basin = porosity [-]
    groundwater_level_above_saturated_zone = normal maximum rise of watertable above H0 [m];
    
    Outputs
    -------
    median_retention_time = minimum required or median retention time in aquifer [d];
    travel_time_distribution_upper_phreatic = travel time in the upper aquifer [d]

    Other parameters
    --------------
    hydraulic_conductivity = hydraulic conductivity of aquifer [m/d];


    Not used but calculated/input in spreadsheet AH check what to do with these?
    ---------------------------
    tIP = maximum bridging period from the dynamic reservoir,
        which allows abstraction during an intake (and spreading) interruption [d];
    VDYN,A, VDYN,B = amount of water in dynamic storage, within the aquifer
        and basin respectively, which can be used to cover a (short) intake stop [Mm3];
    vIN = entry rate [m/d]
    QOUT = pumping rate, approximately equal to QIN prior to intake stop,
        however on a daily basis [m3/d];
    SY = phreatic specific yield [-];
    H0 = saturated thickness of aquifer before spreading [m];

    '''

    # horizontal_distance_basin_gallery = np.sqrt(
    #     hydraulic_conductivity * groundwater_level_above_saturated_zone * median_retention_time / porosity_recharge_basin)
    hydraulic_conductivity = (hydraulic_conductivity_reference_value
                              * ((temperature + 42.5) / (11 + 42.5)) ** 1.5)
    # AH where do these numbers come from? 11, 42.5, 1.5?

    median_retention_time = (porosity_recharge_basin * horizontal_distance_basin_gallery
                             ** 2 / (hydraulic_conductivity * groundwater_level_above_saturated_zone))

    travel_time_array = np.arange(0.50, 4.6, 0.1)

    travel_time_distribution_upper_phreatic = median_retention_time * travel_time_array

    # EPM = exponential piston model, AH percent travel time?
    percent_flux_EPM = 1 - np.exp((travel_time_distribution_upper_phreatic[0] - travel_time_distribution_upper_phreatic) / (
        median_retention_time - travel_time_distribution_upper_phreatic[0]))

    # LMP = linear piston model, AH percent travel time?
    percent_flux_LPM = np.arange(0, 1.1, 0.1)

    # AH same as the percent flux? keep variables consistent...
    # percent_travel_time_distribution = np.arange()
    percent_flux = np.append(percent_flux_LPM[0:9], percent_flux_EPM[9:])

    aquifer_dictionary = {'travel_time_basin': median_retention_time,
                          'travel_time_upper_aquifer': travel_time_distribution_upper_phreatic,
                          'percent_flux': percent_flux,
                          'percent_flux_EPM':  percent_flux_EPM,
                          'percent_flux_LPM': percent_flux_LPM,
                          'aquifer_type': 'BAR_upper_phreatic',
                          # 'travel_time_array': travel_time_array,
                          }

    return aquifer_dictionary


def calculate_travel_time_upper_aquifer_and_basins(length,
                                                   width,
                                                   depth,
                                                   annual_basin_discharge,
                                                   horizontal_distance_basin_gallery,
                                                   hydraulic_conductivity_reference_value,
                                                   temperature,
                                                   porosity_recharge_basin,
                                                   groundwater_level_above_saturated_zone,
                                                   ):

    '''Adds the median travel time in the recharge basins to the travel 
    time calculated in the upper, phreatic aquifer
    
    Inputs
    ------
    width = width of recharge basin [m];
    length = total length of recharge basin and drainage gallery [m];
    depth = mean water-filled depth of recharge basin [m];


    Returns
    -------
    Travel time distribution in the upper aquifer and basin '''

    travel_time_recharge_basin = calculate_travel_time_recharge_basin(
        length=length,
        width=width,
        depth=depth,
        annual_basin_discharge=annual_basin_discharge)

    upper_phreatic_aquifer_dict = calculate_travel_time_upper_aquifer(
        horizontal_distance_basin_gallery=horizontal_distance_basin_gallery,
        hydraulic_conductivity_reference_value=hydraulic_conductivity_reference_value,
        temperature=temperature,
        porosity_recharge_basin=porosity_recharge_basin,
        groundwater_level_above_saturated_zone=groundwater_level_above_saturated_zone,
        )

    travel_time_upper_aquifer_and_basins = travel_time_recharge_basin + upper_phreatic_aquifer_dict['travel_time_upper_aquifer']

    return travel_time_upper_aquifer_and_basins


def calculate_travel_time_deeper_aquifer(horizontal_distance_basin_gallery,
                                         hydraulic_conductivity_reference_value,
                                         temperature,
                                         porosity_recharge_basin,
                                         groundwater_level_above_saturated_zone,
                                         vertical_resistance_aquitard,
                                         thickness_second_aquifer,
                                         porosity_second_aquifer,
                                         vertical_hydraulic_head_difference,
                                         ):
    '''
    For the deeper part we have 2 options: (i) it represents
    the second aquifer, which is semiconfined and (deeply)
    anoxic situated below an aquitard that separates it from
    the upper aquifer, or (ii) it represents the deeper, anoxic
    parts of the upper aquifer, often with lower hydraulic
    conductivity and longer transit times. In the latter
    case, we neglect the contribution from the second aquifer
    if existing at all
    
    parameters
    ---------
    vertical_resistance_aquitard   # [d], c_V
    thickness_second_aquifer: mean thickness [m]
    porosity_second_aquifer (n), [-]
    vertical_hydraulic_head_difference: mean vertical jump 
        in hydraulic head between both aquifers (Δh).
    
    Outputs
    -------
    '''

    travel_time_second_aquifer = (0.3 *  thickness_second_aquifer
        * vertical_resistance_aquitard / vertical_hydraulic_head_difference)
    
    percent_flux_deep_aquifer = 2000 / travel_time_second_aquifer

    multiplication_factor = np.sqrt(travel_time_second_aquifer)

    percent_flux = np.arange(0,110,10)

    # EPM = exponential piston model
    travel_time_distribution_EPM = percent_flux / 100

    travel_time_upper_aquifer_dict = calculate_travel_time_upper_aquifer(
        horizontal_distance_basin_gallery=horizontal_distance_basin_gallery,
        hydraulic_conductivity_reference_value=hydraulic_conductivity_reference_value,
        temperature=temperature,
        porosity_recharge_basin=porosity_recharge_basin,
        groundwater_level_above_saturated_zone=groundwater_level_above_saturated_zone,
        )

    travel_time_distribution_deeper_aquifer = (travel_time_upper_aquifer_dict['travel_time_upper_aquifer']
                                                * multiplication_factor)

    return travel_time_distribution_deeper_aquifer


def travel_time_distribution_BAR(length,
                                 width,
                                 depth,
                                 annual_basin_discharge,
                                 horizontal_distance_basin_gallery,
                                 hydraulic_conductivity_reference_value,
                                 temperature,
                                 porosity_recharge_basin,
                                 groundwater_level_above_saturated_zone,
                                 vertical_resistance_aquitard,
                                 thickness_second_aquifer,
                                 porosity_second_aquifer,
                                 vertical_hydraulic_head_difference,
                                 ):
    
    ''' Calculate travel time distribution for Basin Aquifer Recharge (BAR)
    systems

    Calculates the median travel time distribution (TTD) of the basin, 
    the TTD of the upper aquifer and TTD of the deeper aquifer.

    Assumptions: from \cite{Stuyfzand2020}
    The following assumptions are made: 
    (i) the water inlet is situated at one of the ends of an elongated 
    recharge basin 
    (ii) the infiltration rate is uniformly distributed
    (iii) rainfall and evaporation do not have a significant impact.

    Parameters
    -----------
    width = width of recharge basin [m];
    length = total length of recharge basin and drainage gallery [m];
    depth = mean water-filled depth of recharge basin [m];

    horizontal_distance_basin_gallery = horizontal distance between basin bank and drainage gallery [m];
    hydraulic_conductivity_reference_value = m/d, at specific deg. C
    temperature = deg C
    porosity_recharge_basin = porosity [-]
    groundwater_level_above_saturated_zone = normal maximum rise of watertable above H0 [m];

    '''

    median_retention_time = calculate_travel_time_recharge_basin(
        length=length,
        width=width,
        depth=depth,
        annual_basin_discharge=annual_basin_discharge,
        )

    upper_phreatic_aquifer_dict = calculate_travel_time_upper_aquifer(
        horizontal_distance_basin_gallery=horizontal_distance_basin_gallery,
        hydraulic_conductivity_reference_value=hydraulic_conductivity_reference_value,
        temperature=temperature,
        porosity_recharge_basin=porosity_recharge_basin,
        groundwater_level_above_saturated_zone=groundwater_level_above_saturated_zone,
        )

    travel_time_upper_aquifer_and_basins = calculate_travel_time_upper_aquifer_and_basins(
        length=length,
        width=width,
        depth=depth,
        annual_basin_discharge=annual_basin_discharge,
        horizontal_distance_basin_gallery=horizontal_distance_basin_gallery,
        hydraulic_conductivity_reference_value=hydraulic_conductivity_reference_value,
        temperature=temperature,
        porosity_recharge_basin=porosity_recharge_basin,
        groundwater_level_above_saturated_zone=groundwater_level_above_saturated_zone,
        )

    travel_time_deeper_aquifer = calculate_travel_time_deeper_aquifer(
        horizontal_distance_basin_gallery=horizontal_distance_basin_gallery,
        hydraulic_conductivity_reference_value=hydraulic_conductivity_reference_value,
        temperature=temperature,
        porosity_recharge_basin=porosity_recharge_basin,
        groundwater_level_above_saturated_zone=groundwater_level_above_saturated_zone,
        vertical_resistance_aquitard=vertical_resistance_aquitard,
        thickness_second_aquifer=thickness_second_aquifer,
        porosity_second_aquifer=porosity_second_aquifer,
        vertical_hydraulic_head_difference=vertical_hydraulic_head_difference,
        )

    column_names = ["travel_time_upper_aquifer", "travel_time_upper_aquifer_and_basins",
                        "travel_time_deeper_aquifer", "percent_flux", "well_discharge"]
    
    # AH check the flux calculations for the streamline... Pr*Q? 
    data = [upper_phreatic_aquifer_dict['travel_time_upper_aquifer'],
            travel_time_upper_aquifer_and_basins,
            travel_time_deeper_aquifer,
            upper_phreatic_aquifer_dict['percent_flux'],
            upper_phreatic_aquifer_dict['percent_flux'] * annual_basin_discharge]

    df = pd.DataFrame (data = np.transpose(data), columns=column_names)

    aquifer_dictionary = {'travel_time_basin': median_retention_time,
                          'travel_time_upper_aquifer': upper_phreatic_aquifer_dict['travel_time_upper_aquifer'],
                          'travel_time_upper_aquifer_and_basins': travel_time_upper_aquifer_and_basins,
                          'travel_time_deeper_aquifer': travel_time_deeper_aquifer,
                          'percent_flux': upper_phreatic_aquifer_dict['percent_flux'],
                          'df_output': df,
                          'aquifer_type': 'BAR',
                         }

    return aquifer_dictionary

