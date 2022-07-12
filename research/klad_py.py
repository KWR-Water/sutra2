    def compute_omp_removal(self):
        """ 
        Calculates the concentration in the well of each flowline. Returns
        the values in 'df_flowline' and 'df_particle' as attributes of the object.

        Returns
        -------
        df_flowline: pandas.DataFrame
            Column 'flowline_id': Integer
            Column 'flowline_type': string
            Column 'flowline_discharge': Float
            Column 'particle_release_day': Float
            Column 'input_concentration': float
            Column 'endpoint_id': Integer
            Column 'well_discharge': float
            Column 'substance': string
            Column 'removal_function': string
            Column 'total_breakthrough_travel_time': float
            The breakthrough concentration in the well for the OMP taking into account retardation.
            Column 'breakthrough_concentration': float
            The breakthrough concentration in the well for the OMP taking into account sorption
            and biodegradation.

        df_particle: pandas.DataFrame
            Column 'flowline_id': int
            Column 'zone': string
            Column 'travel_time': float
            Column 'xcoord': float
            Column 'ycoord': float
            Column 'zcoord': float
            Column 'redox': float
            Column 'temp_water': float
            Column 'travel_distance': float
            Column 'porosity': float
            Column 'dissolved_organic_carbon': float
            Column 'pH': float
            Column 'fraction_organic_carbon': float
            Column 'solid_density': float
            Column 'input_concentration': float
            Column 'steady_state_concentration': float
            The steady state concentration at the well of the OMP for the flowline, [mass/L]
            Column 'omp_half_life': float
            Column 'log_Koc': float
            Column 'pKa': float
            Column 'Koc_temperature_correction': float
            The temperature corrected Koc value, only if 'temp_correction_Koc' is 'True' in the HydroChemicalSchematisation.
            Column 'omp_half_life_temperature_corrected': float
            The temperature corrected OMP half-life value, if 'temp_correction_halflife' is 'True' in the HydroChemicalSchematisation.
            Column 'retardation': float
            Column 'breakthrough_travel_time': float
            """

        # load (diffuse) input concentration    
        self.df_flowline.loc[:,'input_concentration'] = self.well.schematisation.diffuse_input_concentration
        self.df_particle.loc[:,'input_concentration'] = None
        if 'steady_state_concentration' not in self.df_particle.columns:
            # Steady state concentration [-]
            self.df_particle.loc[:,'steady_state_concentration'] = None
        # First value of input_concentration/steady_state_concentration equals diffuse input concentration of pathline
        for pid in self.df_flowline.index:
            self.df_particle.loc[pid,"input_concentration"].iloc[0] = self.df_flowline.loc[pid,'input_concentration']
            self.df_particle.loc[pid,"steady_state_concentration"].iloc[0] = self.df_flowline.loc[pid,'input_concentration']

        # self.df_flowline['input_concentration'] = self.well.schematisation.diffuse_input_concentration
        # self.df_particle['input_concentration'] = None
        # self.df_particle['steady_state_concentration'] = None

        # self.df_particle.loc[self.df_particle.zone=='surface', 'input_concentration'] = self.well.schematisation.diffuse_input_concentration
        # self.df_particle.loc[self.df_particle.zone=='surface', 'steady_state_concentration'] = self.well.schematisation.diffuse_input_concentration

        if self.well.schematisation.point_input_concentration:
            ''' point contamination '''
            # need to take into account the depth of the point contamination here....
            # need to change the df_particle and df_flowline to only be the flowlines for the point contamination flowline(s)
            # use a single point contamination for now
            # FIRST recalculate the travel times for the contamination, then initialize the class

            #only for a SINGLE point contamination
            distance = self.well.schematisation.distance_point_contamination_from_well
            depth = self.well.schematisation.depth_point_contamination
            cumulative_fraction_abstracted_water = (math.pi * self.well.schematisation.recharge_rate
                                                        * distance ** 2)/abs(self.well.schematisation.well_discharge)
            ind = self.df_particle.flowline_id.iloc[-1]

            if self.well.schematisation.schematisation_type == 'phreatic':
                head = self.well.schematisation._calculate_hydraulic_head_phreatic(distance=distance)
                self.df_flowline, self.df_particle = self.well._add_phreatic_point_sources(distance=distance,
                                            depth_point_contamination=depth,
                                            cumulative_fraction_abstracted_water=cumulative_fraction_abstracted_water)

            elif self.well.schematisation.schematisation_type == 'semiconfined':
                bottom_vadose_zone = self.well.schematisation.bottom_vadose_zone_at_boundary

                self.df_flowline, self.df_particle = self.well._add_semiconfined_point_sources(distance=distance,
                                        depth_point_contamination=depth,  )

            self.df_particle['flowline_id'] = self.df_particle['flowline_id'] + ind

            self.df_flowline['input_concentration'] = self.well.schematisation.point_input_concentration
            self.df_particle['input_concentration'] = None
            self.df_particle['steady_state_concentration'] = None
            for pid in self.df_particle.flowline_id.unique():
                self.df_particle.loc[pid, 'input_concentration'].iloc[0] = self.well.schematisation.point_input_concentration
                self.df_particle.loc[pid, 'steady_state_concentration'].iloc[0] = self.well.schematisation.point_input_concentration

            self.df_flowline['flowline_id'] = self.df_flowline['flowline_id'] + ind
            self.df_flowline['flowline_type'] = "point_source"
            self.df_flowline['flowline_discharge'] = abs(self.well.schematisation.discharge_point_contamination)

            #AH_todo, something here to loop through different point sources if more than one.

            self.df_particle = self.df_particle.append(df_particle)
            self.df_particle.reset_index(drop=True, inplace=True)

            self.df_flowline = self.df_flowline.append(df_flowline)
            self.df_flowline.reset_index(drop=True, inplace=True)

            self.df_flowline['substance'] = self.removal_parameters['substance_name']