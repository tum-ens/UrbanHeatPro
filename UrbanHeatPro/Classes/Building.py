"""
Building.py
A. Molar-Cruz @ TUM ENS
"""

import os
from random import randint

import numpy as np
import scipy

import UrbanHeatPro.Functions as UrbanHeatPro
from .HotWaterDemand import HotWaterDemand
from .HotWaterDemand_D import HotWaterDemand_D
from .SpaceHeatingDemand import SpaceHeatingDemand


# import ipdb


class Building:
    # --------------------------------------------------------------------------------
    def __init__(self, b, building_stock_stats,
                 # Simulation
                 dt_vectors, resolution, number_of_typ_days, weights,
                 Tamb, I,
                 _space_heating, _hot_water, _energy_only,
                 # Space heating demand
                 Tb0_str, dTset, dT_per_hour, eta, thermal_inertia,
                 _active_population, _workday_weekend, sh_prob,
                 _solar_gains, _internal_gains,
                 _night_set_back, schedule_nsb, T_nsb, power_reduction,
                 # Hot water demand
                 Tw, dhw_prob, hw_tank_limit, hw_flow,
                 day_vector, seasonal_vector, min_vector,
                 # Results
                 result_dir, plot, save, debug):

        # General data
        # --------------------------------------------------
        # Building dataframe (from SynCity)
        self.building = b  # Building dataframe

        # Simulation
        self.dt_vector = dt_vectors[0]  # Vector of time steps as datetime objects
        self.dt_vector_excel = dt_vectors[1]  # Vector of time steps as excel date
        self.resolution = resolution  # Temporal resolution in min
        self.number_of_typ_days = number_of_typ_days  # Number of typical days
        self.weights = weights  # Weights of typical days

        # Flags
        self._space_heating = _space_heating  # Calculate space heating demand?
        self._hot_water = _hot_water  # Calculate hot water demand?
        self._energy_only = _energy_only  # Calculate only aggregated demand?

        # External factors
        self.Tamb = Tamb  # Ambient temperature vector in degC
        self.I = I  # Solar radiation vector in W/m2 [I_Gh, I_Dh, I_ex, hs]
        self.eta = eta  # Heating process efficiency
        self.thermal_inertia = thermal_inertia  # Thermal inertia of the heating system

        # Result directory
        self.result_dir = result_dir + '/Buildings/'  # Directory where results are stored

        # Space heating demand
        # --------------------------------------------------
        # Building stock statistics
        ## Residential
        self.FOOTPRINT = building_stock_stats[0][0]  # Footprint area of building typologies (TABULA)
        self.ARATIO = building_stock_stats[0][1]  # Area ratio [Roof, Wall, Window]
        self.WRATIO_ORIENTATION = building_stock_stats[0][2]  # Window ratio [East, South, West, North]
        self.FLOORS = building_stock_stats[0][3]  # Number of floors
        self.U_RES = building_stock_stats[0][4]  # U-values for residential buildings
        self.V_RES = building_stock_stats[0][5]  # Air flow rate (ventilation losses) for residential
        self.C_RES = building_stock_stats[0][6]  # Thermal cap for residential buildings
        self.CURRENT_REF_RES = building_stock_stats[0][7][0]  # Current percentage of residential refurbished buildings
        self.MAX_REF_RES = building_stock_stats[0][7][1]  # Max percentage of residential refurbished buildings
        self.SINGLE_DWELLING = building_stock_stats[0][8]  # Percentage of single dwellings for SFH and TH
        self.AVG_DWELLING_SIZE = building_stock_stats[0][9]  # Average dwelling size in m2
        self.HOUSEHOLD_SIZE = building_stock_stats[0][10]  # Household size for dwelling size categories
        self.STOCK_RES = building_stock_stats[0][17]  # Building stock statistics for residential
        ## Non-residential
        self.U_NRES = building_stock_stats[0][11]  # U-values for non-residential buildings
        self.V_NRES = building_stock_stats[0][12]  # Air flow rate (ventilation losses) for non-residential
        self.C_NRES = building_stock_stats[0][13]  # Thermal cap for non-residential buildings
        self.CURRENT_REF_NRES = building_stock_stats[0][14][0]  # Current percentage of residential refurbished
        # buildings
        self.MAX_REF_NRES = building_stock_stats[0][14][1]  # Max percentage of residential refurbished buildings
        self.STOCK_NRES = building_stock_stats[0][18]  # Building stock statistics for non-residential
        ## Both
        self.TSET = building_stock_stats[0][15]  # Target temperatures per building use
        self.SCHEDULE = building_stock_stats[0][16]  # Active hours per building use

        # Building thermal properties
        self.Tb0_str = Tb0_str  # Initial building temperature as string: 'ambient' or 'Tset'
        self.dTset = dTset  # Delta temperature (for Tset_min, Tset_max)
        self.dT_per_hour = dT_per_hour  # Maximum dT allowed in building per hour [degC]

        # Activity and occupancy in building
        self._active_population = _active_population  # Consider active population for occupancy vector
        self._workday_weekend = _workday_weekend  # Consider dif between workdays and weekends
        self.sh_prob = sh_prob  # Probability vector of using space heating

        # Heat gains
        self._solar_gains = _solar_gains
        self._internal_gains = _internal_gains

        # DSM
        self._night_set_back = _night_set_back  # Share of buildings with nsb
        self.schedule_nsb = schedule_nsb  # [start, end] of nsb in h
        self.T_nsb = T_nsb  # Night set-back temperature in degC
        self.power_reduction = power_reduction  # Percentage of power reduced (as decimal)

        # Hot water demand
        # --------------------------------------------------
        # Domestic hot water consumption
        self.Tw = Tw  # Hot water temperature in [degC]
        self.dhw_prob = dhw_prob  # Probabilities for dhw-loads
        self.hw_tank_limit = hw_tank_limit  # Hot water tank limit as perc (decimal)
        self.hw_flow = hw_flow  # Flow to refill hot water tank in L/min

        # Seasonality
        self.day_vector = day_vector  # Vector of days in simulation time frame
        self.seasonal_vector = seasonal_vector  # Sinusoidal function for seasonal variations of DHW consumption
        self.min_vector = min_vector  # Vector of simulation time steps in minutes

        # Reporting
        # --------------------------------------------------
        self.plot = plot
        self.save = save
        self.debug = debug

    # Heating energy demand
    # --------------------------------------------------------------------------------
    def calculate_space_heating_demand(self):
        """
        Calculates building space heating demand as timeseries [W] and aggregated value [Wh]
        """

        # Number of time steps
        self.nts = len(self.dt_vector)

        # Activity and Occupancy vector
        if self._active_population:
            self.calculate_building_activity_occupancy_vector()

        # Calculate window areas
        if self._solar_gains:
            self.calculate_building_window_areas()
        else:
            self.window_areas = None

        # Initial building temperature (Tb[0])
        if self.Tb0_str == 'Tamb':
            self.Tb0 = self.Tamb[0]
        elif self.Tb0_str == 'Tset':
            self.Tb0 = self.building.Tset

        # Space heating demand
        # --------------------------------------------------
        self.my_space_heating_demand = SpaceHeatingDemand(self.dt_vector, self.resolution, self.building.heated_area,
                                                          self.Tamb, self.I, self.Tb0, self.dT_per_hour, self.eta,
                                                          self.thermal_inertia,
                                                          self.building.U, self.building.V, self.building.C,
                                                          self.building.Tset, self.dTset,
                                                          self.activity_vector, self.occupancy_vector, self.sh_prob,
                                                          self._solar_gains, self._internal_gains,
                                                          self._night_set_back, self.schedule_nsb, self.T_nsb,
                                                          self.power_reduction,
                                                          self.window_areas, [self.building.lat, self.building.lon],
                                                          self.debug)
        self.my_space_heating_demand.calculate()

        # Power [W]
        self.space_heating_power = self.my_space_heating_demand.sh_power
        self.solar_gains = self.my_space_heating_demand.solar_gains
        self.internal_gains = self.my_space_heating_demand.internal_gains

        # Building temperature [degC]
        self.Tb = self.my_space_heating_demand.Tb

        # Energy [Wh] -> ANNUAL
        if self.number_of_typ_days < 365:
            self.space_heating_energy = self.calculate_annual_demand(self.space_heating_power)
            self.solar_gains_energy = self.calculate_annual_demand(self.solar_gains)
            self.internal_gains_energy = self.calculate_annual_demand(self.internal_gains)

        else:
            self.space_heating_energy = (self.space_heating_power.sum() * self.resolution * 1 / 60)
            self.solar_gains_energy = (self.solar_gains.sum() * self.resolution * 1 / 60)
            self.internal_gains_energy = (self.internal_gains.sum() * self.resolution * 1 / 60)

        # Energy per (conditioned) area [kWh/m2]
        self.space_heating_energy_per_area = (self.space_heating_energy / self.building.heated_area)
        self.solar_gains_per_area = (self.solar_gains_energy / self.building.heated_area)
        self.internal_gains_per_area = (self.internal_gains_energy / self.building.heated_area)

        result = [self.space_heating_energy, self.solar_gains_energy, self.internal_gains_energy,
                  self.space_heating_energy_per_area, self.solar_gains_per_area, self.internal_gains_per_area]

        result = np.resize(result, (1, len(result)))[0]

        return result

    #
    def calculate_hot_water_demand(self, save_debug=False):
        """
        Calculates building hot water demand as timeseries [W] and [m3] and aggregated value [Wh]
        Only for residential buildings

        args:
            save_debug	<boolean>	Is debug file saved?

        returns:
            self.hot_water_m3
            self.hot_water_power
            self.dhw_energy
        """

        # Number of time steps
        self.nts = len(self.dt_vector)

        if self.building.use == 3:  # residential

            if not self._space_heating:
                # self.activity_vector, self.occupancy_vector
                self.calculate_building_activity_occupancy_vector()

            # Hot water demand
            # --------------------------------------------------

            # Energy only: simplified model
            if self._energy_only:
                self.my_hot_water_demand = HotWaterDemand_D(self.resolution,
                                                            self.day_vector,
                                                            self.Tw, self.building.hw_demand_day)
                self.my_hot_water_demand.calculate()

                self.dhw_m3 = self.my_hot_water_demand.dhw_m3
                self.dhw_energy = self.my_hot_water_demand.dhw_energy

                #
                # scale up to year
                if self.number_of_typ_days < 365:
                    self.dhw_m3 = self.dhw_m3 / self.number_of_typ_days * 365.
                    self.dhw_energy = self.dhw_energy / self.number_of_typ_days * 365.

            # Timeseries: probabilistic model
            else:

                # define size and state of hot water tank
                self.set_hot_water_tank_initial_state()

                # Hot water demand
                # --------------------------------------------------
                self.my_hot_water_demand = HotWaterDemand(self.dt_vector, self.resolution,
                                                          self.day_vector, self.seasonal_vector, self.activity_vector,
                                                          self.Tw, self.building.hw_demand_day, self.dhw_prob,
                                                          self.building.hw_tank_cap, self.hw_tank_limit,
                                                          self.hw_tank_volume_t0, self.hw_flow,
                                                          self.result_dir, self.building.use,
                                                          int(self.building.year_class), int(self.building.size_class),
                                                          self.building.bid,
                                                          self.debug, save_debug)
                self.my_hot_water_demand.calculate()

                if save_debug:
                    self.dhw_debug = self.my_hot_water_demand.dhw_debug

                # Power [W]
                self.hot_water_power = self.my_hot_water_demand.dhw_power

                # Flow rate [m3]
                self.hot_water_m3 = self.my_hot_water_demand.dhw_m3
                self.hot_water_tank_m3 = self.my_hot_water_demand.dhw_tank_m3

                # Energy [Wh]
                self.dhw_m3 = self.hot_water_m3.sum()
                self.dhw_energy = (self.hot_water_power.sum() * self.resolution * 1 / 60)

                #
                # scale up to year
                if self.number_of_typ_days < 365:
                    self.dhw_energy = self.calculate_annual_demand(self.hot_water_power)
                    self.dhw_m3 = self.calculate_annual_demand(self.hot_water_m3)

            result = np.array([self.dhw_energy, self.dhw_m3])


        else:  # Non-residential
            # No hot water demand
            self.hot_water_power = np.zeros([self.nts])
            self.hot_water_m3 = np.zeros([self.nts])
            self.hot_water_tank_m3 = np.zeros([self.nts])
            self.dhw_energy = 0.
            self.dhw_m3 = 0.

        return np.array([self.dhw_energy, self.dhw_m3])

    #
    def calculate_total_heat_demand(self):
        """
        Agrregate the space heating and/or hot water demand time series.
        The total time series is delayed depending on the distance to the heat plant.
        """

        # initialize matrices
        self.total_power = np.zeros([self.nts])  # Total heat demand in W
        self.total_energy = 0.  # Aggregated total heat demand in Wh

        # Space heating demand
        if self._space_heating:
            if not self._energy_only:
                self.total_power += self.space_heating_power
            self.total_energy += self.space_heating_energy

        # Hot water demand
        if self._hot_water:
            if not self._energy_only:
                self.total_power += self.hot_water_power
            self.total_energy += self.dhw_energy

        # Delayed time series
        if not self._energy_only:
            if self.building.dist_to_heat_source > 0.:
                self.calculate_delayed_timeseries()

        return np.array([self.total_energy, self.total_energy / self.building.heated_area])

    #
    def calculate_delayed_timeseries(self, flow_vel=1.):
        """
        Delays vector of heat demand depending on the distance of the building
        centroid to the (geothermal) heat plant. A flow velocity of 1 m/s in the district heating network is
        considered.
        """

        self.delayed_min_vector = np.zeros([self.nts])  # Vector of delayed time steps in min
        self.total_power_delayed = np.zeros([self.nts])  # Delayed total heat demand in W

        # calculate flow velocity in minutes
        flow_vel_min = flow_vel * 60.  # m/min

        # calculate delay
        # Delay is two times the distance between the building and the power plant (round trip: flow
        # from the power plant to the building and back)
        self.delay = 2. * self.building.dist_to_heat_source / flow_vel_min  # min

        # add delay to datetime minute vector
        for iii, dt in enumerate(self.min_vector):
            self.delayed_min_vector[iii] = self.min_vector[iii] - self.delay

        # interpolate to calculate delayed time series in base datetime vector
        f = scipy.interpolate.interp1d(self.delayed_min_vector, self.total_power, fill_value="extrapolate")
        self.total_power_delayed = f(self.min_vector)

        # calculate delayed energy demand in Wh
        self.total_energy_delayed = (self.total_power_delayed.sum() * self.resolution * 1 / 60)

    # Building parametrization
    # --------------------------------------------------------------------------------
    def parametrize_building(self):
        """
        Calculates missing building properties necessary for the heat demand calculation.
        """

        # MIN PARAMETERS
        self.bid = self.building.bid
        self.footprint_area = self.building.footprint_area
        self.use = int(self.building.use)
        self.free_walls = self.building.free_walls
        self.lat = self.building.lat
        self.lon = self.building.lon

        # RESIDENTIAL
        # -----------
        if self.building.use == 3:

            #
            # Construction year class and size class
            try:
                self.year_class = int(self.building.year_class)
                self.size_class = int(self.building.size_class)
            except:
                # categorize building according TABULA typology -> self.year_class, self.size_class
                self.categorize_building_residential()

            #
            # Number of floors
            try:
                self.floors = self.building.floors
            except:
                self.floors = False

            #
            # Area correction factor
            try:
                self.area_corr_factor = self.building.area_corr_factor
            except:
                self.area_corr_factor = False

            #
            # Heated area
            try:
                self.heated_area = self.building.heated_area
            except:
                self.heated_area = False

            # calculate floors and reference areas -> self.area_corr_factor, self.floors, self.storey_area,
            # self.heated_area
            self.calculate_areas_residential()

            #
            # Number of dwellings
            try:
                self.dwellings = self.building.dwellings
            except:
                self.dwellings = False
            # calculate number of dwellings per building -> self.dwellings, self.dwelling_size
            self.calculate_number_of_dwellings()

            #
            # Number of occupants
            try:
                self.occupants_in_building = self.building.occupants
            except:
                # calculate number of occupants per building -> self.occupants
                self.calculate_number_of_occupants_residential()

            #
            # Refurbishment level per element
            try:
                self.ref_level_roof = self.building.ref_level_roof
                self.ref_level_wall = self.building.ref_level_wall
                self.ref_level_floor = self.building.ref_level_floor
                self.ref_level_window = self.building.ref_level_window
                self.ref_level = [self.ref_level_roof, self.ref_level_wall, self.ref_level_floor, self.ref_level_window]
            except:
                # compute refurbishment level -> self.ref_level
                self.compute_current_refurbishment_level_residential()

            #
            # Envelope areas
            self.calculate_building_envelope_areas_residential()

            #
            # Thermal properties per unit area
            self.get_building_thermal_properties_per_unit_area_residential()

            #
            # Domestic hot water demand and tank capacity -> self.hw_demand_day, self.hw_tank_cap
            self.calculate_daily_hot_water_demand()
            self.parametrize_hot_water_tank()


        # NON-RESIDENTIAL
        # ---------------
        elif self.building.use == 0 or self.building.use == 1 or self.building.use == 2:

            #
            # Construction year class and size class
            try:
                self.year_class = int(self.building.year_class)
                self.size_class = int(self.building.size_class)
            except:
                # categorize building -> self.year_class, self.size_class
                self.categorize_building_non_residential()

            #
            # Number of floors
            try:
                self.floors = self.building.floors
            except:
                self.floors = False
            # calculate floors and reference areas -> self.area_corr_factor, self.floors, self.storey_area, self.heated_area
            self.calculate_areas_non_residential()

            #
            # Number of occupants
            try:
                self.occupants_in_building = self.building.occupants
            except:
                # calculate number of occupants per building -> self.occupants
                self.calculate_number_of_occupants_non_residential()
            self.dwellings = -1
            self.dwelling_size = -1.

            #
            # Refurbishment level per element
            try:
                self.ref_level_roof = self.building.ref_level_roof
                self.ref_level_wall = self.building.ref_level_wall
                self.ref_level_floor = self.building.ref_level_floor
                self.ref_level_window = self.building.ref_level_window
                self.ref_level = [self.ref_level_roof, self.ref_level_wall, self.ref_level_floor, self.ref_level_window]
            except:
                # compute refurbishment level -> self.ref_level
                self.compute_current_refurbishment_level_non_residential()

            #
            # Envelope areas
            self.calculate_building_envelope_areas_non_residential()

            #
            # Thermal properties per unit area
            self.get_building_thermal_properties_per_unit_area_non_residential()

            #
            # Hot water demand
            self.hw_demand_day = -1.
            self.hw_tank_cap = -1.

        else:
            raise ValueError('Building use does not exist')

        # RES + NRES
        # ----------

        #
        # Distance to heat source
        try:
            self.dist_to_heat_source = self.building.dist_to_heat_source
        except:
            self.dist_to_heat_source = -1.

        #
        # Equivalent thermal properties -> self.U, self.V, self.C
        self.calculate_building_thermal_properties()

        #
        # Set temperature -> self.Tset
        try:
            self.Tset = self.building.Tset
        except:
            self.calculate_building_Tset()

        #
        # Active hours -> self.active_hours
        self.calculate_building_active_hours()

        # ==============================================================
        # add parameters to building dataframe
        self.building.is_copy = False
        self.building.loc['dist_to_heat_source'] = self.dist_to_heat_source
        self.building.loc['year_class'] = self.year_class
        self.building.loc['size_class'] = self.size_class
        self.building.loc['floors'] = self.floors
        self.building.loc['storey_area'] = self.storey_area
        self.building.loc['area_corr_factor'] = self.area_corr_factor
        self.building.loc['heated_area'] = self.heated_area
        self.building.loc['window_area'] = self.env_areas[3]
        self.building.loc['dwellings'] = self.dwellings
        self.building.loc['dwelling_size'] = self.dwelling_size
        self.building.loc['occupants'] = self.occupants_in_building
        self.building.loc['ref_level_roof'] = self.ref_level[0]
        self.building.loc['ref_level_wall'] = self.ref_level[1]
        self.building.loc['ref_level_floor'] = self.ref_level[2]
        self.building.loc['ref_level_window'] = self.ref_level[3]
        self.building.loc['U'] = self.U
        self.building.loc['V'] = self.V
        self.building.loc['C'] = self.C
        self.building.loc['Tau'] = self.Tau
        self.building.loc['Tset'] = self.Tset
        self.building.loc['act_start'] = self.active_hours[0]
        self.building.loc['act_end'] = self.active_hours[1]
        self.building.loc['hw_demand_day'] = self.hw_demand_day
        self.building.loc['hw_tank_cap'] = self.hw_tank_cap

        return self.building

    #
    def update_building_refurbishment_level(self, ref_matrix_res, ref_matrix_nres):
        """
        Updates building refurbishment level based on desired refurbishment level
        from scenario and max percentage of refurbished buildings in region
        (MAX_REF_RES and MAX_REF_NRES).

        Refurbishment levels according to TABULA typology:
            1	National minimum requirement
            2	Improved standard
            3	Ambitious standard
        """

        # Parameters
        self.use = int(self.building.use)
        self.size_class = int(self.building.size_class)
        self.year_class = int(self.building.year_class)
        self.free_walls = self.building.free_walls
        self.floors = self.building.floors
        self.heated_area = self.building.heated_area

        # RESIDENTIAL
        # -----------
        if self.use == 3:
            flag = self.compute_scenario_refurbishment_level_residential(ref_matrix_res)

            if flag:
                #
                # Envelope areas
                self.calculate_building_envelope_areas_residential()

                #
                # Thermal properties per unit area
                self.get_building_thermal_properties_per_unit_area_residential()


        # NON-RESIDENTIAL
        # ---------------
        elif self.use == 0 or self.use == 1 or self.use == 2:
            flag = self.compute_scenario_refurbishment_level_non_residential(ref_matrix_nres)

            if flag:
                #
                # Envelope areas
                self.calculate_building_envelope_areas_non_residential()

                #
                # Thermal properties per unit area
                self.get_building_thermal_properties_per_unit_area_non_residential()


        else:
            raise ValueError('Building use does not exist')

        if flag:
            # Equivalent thermal properties -> self.U, self.V, self.C
            self.calculate_building_thermal_properties()

            # update values in building df
            self.building.loc['ref_level_roof'] = self.ref_level[0]
            self.building.loc['ref_level_wall'] = self.ref_level[1]
            self.building.loc['ref_level_floor'] = self.ref_level[2]
            self.building.loc['ref_level_window'] = self.ref_level[3]
            self.building.loc['U'] = self.U
            self.building.loc['V'] = self.V
            self.building.loc['C'] = self.C
            self.building.loc['Tau'] = self.Tau

        return self.building

    ## Building parameters: Space heating demand
    #
    def categorize_building_residential(self):
        """
        Probabilistic categorization building according to TABULA typologies.
        Construction year class and building type are calculated by comparing the residential
        building gross floor area (footprint_area) with the FOOTPRINT of typical buildings (from TABULA).
        Values are adapted to fit building stock statistics.
        """

        # get size of FOOTPRINT matrix
        rows = len(self.FOOTPRINT)  # construction year class
        cols = len(self.FOOTPRINT[0])  # building type

        # initialize distance vectors
        distance = np.zeros(rows * cols)
        distance_inv = np.zeros(rows * cols)  # 1/d
        row_col = [[] for _ in range(rows * cols)]
        t_distance_inv = 0  # total inv distance, sum(1/d)
        kkk = 0  # typology counter

        #
        for iii in range(0, rows):
            for jjj in range(0, cols):
                # calculate distance between building footprint_area and FOOTPRINT (TABULA) matrix
                distance[kkk] = abs(
                    self.FOOTPRINT[iii][jjj] - self.footprint_area) + 0.001  # "+ 0.001" to avoid having x/0
                # calculate inverse and weight by building stock statistics
                distance_inv[kkk] = (1.0 / distance[kkk]) * self.STOCK_RES[iii, jjj]
                row_col[kkk] = [iii, jjj]
                t_distance_inv += distance_inv[kkk]
                kkk += 1

        # normalize inverse distance vector
        norm_distance = distance_inv / t_distance_inv

        # calculate cumulative density function
        cdf = np.cumsum(norm_distance)

        # Probabilistic categorization
        # test categorization X times or until the typology exists
        value = 0
        counter = 0
        while value == 0 and counter < 5:
            rnd = np.random.uniform(0, 1, 1)
            index = np.argmax(rnd < cdf)
            ### Construction year class as int
            self.year_class = UrbanHeatPro.year_class_to_tuple(self.use, row_col[index][0])[0]
            ### Building size class as int
            if self.free_walls == 2:
                self.size_class = 1
            else:
                self.size_class = UrbanHeatPro.size_class_to_tuple(self.use, row_col[index][1])[0]
            ### test if typology exists
            value = self.FOOTPRINT[self.year_class][self.size_class]
            counter += 1

    #
    def categorize_building_non_residential(self):
        """
        Probabilistic categorization of non-residential buildings according to the
        building statistics and the following construction year classes:

        int	construction year class
         0	 < 1918
         1	 1919 - 1976
         2	 1977 - 1983
         3	 1984 - 1994
         4	 > 1995
         
         >> Source for building stock missing
        """

        # calculate cumulative density function
        cdf = np.cumsum(self.STOCK_NRES) / np.sum(self.STOCK_NRES)

        # Probabilistic categorization
        rnd = np.random.uniform(0, 1, 1)
        index = np.argmax(rnd < cdf)
        self.year_class = UrbanHeatPro.year_class_to_tuple(self.use, index)[0]

        # set size class to Nan
        self.size_class = -1

    #
    def compute_current_refurbishment_level_residential(self):
        """
        Computes the refurbishment level for the different building elements
        [roof, wall, floor, window] according to the current refurbishment statistics.
        Refurbishment levels according to TABULA typology:
            1	National minimum requirement
            2	Improved standard
            3	Ambitious standard
        """

        # check current share of refurbished buildings for building size and year class
        perc_refurbished = self.CURRENT_REF_RES[self.size_class][self.year_class]

        # assign refurbishment level for every building element
        # if random number is smaller than percentage, then refurbishment level of
        # element is 2, otherwise is 1. No building is considered as level 3.
        self.ref_level = np.ones(4, dtype=int)
        for element in range(4):  # [roof, wall, floor, window]
            rand_num = np.random.uniform(0, 1, 1)
            if rand_num <= perc_refurbished[element]:
                self.ref_level[element] = 2

    #
    def compute_current_refurbishment_level_non_residential(self):
        """
        Computes the refurbishment level for the different building components
        [roof, wall, floor, window] according to the refurbishment statistics.
        Refurbishment levels according to TABULA typology:
            1	National minimum requirement
            2	Improved standard
            3	Ambitious standard

        >> Statistics on refurbishment in non-residential buildings missing
        """

        # check current share of refurbished buildings for building size and year class
        perc_refurbished = self.CURRENT_REF_NRES[self.year_class]

        # assign refurbishment level for every building element
        # if random number is smaller than percentage, then refurbishment level of
        # element is 2, otherwise is 1. No building is considered as level 3.
        self.ref_level = np.ones(4, dtype=int)
        for element in range(4):  # [roof, wall, floor, window]
            rand_num = np.random.uniform(0, 1, 1)
            if rand_num <= perc_refurbished[element]:
                self.ref_level[element] = 2

    #
    def compute_scenario_refurbishment_level_residential(self, ref_matrix_res):
        """
        Computes the refurbishment level for the different building elements
        [roof, wall, floor, window] according to the scenario refurbishment level
        per typology and the maximum share of refurbished buildings (MAX_REF_RES).

        Refurbishment levels according to TABULA typology:
            1	National minimum requirement
            2	Improved standard
            3	Ambitious standard
        """

        # check max share of refurbished buildings per building size class and year class
        perc_max_refurbished = self.MAX_REF_RES[int(self.building.size_class)][int(self.building.year_class)]

        # check current refurbishment level
        self.ref_level = [self.building.ref_level_roof,
                          self.building.ref_level_wall,
                          self.building.ref_level_floor,
                          self.building.ref_level_window]

        # check (desired) set refurbishment level
        ref_level_set = ref_matrix_res[int(self.building.size_class)][int(self.building.year_class)]

        flag = False
        # if current level is equal or greater than desired level, then no changes
        # otherwise, the refurbishment level is likely to be improved
        for element in range(4):  # [roof, wall, floor, window]
            if self.ref_level[element] < ref_level_set[element]:

                # improve element refurbishment to desired level if random number
                # is smaller than max refurbishment percentage
                rand_num = np.random.uniform(0, 1, 1)
                if rand_num <= perc_max_refurbished[element]:
                    self.ref_level[element] = ref_level_set[element]
                    flag = True

        return flag

    #
    def compute_scenario_refurbishment_level_non_residential(self, ref_matrix_nres):
        """
        Computes the refurbishment level for the different building elements
        [roof, wall, floor, window] according to the scenario refurbishment level
        per typology and the maximum share of refurbished buildings (MAX_REF_NRES).

        Refurbishment levels according to TABULA typology:
            1	National minimum requirement
            2	Improved standard
            3	Ambitious standard
        """

        # check max share of refurbished buildings per year class
        perc_max_refurbished = self.MAX_REF_NRES[int(self.building.year_class)]

        # check current refurbishment level
        self.ref_level = [self.building.ref_level_roof,
                          self.building.ref_level_wall,
                          self.building.ref_level_floor,
                          self.building.ref_level_window]

        # check (desired) set refurbishment level
        ref_level_set = ref_matrix_nres[int(self.building.year_class)]

        flag = False
        # if current level is equal or greater than desired level, then no changes
        # otherwise, the refurbishment level is likely to be improved
        for element in range(4):  # [roof, wall, floor, window]

            if self.ref_level[element] < ref_level_set[element]:
                # improve element refurbishment to desired level if random number
                # is smaller than max refurbishment percentage
                rand_num = np.random.uniform(0, 1, 1)

                if rand_num <= perc_max_refurbished[element]:
                    self.ref_level[element] = ref_level_set[element]

                    # update thermal properties flag
                    flag = True

        return flag

    #
    def calculate_areas_residential(self):
        """
        Calculate storey area and heated/conditioned area based on definitions
        from VDI 3807.
        """

        # Area correction factor (VDI 3807-2 Section 7.1.1 Table 1)
        if not self.area_corr_factor:
            if not self.heated_area:
                self.area_corr_factor = 0.84
            else:
                self.area_corr_factor = self.heated_area / self.footprint_area

        # Floors
        if not self.floors:
            self.calculate_number_of_floors_residential()

        # Storey area (VDI 3807-1 Section 5.4.2)
        self.storey_area = self.footprint_area * self.floors

        # Heated area (VDI 3807-1 Section 5.4.2) or reference area in TABULA
        if not self.heated_area:
            self.heated_area = self.footprint_area * self.area_corr_factor * self.floors

    #
    def calculate_areas_non_residential(self):
        """
        Calculate storey area and heated/conditioned area based on definitions
        from VDI 3807.
        """

        # Area correction factor (VDI 3807-2 Section 7.1.1 Table 1)
        self.area_corr_factor = np.random.randint(75, 85) / 100.

        # Floors
        if not self.floors:
            self.calculate_number_of_floors_non_residential()

        # Storey area (VDI 3807-1 Section 5.4.2)
        self.storey_area = self.footprint_area * self.floors

        # Heated area (VDI 3807-1 Section 5.4.2)
        self.heated_area = self.footprint_area * self.area_corr_factor * self.floors

    #
    def calculate_number_of_floors_residential(self):
        """
        Calculates number of floors based on the TABULA typology.
        The number of floors calculated from TABULA are referenced to the conditioned
        or heated area but the number of floors are calculated using the storey area.
        """

        # get number of floors
        self.floors = int(np.ceil(self.FLOORS[self.year_class][self.size_class] / (1 + (1 - self.area_corr_factor))))

        # set minimum to 1
        if self.floors == 0:
            self.floors = 1

    #
    def calculate_number_of_floors_non_residential(self, left=1, mode=2, right=3):
        """
        Calculates number of floors as random sample number from the triangular
        distribution with lower limit left, peak at mode and upper limit right.

        >> Non-residential buildings are assumed to have two floors as mode and
            a maximum of three floors
            Source missing
        """

        # calculate number of floors
        self.floors = int(np.random.triangular(left, mode, right))

    #
    def calculate_number_of_dwellings(self):
        """
        Calculates number of dwellings based on the building living area and mean dwelling size.
        It is assumed that SFH and TH have only 1 or 2 dwellings which is determined using the
        single-dwelling buildings statistics. For MFH and AB, the number of dweelings is calculated
        based on the average dwelling size.
        """

        if not self.dwellings:
            # For SFH and TH
            # compare random number to single-dwelling buildings statistics
            if self.size_class == 0 or self.size_class == 1:  # SFH or TH
                perc_single_dwelling = self.SINGLE_DWELLING[self.year_class][self.size_class]
                rand_num = np.random.uniform(0, 1, 1)
                if rand_num <= perc_single_dwelling:
                    self.dwellings = 1
                else:
                    self.dwellings = 2

            # For MFH and AB
            # divide total living area by dwelling size
            else:
                # calculate random dwelling size from normal dist with mean AVG_DWELLING_SIZE
                # and sigma AVG_DWELLING_SIZE/10
                self.dwelling_size = 0
                while self.dwelling_size <= 25:  # min dwelling size [m2]
                    self.dwelling_size = np.random.normal(self.AVG_DWELLING_SIZE, self.AVG_DWELLING_SIZE / 10, 1)
                # calculate number of dwellings
                self.dwellings = int(np.ceil(self.heated_area / self.dwelling_size))

        # Dwelling size
        self.dwelling_size = self.heated_area / self.dwellings

    #
    def determine_dwelling_size_category(self):
        """
        Determine dwelling size category based on statistics
        https://ergebnisse.zensus2011.de/#StaticContent:091840148148,GWZ_4_3_2,m,table
        """

        # determine size category
        if self.dwelling_size > 0 and self.dwelling_size <= 40:
            dwelling_size_cat = 0
        elif self.dwelling_size > 40 and self.dwelling_size <= 60:
            dwelling_size_cat = 1
        elif self.dwelling_size > 60 and self.dwelling_size <= 80:
            dwelling_size_cat = 2
        elif self.dwelling_size > 80 and self.dwelling_size <= 100:
            dwelling_size_cat = 3
        elif self.dwelling_size > 100 and self.dwelling_size <= 120:
            dwelling_size_cat = 4
        elif self.dwelling_size > 120 and self.dwelling_size <= 140:
            dwelling_size_cat = 5
        elif self.dwelling_size > 140 and self.dwelling_size <= 160:
            dwelling_size_cat = 6
        elif self.dwelling_size > 160 and self.dwelling_size <= 180:
            dwelling_size_cat = 7
        elif self.dwelling_size > 180 and self.dwelling_size <= 200:
            dwelling_size_cat = 8
        else:
            dwelling_size_cat = 9

        return dwelling_size_cat

    #
    def calculate_number_of_occupants_residential(self):
        """
        Calculates number of occupants based on household size and number of dwellings
        statistics.
        """

        # calculate dwelling size category
        dwelling_size_cat = self.determine_dwelling_size_category()

        # CFD for household size
        ### x-values
        household_size = range(1, 7)
        ### CFD
        cfd = np.cumsum(self.HOUSEHOLD_SIZE[dwelling_size_cat])

        # calculate number of occupants
        self.occupants_in_building = 0
        for dwelling in range(self.dwellings):
            # get random_number and obtain x-value from cfd
            rand_num = np.random.uniform(0, 1, 1)
            index = np.argmax(rand_num <= cfd)
            # get household size (number of occupants)
            occupants_in_dwelling = household_size[index]

            # add occupants to building
            self.occupants_in_building += occupants_in_dwelling

    #
    def calculate_number_of_occupants_non_residential(self, capacity=0.1):
        """
        Calculates random number of occupants in the building based on the
        recommended area per person for different building types from
        https://www.engineeringtoolbox.com/number-persons-buildings-d_118.html.
        """

        if self.use == 0:  # commercial
            # retail, supermarket, department stores
            area_per_person = np.random.randint(3, 10)  # in m2

        elif self.use == 1:  # industrial
            # light manufacturing, heavy manufacturing
            area_per_person = np.random.randint(10, 30)  # in m2

        elif self.use == 2:  # public
            # municipal buildings, library, museum
            area_per_person = np.random.randint(3, 10)  # in m2

        # calculate number of occupants
        # buildings are assumed to be occupied at a perc of total capacity
        self.occupants_in_building = int(self.heated_area / area_per_person * capacity)

        return self.occupants_in_building

    #
    def get_building_thermal_properties_per_unit_area_residential(self):
        """
        Gets building thermal properties from TABULA Web Tool data based on the
        building typology [year_class, btype]

        returns:
            u	<list>	[u_roof, u_wall, u_floor, u_window] in W/(K m2)
            v	<list>  [v_usage, v_infiltration] in 1/h
            c	<list>	[c_roof, c_wall, c_floor] in J/(K m2)
        """

        # Transmission losses (U-values) [W/(K m2)]

        ### ROOF
        if self.ref_level[0] == 1:
            u_roof = self.U_RES[0][0][self.year_class][self.size_class]
        elif self.ref_level[0] == 2:
            u_roof = self.U_RES[0][1][self.year_class][self.size_class]
        elif self.ref_level[0] == 3:
            u_roof = self.U_RES[0][2][self.year_class][self.size_class]

        ### WALL
        if self.ref_level[1] == 1:
            u_wall = self.U_RES[1][0][self.year_class][self.size_class]
        elif self.ref_level[1] == 2:
            u_wall = self.U_RES[1][1][self.year_class][self.size_class]
        elif self.ref_level[1] == 3:
            u_wall = self.U_RES[1][2][self.year_class][self.size_class]

        ### FLOOR
        if self.ref_level[2] == 1:
            u_floor = self.U_RES[2][0][self.year_class][self.size_class]
        elif self.ref_level[2] == 2:
            u_floor = self.U_RES[2][1][self.year_class][self.size_class]
        elif self.ref_level[2] == 3:
            u_floor = self.U_RES[2][2][self.year_class][self.size_class]

        ### WINDOW
        if self.ref_level[3] == 1:
            u_window = self.U_RES[3][0][self.year_class][self.size_class]
        elif self.ref_level[3] == 2:
            u_window = self.U_RES[3][1][self.year_class][self.size_class]
        elif self.ref_level[3] == 3:
            u_window = self.U_RES[3][2][self.year_class][self.size_class]

        # Ventilation losses (air exchange rate) [1/h]
        overall_ref_level = np.ceil(np.mean(self.ref_level))
        if overall_ref_level == 1:
            v_usage = self.V_RES[0][0][self.year_class][self.size_class]
            v_inf = self.V_RES[1][0][self.year_class][self.size_class]
        elif overall_ref_level == 2:
            v_usage = self.V_RES[0][1][self.year_class][self.size_class]
            v_inf = self.V_RES[1][1][self.year_class][self.size_class]
        elif overall_ref_level == 3:
            v_usage = self.V_RES[0][2][self.year_class][self.size_class]
            v_inf = self.V_RES[1][2][self.year_class][self.size_class]

        # Thermal capacitance [J/(K m2)]
        # Impact of refurbishment is not considered
        c_roof = self.C_RES[0][self.year_class][self.size_class]
        c_wall = self.C_RES[1][self.year_class][self.size_class]
        c_floor = self.C_RES[2][self.year_class][self.size_class]

        self.u = [u_roof, u_wall, u_floor, u_window]
        self.v = [v_usage, v_inf]
        self.c = [c_roof, c_wall, c_floor]

    #
    def get_building_thermal_properties_per_unit_area_non_residential(self):
        """
        Gets building thermal properties based on:
        >> source missing

        returns:
            u	<list>	[u_roof, u_wall, u_floor, u_window] in W/(K m2)
            v	<list>  [v_usage, v_infiltration] in 1/h
            c	<list>	[c_roof, c_wall, c_floor] in J/(K m2)
        """

        # Transmission losses (U-values) [W/(K m2)]

        ### ROOF
        if self.ref_level[0] == 1:
            u_roof = self.U_NRES[0][0][self.year_class]
        elif self.ref_level[0] == 2:
            u_roof = self.U_NRES[0][1][self.year_class]
        elif self.ref_level[0] == 3:
            u_roof = self.U_NRES[0][2][self.year_class]

        ### WALL
        if self.ref_level[1] == 1:
            u_wall = self.U_NRES[1][0][self.year_class]
        elif self.ref_level[1] == 2:
            u_wall = self.U_NRES[1][1][self.year_class]
        elif self.ref_level[1] == 3:
            u_wall = self.U_NRES[1][2][self.year_class]

        ### FLOOR
        if self.ref_level[2] == 1:
            u_floor = self.U_NRES[2][0][self.year_class]
        elif self.ref_level[2] == 2:
            u_floor = self.U_NRES[2][1][self.year_class]
        elif self.ref_level[2] == 3:
            u_floor = self.U_NRES[2][2][self.year_class]

        ### WINDOW
        if self.ref_level[3] == 1:
            u_window = self.U_NRES[3][0][self.year_class]
        elif self.ref_level[3] == 2:
            u_window = self.U_NRES[3][1][self.year_class]
        elif self.ref_level[3] == 3:
            u_window = self.U_NRES[3][2][self.year_class]

        # Ventilation losses (air exchange rate) [1/h]
        # Impact of refurbishment is not considered
        v_usage = self.V_NRES[0][self.year_class]
        v_inf = self.V_NRES[1][self.year_class]

        # Thermal capacitance [J/(K m2)]
        # Impact of refurbishment is not considered
        c_roof = self.C_NRES[0][self.year_class]
        c_wall = self.C_NRES[1][self.year_class]
        c_floor = self.C_NRES[2][self.year_class]

        self.u = [u_roof, u_wall, u_floor, u_window]
        self.v = [v_usage, v_inf]
        self.c = [c_roof, c_wall, c_floor]

    #
    def calculate_building_envelope_areas_residential(self):
        """
        Calculates building envelope areas (wall, roof and window).
        Residential: areas are calculated according to building typologies in TABULA.
        Only the heated area is considered.
        """

        # Floor-to-floor height
        # from http://www.ctbuh.org/HighRiseInfo/TallestDatabase/Criteria/HeightCalculator/tabid/1007/language/en-GB/Default.aspx
        # self.f2f_height = (randint(30, 32) / 10.0) # m
        # from TABULA
        self.f2f_height = 2.5  # m

        self.env_areas = np.zeros(4)

        # Building areas
        ### Floor
        self.env_areas[2] = self.heated_area / self.floors
        ### Roof
        self.env_areas[0] = self.ARATIO[0][self.year_class][self.size_class] * self.env_areas[2]
        ### Wall
        self.env_areas[1] = self.ARATIO[1][self.year_class][self.size_class] * self.env_areas[2] * self.free_walls
        ### Window
        self.env_areas[3] = self.ARATIO[2][self.year_class][self.size_class] * self.env_areas[2]

    #
    def calculate_building_envelope_areas_non_residential(self):
        """
        Calculates building envelope areas (wall, roof and window.
        Non-residential: number of floors and window-to-wall ratio are
                     derived from statistics and used to calculate the building areas.
        Only the heated area is considered.
        """

        self.env_areas = np.zeros(4)

        # Floor-to-floor height
        # from http://www.ctbuh.org/HighRiseInfo/TallestDatabase/Criteria/HeightCalculator/tabid/1007/language/en-GB/Default.aspx
        self.f2f_height = (randint(30, 39) / 10.0)  # m

        # Building areas
        ### Floor
        self.env_areas[2] = self.heated_area / self.floors
        ### Roof
        self.env_areas[0] = self.env_areas[2]
        ### Wall
        width = np.sqrt(self.env_areas[2])  # [m] building is assumed to be a cube
        self.env_areas[1] = self.free_walls * self.floors * self.f2f_height * width
        ### Window
        # >> Add Window-to-Wall ratio for the different building types
        # >> Source missing
        self.env_areas[3] = (randint(1, 4) / 10.0) * self.env_areas[1]

    #
    def calculate_building_window_areas(self, ):
        """
        Calculate window areas in each direction to calculate solar gains.
        """

        # Window areas oriented to [east, south, west, north]
        self.window_areas = np.zeros(4)

        # RESIDENTIAL
        if self.building.use == 3:
            for j in range(4):
                self.window_areas[j] = self.building.window_area * self.WRATIO_ORIENTATION[
                    5 * j + int(self.building.size_class)]

        # NON-RESIDENTIAL
        else:
            for j in range(4):
                self.window_areas[j] = self.building.window_area * self.WRATIO_ORIENTATION[5 * j + 4]

    #
    def calculate_building_thermal_properties(self):
        """
        Calculates equivalent U-value, thermal capacitance (C) and time constant (Tau)
        for the building. These properties are used in the first order thermal model.
        """

        # multiply u and c by envelope areas
        # consider adjustment factor for border situation acording to TABULA
        b_tr = 0.5  # adjustment factor for surface bordering on soil
        self.U_roof = self.u[0] * self.env_areas[0]
        self.U_wall = self.u[1] * self.env_areas[1]
        self.U_floor = self.u[2] * self.env_areas[2] * b_tr
        self.U_window = self.u[3] * self.env_areas[3]

        self.C_roof = self.c[0] * self.env_areas[0]
        self.C_wall = self.c[1] * self.env_areas[1]
        self.C_floor = self.c[2] * self.env_areas[2]

        # calculate ventilation losses, V [W/K]
        # air_specific_heat_capacity * air_density = 0.34 Wh/(m3 K)
        # air flow rate related to usage is considered as half the value in TABULA (0.4)
        self.V = (self.v[0] + self.v[1]) * self.env_areas[2] * self.f2f_height * 0.34

        # Equivalent building thermal properties
        ### Transmission losses (U-value), U [W/K]
        self.U = self.U_roof + self.U_wall + self.U_floor + self.U_window

        # self.adjust_building_thermal_properties(ref_level, size_class, U)

        ### Thermal capacitance, C [J/K]
        self.C = self.C_floor + self.C_wall + self.C_roof

        # Time constant, Tau [s]
        if self.U != 0:
            self.Tau = self.C / (self.U + self.V)
        else:
            self.Tau = np.nan

    #
    def adjust_building_thermal_properties(self):
        """
        Empirical adjustment of U-values to match TABULA results
        """

        if np.mean(self.ref_level) >= 1. and np.mean(self.ref_level) < 2.:
            if self.size_class == 0:  # SFH
                self.U = self.U * 1.
            elif self.size_class == 1:  # TH
                self.U = self.U * 1.2
            elif self.size_class == 2:  # MFH
                self.U = self.U * 1.18
            elif self.size_class == 3:  # AB
                self.U = self.U * 1.18
        elif np.mean(self.ref_level) >= 2. and np.mean(self.ref_level) < 3.:
            if self.size_class == 0:  # SFH
                self.U = self.U * 1.35
            elif self.size_class == 1:  # TH
                self.U = self.U * 1.31
            elif self.size_class == 2:  # MFH
                self.U = self.U * 1.53
            elif self.size_class == 3:  # AB
                self.U = self.U * 1.5
        elif np.mean(self.ref_level) >= 3.:
            if self.size_class == 0:  # SFH
                self.U = self.U * 1.6
            elif self.size_class == 1:  # TH
                self.U = self.U * 1.65
            elif self.size_class == 2:  # MFH
                self.U = self.U * 1.7
            elif self.size_class == 3:  # AB
                self.U = self.U * 1.55

    #
    def calculate_building_Tset(self):
        """
        Derives a target temperature by choosing a random temperature from Tset_mean
        +/- dT. Values differ for different building types.
        From http://tc76.org/spc100/docs/IBP%2018599/18599-10.pdf
        """

        if self.use == 0:  # commercial
            Tset_mean = self.TSET[0][0]
            dT = self.TSET[1][0]
        elif self.use == 1:  # industrial
            Tset_mean = self.TSET[0][1]
            dT = self.TSET[1][1]
        elif self.use == 2:  # public
            Tset_mean = self.TSET[0][2]
            dT = self.TSET[1][2]
        elif self.use == 3:  # residential
            Tset_mean = self.TSET[0][3]
            dT = self.TSET[1][3]

        self.Tset = float(randint(Tset_mean - dT, Tset_mean + dT))

    #
    def calculate_building_active_hours(self):
        """
        Assigns random start and end hours for building active hours.
        Values differ for different building types.

        returns:
            self.active_hours	<list>		[(start0, end0), (start1, end1)] in h
        """

        ### Residential
        if self.use == 3:
            self.active_hours = [0, 23]

        ### Non-residential
        else:
            if self.use == 0:  # commercial
                start_mean = self.SCHEDULE[0][0][0]
                end_mean = self.SCHEDULE[0][0][1]
                dt = self.SCHEDULE[1][0]
            elif self.use == 1:  # industrial
                start_mean = self.SCHEDULE[0][1][0]
                end_mean = self.SCHEDULE[0][1][1]
                dt = self.SCHEDULE[1][1]
            elif self.use == 2:  # public
                start_mean = self.SCHEDULE[0][2][0]
                end_mean = self.SCHEDULE[0][2][1]
                dt = self.SCHEDULE[1][2]

            start = randint(start_mean - dt, start_mean + dt)
            end = randint(end_mean - dt, end_mean + dt)
            self.active_hours = [start, end]

    ## Building parameters: Hot water demand
    #
    def calculate_daily_hot_water_demand(self):
        """
        Returns the daily hot water demand by getting a random value
        from the cdf function based on the statistics from VDI 3807-3
        (specific dhw demand in m3/m2 of living area)
        """

        # calculate specific dhw demand
        rand_num = np.random.uniform(0, 1, 1)
        specific_dhw = self.dhw_prob[0](rand_num)[0]

        # calculate building dhw demand
        self.hw_demand_day = specific_dhw * self.heated_area

    #
    def parametrize_hot_water_tank(self, X=1.5):
        """
        Calculates size and initial state of hot water tank.
        Size is X times the calculated daily demand.
        """

        # set tank capacity as function of daily hot water demand limit
        self.hw_tank_cap = X * self.hw_demand_day

    #
    def set_hot_water_tank_initial_state(self):
        """
        """
        # define initial state of hot water tank as random percentage of capacity
        self.hw_tank_volume_t0 = np.random.uniform(0, 1, 1)[0] * self.building.hw_tank_cap

    # Building occupancy
    # --------------------------------------------------------------------------------
    def calculate_building_activity_occupancy_vector(self):
        """
        Calculate vector of activity in building, i.e. percentage of occupied dwellings (for space heating)
        Active_hours (scheduled), building occupancy and weekends are considered
        """

        # initialize matrices
        self.activity_vector = np.ones([self.nts]) * 1.  # Building activity vector
        self.occupancy_vector = np.ones([self.nts]) * 1.  # Number of occupants in building per timestep

        # calculate occupant schedule
        if self.building.use == 3:  # residential
            self.calculate_occupants_schedule()

        # calculate activity vector
        for iii in range(0, self.nts):
            hour = self.dt_vector[iii].hour
            wd = self.dt_vector[iii].weekday()  # 0 Monday ... 6 Sunday

            if self.building.use == 0:  # commercial
                if self._workday_weekend:
                    if wd == 6:  # no activities on Sunday
                        self.activity_vector[iii] = 0.
                else:
                    if hour < self.building.act_start or hour > self.building.act_end:
                        self.activity_vector[iii] = 0.

                # calculate occupancy vector
                self.occupancy_vector[iii] = self.activity_vector[iii] * self.building.occupants

            elif self.building.use == 1 or self.building.use == 2:  # industrial/public
                if self._workday_weekend:
                    if wd >= 5:  # no activities on Saturday/Sunday
                        self.activity_vector[iii] = 0.
                else:
                    if hour < self.building.act_start or hour > self.building.act_end:
                        self.activity_vector[iii] = 0.

                # calculate occupancy vector
                self.occupancy_vector[iii] = self.activity_vector[iii] * self.building.occupants

            elif self.building.use == 3:  # residential

                # calculate building occupancy considering active population if the flag active_pop is True
                if self._active_population:
                    # initialize counters
                    not_in_building = 0
                    if self._workday_weekend:
                        if wd < 5:  # weekends are considered as active
                            # occupant loop
                            for occupant in self.occupant_vector:
                                # check if occupant is in the building
                                if hour > occupant[1][0][1] and hour < occupant[1][1][0]:  # not at home
                                    not_in_building += 1
                    else:
                        # occupant loop
                        for occupant in self.occupant_vector:
                            # check if occupant is in the building
                            if hour > occupant[1][0][1] and hour < occupant[1][1][0]:  # not at home
                                not_in_building += 1
                    # check occupancy share
                    self.occupancy_vector[iii] = self.building.occupants - not_in_building
                    self.activity_vector[iii] = (self.building.occupants - not_in_building) / self.building.occupants

                else:
                    # calculate occupancy vector
                    self.occupancy_vector[iii] = self.activity_vector[iii] * self.building.occupants

    #
    def calculate_occupants_schedule(self):
        """
        A schedule is assigned to every occupant based on studying/working schedule.

        returns:
            occupant_vector	<list>		[dwelling, [occupant, [schedule]]] for occupant and dwelling in building
        """

        # calculate occupancy based on percentage of active population
        perc_working_pop = 0.583  # percentage of population that works or studies
        # from https://ergebnisse.zensus2011.de/#StaticContent:091840148148,BEG_1_6_1,m,table

        # occupant vector
        O = int(self.building.occupants)
        self.occupant_vector = [[[], []] for occupant in range(O)]

        for iii in range(O):

            # add occupant to vector
            self.occupant_vector[iii][0] = iii + 1

            # compare random number with percentage of working/student population
            # to determine if occupants are at home
            rand_num = np.random.uniform(0, 1, 1)
            if rand_num > perc_working_pop:  # occupant leaves the building to work/study
                # calculate working/studying start/end time from "public" building usage
                start_mean = self.SCHEDULE[0][2][0]
                end_mean = self.SCHEDULE[0][2][1]
                dt = self.SCHEDULE[1][2]

                start = randint(start_mean - dt, start_mean + dt)
                end = randint(end_mean - dt, end_mean + dt)

                # define active_hours per occupant
                self.occupant_vector[iii][1] = [[0, start], [end, 23]]
            else:
                self.occupant_vector[iii][1] = [[0, 23], [0, 0]]

    # Results
    # --------------------------------------------------------------------------------
    def calculate_annual_demand(self, data):
        """
        Calculate the annual energy demand by weighting the heating demand of typical days
        """

        # Divide data into days
        ## initialize variables
        days_in_sim = self.number_of_typ_days
        data_in_days = np.zeros((days_in_sim, (24 * 60) // self.resolution), dtype=float)
        energy = np.zeros(days_in_sim, dtype=float)
        row_start = np.zeros(days_in_sim, dtype=int)
        row_end = np.zeros(days_in_sim, dtype=int)

        ## arrange data
        for day in range(days_in_sim):
            if day == 0:
                row_start[day] = 0
                row_end[day] = 24 * 60 / self.resolution
            else:
                row_start[day] = row_end[day - 1]
                row_end[day] = 24 * 60 / self.resolution + row_start[day]

            data_in_days[day, :] = np.reshape(data[row_start[day]:row_end[day]], ((24 * 60) // self.resolution,))

            # add energy per day and multiply by weight of typical day
            energy[day] = data_in_days[day, :].sum() * self.resolution * 1 / 60 * self.weights[day]

        # calculate annual demand
        annual_demand = energy.sum()

        return annual_demand

    #
    def plot_timeseries(self,
                        space_heating=True, Tb=False,
                        hot_water=True,
                        total=True,
                        xticks=('month', 3)):
        """
        Plots heat demand timeseries
        """

        if not os.path.exists(self.result_dir):
            os.makedirs(self.result_dir)

        if space_heating:
            fig_name = '{}/SpaceHeatingDemand_{}_{}_{}_{}.png'.format(self.result_dir,
                                                                      self.building.use, self.building.year_class,
                                                                      self.building.size_class, self.building.bid)

            UrbanHeatPro.plot_timeseries(self.dt_vector,
                                         [self.space_heating_power], ['Space heating demand'], fig_name,
                                         ynumticks='auto', ylabel='Power [kW]', ylim0=True, yfactor=1e3,
                                         xticks=xticks)

        if Tb:
            fig_name = '{}/BuildingTemperature_{}_{}_{}_{}.png'.format(self.result_dir,
                                                                       self.building.use, self.building.year_class,
                                                                       self.building.size_class, self.building.bid)

            UrbanHeatPro.plot_timeseries(self.dt_vector,
                                         [self.Tamb, self.Tb], ['T_amb', 'T_b'], fig_name,
                                         ynumticks='auto', ylabel='Temperature [degC]', ylim0=False, yfactor=1,
                                         xticks=xticks)

        if hot_water:
            fig_name = '{}/HotWaterDemand_{}_{}_{}_{}.png'.format(self.result_dir,
                                                                  self.building.use, self.building.year_class,
                                                                  self.building.size_class, self.building.bid)

            UrbanHeatPro.plot_timeseries(self.dt_vector,
                                         [self.hot_water_power], ['Hot water demand'], fig_name,
                                         ynumticks='auto', ylabel='Power [kW]', ylim0=True, yfactor=1e3,
                                         xticks=xticks)

        if total:
            fig_name = '{}/TotalHeatDemand_{}_{}_{}_{}.png'.format(self.result_dir,
                                                                   self.building.use, self.building.year_class,
                                                                   self.building.size_class, self.building.bid)

            UrbanHeatPro.plot_stacked_timeseries(self.dt_vector,
                                                 [self.hot_water_power, self.space_heating_power],
                                                 ['Hot water', 'Space heating'], fig_name,
                                                 ynumticks='auto', ylabel='Power [kW]', ylim0=True, yfactor=1e3,
                                                 xticks=xticks)

    #
    def save_csv(self):
        """
        Saves key building parameters and heat demand (space heating, hot water and
        total) as timeseries.
        """

        if not os.path.exists(self.result_dir):
            os.makedirs(self.result_dir)

        #
        try:
            self.hot_water_power *= 1.
        except:
            self.hot_water_power = np.zeros([self.nts])
            self.hot_water_tank_m3 = np.zeros([self.nts])
        try:
            self.total_power_delayed *= 1.
        except:
            self.total_power_delayed = np.zeros([self.nts])

        array_to_save = np.array([self.dt_vector_excel, self.activity_vector,
                                  self.Tamb,
                                  self.Tb,
                                  self.occupancy_vector,
                                  self.space_heating_power, self.internal_gains, self.solar_gains,
                                  self.hot_water_power, self.hot_water_tank_m3,
                                  self.total_power,
                                  self.total_power_delayed]).transpose()

        filename = '{}/HeatDemand_{}_{}_{}_{}.csv'.format(self.result_dir,
                                                          self.building.use, self.building.year_class,
                                                          self.building.size_class, self.building.bid)

        # Column names
        with open(filename, 'w') as text_file:
            text_file.write(
                'datenum;active?;Tamb_degC;Tb_degC;Occupants;SpaceHeatingD_W;int_gains_W;sol_gains_W;HotWaterD_W'
                ';HotWaterTank_m3;TotalD_W;DelayedTotalD_W\n')

        # Time series
        with open(filename, 'a') as f:
            np.savetxt(f, array_to_save, delimiter=';', fmt='%.2f')

    #
    def save_load_duration_curve(self):
        """
        Save sorted demand
        """

        filename = '{}/LoadDurationCurve_{}_{}_{}_{}.csv'.format(self.result_dir,
                                                                 self.use, self.year_class, self.size_class, self.bid)

        # sort demand to save
        array_to_save = np.sort(self.total_power_delayed)[::-1].transpose()

        with open(filename, 'a') as f:
            np.savetxt(f, array_to_save, delimiter=';', fmt='%.4f')

    #
    def save_dhw_debug_csv(self):
        """
        Saves debug values for dhw demand.
        """

        array_to_save = np.array(self.dhw_debug)

        filename = '{}/DHWDemand_{}_{}_{}_{}.csv'.format(self.result_dir,
                                                         self.use, self.year_class, self.size_class, self.bid)

        # Building parameters
        with open(filename, 'w') as text_file:
            text_file.write(
                'timstep;hour;day;day_in_year;daily_dhw;daily_limit;activity;shower?;bath?;medium?;small?;fr_shower'
                ';fr_bath;fr_medium;fr_small;d_shower;d_bath;d_medium;d_small;c_shower;c_bath;c_medium;c_small'
                ';daily_agg;daily_convergence\n')

        # Time series
        with open(filename, 'a') as f:
            np.savetxt(f, array_to_save, delimiter=';', fmt='%.4f')
