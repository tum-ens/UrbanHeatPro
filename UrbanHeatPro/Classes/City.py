"""
City.py
A. Molar-Cruz @ TUM ENS
"""

import multiprocessing
import sys

import numpy as np
from scipy import interpolate

import UrbanHeatPro.Functions as UrbanHeatPro
from .Building import Building


# import ipdb


class City:
    # --------------------------------------------------------------------------------
    def __init__(self, NAME, SIMULATION, CITY, SPACE_HEATING, HOT_WATER, REPORTING):

        # SIMULATION
        self.name = NAME  # Simulation name
        self.region = SIMULATION[0][0]  # Region
        self.sce = SIMULATION[0][1]  # Scenario
        self.dt_vector = SIMULATION[1][0]  # Vector of time steps as datetime objects
        self.dt_vector_excel = SIMULATION[1][1]  # Vector of time steps as excel date
        self.nts = len(self.dt_vector)  # Number of time steps
        self.resolution = SIMULATION[2][0]  # Temporal resolution in min
        self.number_of_typ_days = SIMULATION[4][0]  # Number of typical days
        self.weights = SIMULATION[4][1]  # Weights of typical days

        # PARALELLIZATION
        self.processes = SIMULATION[3][0]  # Number of parallel processes
        self.chunk_size = SIMULATION[3][1]  # Number of buildings in chunk to save
        self.b_to_save_syncity = np.zeros([self.chunk_size, 30])
        self.b_to_save_heat = np.zeros([self.chunk_size, 40])
        self.counter_syncity = 0
        self.counter_heat = 0

        # CITY
        ### Ambient conditions
        self.Tamb = CITY[0][0]  # Ambient temperature vector in degC
        self.I = CITY[0][1]  # Solar radiation vector in W/m2 [I_Gh, I_Dh, I_ex, hs]
        ### Building data
        self.buildings = CITY[1][0]  # Building data
        self.building_stock_stats = CITY[1][1]  # Building data from statistics
        self.nb = len(self.buildings)  # Number of buildings
        self.connection_factor = CITY[2][0]  # Share of buildings connected to the network
        ### Flags
        self._space_heating = CITY[3][0]  # calculate space heating demand?
        self._hot_water = CITY[3][1]  # calculate hot water demand?
        self._energy_only = CITY[3][2]  # calculate only aggregated demand?
        ### Base load
        self.base_load = np.ones([self.nts]) * CITY[4][0]  # Vector of base load in W

        # SPACE HEATING DEMAND
        ### Flags
        self._internal_gains = SPACE_HEATING[0][0]  # consider internal gains?
        self._solar_gains = SPACE_HEATING[0][1]  # consider solar gains?
        self._active_population = SPACE_HEATING[0][2]  # consider active population for occupancy vector?
        self._workday_weekend = SPACE_HEATING[0][3]  # consider difference between workdays and weekends?
        self._monthly_sh_prob = SPACE_HEATING[0][4]  # consider monthly probability of using heating?
        ### Transmission losses
        self.refurbishment_level = SPACE_HEATING[1][0]  # Refurbishment level for all buildings
        self.Tb0_str = SPACE_HEATING[2][0]  # Initial building temperature as 'ambient' or 'Tset'?
        self.dTset = SPACE_HEATING[2][1]  # Delta temperature (for Tset_min, Tset_max)
        ### Heating system
        self.eta = SPACE_HEATING[3][0]  # Heating process efficiency
        self.dT_per_hour = SPACE_HEATING[3][1]  # Maximum dT allowed in building per hour [degC]
        self.thermal_inertia = SPACE_HEATING[3][2]  # Thermal inertia of the heating system
        ### DSM
        self._night_set_back = SPACE_HEATING[4][0]  # Share of buildings with night set-back
        self.schedule_nsb = SPACE_HEATING[4][1]  # [start, end] of nsb in h
        self.T_nsb = SPACE_HEATING[4][2]  # Night set-back temperature in degC
        self.power_reduction = SPACE_HEATING[4][3]  # Percentage of power reduced (as decimal)

        # HOT WATER DEMAND
        ### Calculation
        self.Tw = HOT_WATER[0][0]  # Hot water temperature in degC
        self.hw_tank_limit = HOT_WATER[1][0]  # Hot water tank limit as perc (decimal)
        self.hw_flow = HOT_WATER[1][1]  # Flow to refill hot water tank in L/min
        self.dhw_prob = self.initialize_dhw_probabilities()  # Probabilities for calculation od dhw demand

        # RESULTS
        ### Runs
        self.rid = SIMULATION[5][0]  # Run id
        self.result_dir = SIMULATION[5][1]  # Directory where results are stored

        # REPORTING
        self.plot = REPORTING[0]  # Plot level  [0, 1, 2]
        self.save = REPORTING[1]  # Save level  [0, 1, 2]
        self.debug = REPORTING[2]  # Debug level [0, 1, 2]

    # Synthetic city
    # --------------------------------------------------------------------------------
    def create_synthetic_city(self):
        """
        Create a synthetic city representing the building stock based on statistics.
        """

        # extract building stock statistics
        ## Residential
        self.FOOTPRINT = self.building_stock_stats[0][0]  # Footprint area of building typologies (TABULA)
        self.FLOORS = self.building_stock_stats[0][3]  # Number of floors
        self.CURRENT_REF_RES = self.building_stock_stats[0][7][0]  # Percentage of residential refurbished buildings
        self.SINGLE_DWELLING = self.building_stock_stats[0][8]  # Percentage of single dwellings for SFH and TH
        self.AVG_DWELLING_SIZE = self.building_stock_stats[0][9]  # Average dwelling size in m2
        self.HOUSEHOLD_SIZE = self.building_stock_stats[0][10]  # Household size for dwelling size categories
        self.STOCK_RES = self.building_stock_stats[0][17]  # Building stock statistics for residential
        ## Non-residential
        self.CURRENT_REF_NRES = self.building_stock_stats[0][14][0]  # Percentage of residential refurbished buildings
        self.STOCK_NRES = self.building_stock_stats[0][18]  # Building stock statistics for non-residential

        """
        # DEBUG
        for iii in range(1000):
            building = self.buildings.loc[iii, :]
            result = self.create_synthetic_building(building)
            print(iii)
            print(result)
        input()
        """

        # save synthetic city file
        # filename = '{}/SynCity_{}_{}_{}.csv'.format(self.result_dir, self.region, self.sce, self.rid)
        filename = '{}/SynCity_{}_{}.csv'.format(self.result_dir, self.name, self.rid)
        self.save_csv_syn_city_header(filename)

        # Multiprocessing
        if self.debug >= 1:
            print('      ***')
            print('      Creating synthetic city with {} buildings...'.format(len(self.buildings)))
            print('         Starting multiprocessing with {}/{} processors...'.format(self.processes,
                                                                                      multiprocessing.cpu_count()))

        ### Writer queue
        if self.debug >= 1:
            print('         Starting queues...')
        writerQueue = multiprocessing.Queue()
        writeProc = multiprocessing.Process(target=self.write_to_synthetic_city, args=(writerQueue, filename))

        ### Feeder queue
        buildings_list = range(len(self.buildings))
        feederQueue = multiprocessing.Queue()
        feedProc = multiprocessing.Process(target=self.feed_building_to_process, args=(feederQueue, buildings_list))

        ### Processes
        calcProc = [multiprocessing.Process(target=self.call_create_synthetic_building,
                                            args=(feederQueue, writerQueue)) for iii in range(self.processes)]

        ### start multiprocessing
        if self.debug >= 1:
            print('         Calculating synthetic building...')
        feedProc.start()
        for p in calcProc:
            p.start()
        if self.debug >= 1:
            print('\n         Saving to file...')
        writeProc.start()

        feedProc.join()
        for p in calcProc:
            p.join()
        writeProc.join()

    #
    def update_synthetic_city(self, ref_matrix_res, ref_matrix_nres):
        """
        Create a synthetic city representing the building stock based on statistics.
        """

        # values to update
        self.ref_matrix_res = ref_matrix_res
        self.ref_matrix_nres = ref_matrix_nres

        """
        # DEBUG
        for iii in range(1000):
            building = self.buildings.loc[iii, :]
            result = self.update_synthetic_building(building)
            print(iii)
            print(result)
            input()
        input()
        """

        # save updated synthetic city file
        # filename = '{}/SynCity_{}_{}_{}.csv'.format(self.result_dir, self.region, self.sce, self.rid)
        filename = '{}/SynCity_{}_{}.csv'.format(self.result_dir, self.name, self.rid)
        self.save_csv_syn_city_header(filename)

        # Multiprocessing
        if self.debug >= 1:
            print('      ***')
            print('      Updating synthetic city...')
            print('         Starting multiprocessing with {}/{} processors...'.format(self.processes,
                                                                                      multiprocessing.cpu_count()))

        ### Writer queue
        if self.debug >= 1:
            print('         Starting queues...')
        writerQueue = multiprocessing.Queue()
        writeProc = multiprocessing.Process(target=self.write_to_synthetic_city, args=(writerQueue, filename))

        ### Feeder queue
        buildings_list = range(len(self.buildings))
        feederQueue = multiprocessing.Queue()
        feedProc = multiprocessing.Process(target=self.feed_building_to_process, args=(feederQueue, buildings_list))

        ### Processes
        calcProc = [multiprocessing.Process(target=self.call_update_synthetic_building,
                                            args=(feederQueue, writerQueue)) for iii in range(self.processes)]

        ### start multiprocessing
        if self.debug >= 1:
            print('         Updating synthetic building...')
        feedProc.start()
        for p in calcProc:
            p.start()
        if self.debug >= 1:
            print('\n         Saving to file...')
        writeProc.start()

        feedProc.join()
        for p in calcProc:
            p.join()
        writeProc.join()

    #
    def feed_building_to_process(self, feederQueue, buildings_list):
        """
        Feeds building data to Queue
        """

        for iii in buildings_list:
            building = self.buildings.loc[iii, :]
            feederQueue.put((building, iii))

    #
    def call_create_synthetic_building(self, feederQueue, writerQueue):
        """
        Calls function to create synthetic building
        """

        while True:
            try:
                # get from Queue
                (building, iii) = feederQueue.get(block=True, timeout=10)

                # calculate
                result = self.create_synthetic_building(building)

                # write to Queue
                writerQueue.put((result, iii))
                # print('         {} | bid: {}'.format(multiprocessing.current_process().name, result[0][0]))

                if self.debug >= 1:
                    # Progress bar
                    sys.stdout.write('\r')
                    i = float((iii + 1) / len(self.buildings))
                    sys.stdout.write("         Calculation:   [%-20s] %d%%" % ('|' * int(20 * i), 100 * i))
                    sys.stdout.flush()

            except Exception as e:
                # print('>>> Error: {}'.format(e))
                break

    #
    def call_update_synthetic_building(self, feederQueue, writerQueue):
        """
        Calls function to update synthetic building.
        To update:
            - Refurbishment level
        """

        while True:
            try:
                # get from Queue
                (building, iii) = feederQueue.get(block=True, timeout=10)

                # calculate
                result = self.update_synthetic_building(building)

                # write to Queue
                writerQueue.put((result, iii))
                # print('         {} | bid: {}'.format(multiprocessing.current_process().name, result[0][0]))

                if self.debug >= 1:
                    # Progress bar
                    sys.stdout.write('\r')
                    i = float((iii + 1) / len(self.buildings))
                    sys.stdout.write("         Calculation:   [%-20s] %d%%" % ('|' * int(20 * i), 100 * i))
                    sys.stdout.flush()

            except Exception as e:
                # print('>>> Error: {}'.format(e))
                break

    #
    def create_synthetic_building(self, building):
        """
        Creates building object and calculates missing building properties
        """

        my_building = Building(
            # Building data
            building, self.building_stock_stats,
            # Simulation
            [None, None], None,
            None, None,
            None, None,
            None, None, None,
            # Space heating demand
            None, None, None, None, None,
            None, None, None,
            None, None,
            None, None, None, None,
            # Hot water demand
            None, self.dhw_prob, None, None,
            None, None, None,
            # Results
            self.result_dir, None, None, None)

        r = my_building.parametrize_building()

        result = [r.bid, r.footprint_area, r.use, r.free_walls,
                  r.lat, r.lon, r.dist_to_heat_source,
                  r.year_class, r.size_class, r.floors,
                  r.storey_area, r.area_corr_factor, r.heated_area, r.window_area,
                  r.dwellings, r.dwelling_size, r.occupants,
                  r.ref_level_roof, r.ref_level_wall, r.ref_level_floor, r.ref_level_window,
                  r.U, r.V, r.C, r.Tau,
                  r.Tset, r.act_start, r.act_end,
                  r.hw_demand_day, r.hw_tank_cap]
        result = np.resize(result, (1, len(result)))

        return result

    #
    def update_synthetic_building(self, building):
        """
        Updates synthetic building.
        To update:
            - Refurbishment level
        """

        my_building = Building(
            # Building data
            building, self.building_stock_stats,
            # Simulation
            [None, None], None,
            None, None,
            None, None,
            None, None, None,
            # Space heating demand
            None, None, None, None, None,
            None, None, None,
            None, None,
            None, None, None, None,
            # Hot water demand
            None, None, None, None,
            None, None, None,
            # Results
            self.result_dir, None, None, None)

        r = my_building.update_building_refurbishment_level(self.ref_matrix_res, self.ref_matrix_nres)

        result = [r.bid, r.footprint_area, r.use, r.free_walls,
                  r.lat, r.lon, r.dist_to_heat_source,
                  r.year_class, r.size_class, r.floors,
                  r.storey_area, r.area_corr_factor, r.heated_area, r.window_area,
                  r.dwellings, r.dwelling_size, r.occupants,
                  r.ref_level_roof, r.ref_level_wall, r.ref_level_floor, r.ref_level_window,
                  r.U, r.V, r.C, r.Tau,
                  r.Tset, r.act_start, r.act_end,
                  r.hw_demand_day, r.hw_tank_cap]
        result = np.resize(result, (1, len(result)))

        return result

    #
    def write_to_synthetic_city(self, writerQueue, filename):
        """
        Write synthetic building results to file
        """

        while True:
            try:
                # get from Queue
                (result, iii) = writerQueue.get(block=True, timeout=20)

                # add building to matrix to save
                if self.counter_syncity < self.chunk_size:
                    self.b_to_save_syncity[self.counter_syncity, :] = result
                    self.counter_syncity += 1

                # save results once chunk size is reached
                else:

                    # save
                    f = open(filename, 'ab')
                    fmt = ['%d', '%.2f', '%d', '%d',
                           '%.3f', '%.3f', '%.2f',
                           '%d', '%d', '%d',
                           '%.2f', '%.2f', '%.2f', '%.2f',
                           '%d', '%.2f', '%d',
                           '%d', '%d', '%d', '%d',
                           '%.2f', '%.2f', '%.2f', '%.2f',
                           '%d', '%d', '%d',
                           '%.2f', '%.2f']
                    np.savetxt(f, self.b_to_save_syncity, delimiter=';', fmt=fmt)
                    f.close()

                    # restart counter
                    self.counter_syncity = 0
                    self.b_to_save_syncity = np.zeros([self.chunk_size, 30])

                    # Progress bar
                    # print('      bid: {}'.format(result[0][0]))
                    if self.debug >= 1:
                        sys.stdout.write('\r')
                        i = float((iii + 1) / len(self.buildings))
                        sys.stdout.write("\n            Saving:        [%-20s] %d%%" % ('|' * int(20 * i), 100 * i))
                        sys.stdout.flush()

            except Exception as e:
                # print('>>> Error: {}'.format(e))

                # save again to not loose data
                f = open(filename, 'ab')
                fmt = ['%d', '%.2f', '%d', '%d',
                       '%.3f', '%.3f', '%.2f',
                       '%d', '%d', '%d',
                       '%.2f', '%.2f', '%.2f', '%.2f',
                       '%d', '%.2f', '%d',
                       '%d', '%d', '%d', '%d',
                       '%.2f', '%.2f', '%.2f', '%.2f',
                       '%d', '%d', '%d',
                       '%.2f', '%.2f']

                # drop zero rows
                self.b_to_save_syncity = self.b_to_save_syncity[~np.all(self.b_to_save_syncity == 0, axis=1)]
                np.savetxt(f, self.b_to_save_syncity, delimiter=';', fmt=fmt)
                f.close()

                break

    # City heat demand
    # --------------------------------------------------------------------------------
    def calculate_city_heat_demand(self):
        """
        Paralellizes the calculation of heating energy demand per building using a
        given number of processes. Every process modifies a shared dictionary where
        the heat demand is stored as power and energy.
        """

        if self.debug >= 1:
            print('      ***')
            print('      Calculating city heat demand...')
            print('         Starting multiprocessing with {}/{} available processors'.format(self.processes,
                                                                                             multiprocessing.cpu_count()))

        # save heat demand file
        # filename = '{}/EnergyPerBuilding_{}_{}_{}.csv'.format(self.result_dir, self.region, self.sce, self.rid)
        filename = '{}/EnergyPerBuilding_{}_{}.csv'.format(self.result_dir, self.name, self.rid)
        self.save_csv_energy_header(filename)

        # Initialize probabilities for space heating and dhw
        self.sh_prob = None
        self.day_vector = None
        self.seasonal_vector = None
        self.min_vector = self.calculate_min_vector()  # Vector of minutes in time frame

        #
        if self._space_heating:
            if self._monthly_sh_prob:
                SH_PROB = self.building_stock_stats[0][19]  # Monthly probability of using heating
                self.sh_prob = [SH_PROB[dt.month - 1] for dt in self.dt_vector]
            else:
                self.sh_prob = [1. for dt in self.dt_vector]
        #
        if self._hot_water:
            self.day_vector = self.calculate_day_vector()  # Vector of days in time frame
            if not self._energy_only:
                self.seasonal_vector = self.calculate_seasonal_variation_vector()  # Sine vector to represent seasonal variations

        """
        # DEBUG
        self.sh_power = np.zeros([self.nts])
        self.sh_energy = 0.
        self.dhw_power = np.zeros([self.nts])
        self.dhw_energy = 0.
        self.total_power = np.zeros([self.nts])
        self.total_energy = 0.
        for iii in range(1000):
            print(iii)
            building = self.buildings.loc[iii, :]
            result = self.calculate_building_heat_demand(building)
            print(result)
        input()	
        """

        # Calculation
        # ----------------------------------------------------------------
        # City manager
        # this dictionary is shared with all the processes
        if self.debug >= 1:
            print('         Starting multiprocessing manager...')
        manager = multiprocessing.Manager()

        if self._energy_only:
            self.results = manager.dict(energy_sh=0.,
                                        energy_hw=0.,
                                        energy_total=0.)
        else:
            self.results = manager.dict(power_sh=np.zeros([self.nts]),
                                        energy_sh=0.,
                                        power_hw=np.zeros([self.nts]),
                                        energy_hw=0.,
                                        power_total=np.zeros([self.nts]),
                                        power_total_delayed=np.zeros([self.nts]),
                                        energy_total=0.)

        # Queues
        if self.debug >= 1:
            print('         Starting queues...')

        ## Feeder queue
        buildings_list = range(len(self.buildings))
        feederQueue = multiprocessing.Queue()
        feedProc = multiprocessing.Process(target=self.feed_building_to_process, args=(feederQueue, buildings_list))

        ## Writer queue
        writerQueue = multiprocessing.Queue()
        writeProc = multiprocessing.Process(target=self.write_to_city_heat_demand, args=(writerQueue, filename))

        # Processes
        calcProc = [multiprocessing.Process(target=self.call_calculate_building_heat_demand,
                                            args=(feederQueue, writerQueue)) for iii in range(self.processes)]

        ## start
        # if (self.debug >= 1):
        # print('      Calculating building heat demand...')
        feedProc.start()
        for p in calcProc:
            p.start()
        # if (self.debug >= 1):
        # print('\n      Saving to file...')
        writeProc.start()

        ## join
        feedProc.join()
        for p in calcProc:
            p.join()
        writeProc.join()

        # add base load to timeseries
        if not self._energy_only:
            self.results['power_total'] += self.base_load
            self.results['energy_total'] += np.sum(self.base_load)

        # save results per city
        if self.save >= 1:
            if not self._energy_only:
                self.save_csv_power()
            self.save_csv_energy()

    #
    def call_calculate_building_heat_demand(self, feederQueue, writerQueue):
        """
        Calls function to calculate building heat demand
        """

        while True:
            try:
                # get from Queue
                (building, iii) = feederQueue.get(block=True, timeout=10)

                # calculate
                result = self.calculate_building_heat_demand(building)

                # write to Queue
                writerQueue.put((result, iii))
                # print('         {} | bid: {}'.format(multiprocessing.current_process().name, result[0][0]))

                if self.debug >= 1:
                    # Progress bar
                    sys.stdout.write('\r')
                    i = float((iii + 1) / len(self.buildings))
                    sys.stdout.write("         Calculation:   [%-20s] %d%%" % ('|' * int(20 * i), 100 * i))
                    sys.stdout.flush()

            except Exception as e:
                # print('>>> Error: {}'.format(e))
                break

    #
    def calculate_building_heat_demand(self, building):
        """

        Extracts building information needed to create a Building object.
        If the building is connected to the district heating network, then a Building
        object is created and the heat demand is calculted. If it is not, then the
        heat demand is set to zero.

        args:
            building	dataframe with building information
            iii			building counter
        """

        # initialize results
        sh = np.zeros(6)
        hw = np.zeros(2)
        total = np.zeros(2)

        # check if building is connected to the network
        rand_num = np.random.uniform(0, 1, 1)[0]
        if rand_num < self.connection_factor:

            my_building = self.create_building_object(building)

            # Calculations
            ## Space heating demand
            if self._space_heating:
                sh = my_building.calculate_space_heating_demand()

                ### add building space heating demand to city space heating demand
                if not self._energy_only:
                    # self.sh_power += my_building.space_heating_power
                    self.results['power_sh'] += my_building.space_heating_power
                # self.sh_energy += my_building.space_heating_energy
                self.results['energy_sh'] += my_building.space_heating_energy

            ## Hot water demand
            if self._hot_water:
                hw = my_building.calculate_hot_water_demand()

                ### add building hot water demand to city hot water demand
                if not self._energy_only:
                    # self.dhw_power += my_building.hot_water_power
                    self.results['power_hw'] += my_building.hot_water_power
                # self.dhw_energy += my_building.dhw_energy
                self.results['energy_hw'] += my_building.dhw_energy

            ## Total heat demand
            total = my_building.calculate_total_heat_demand()

            ## add building demand to city demand
            if not self._energy_only:
                # self.total_power += my_building.total_power
                self.results['power_total'] += my_building.total_power
            # self.total_energy += my_building.total_energy
            self.results['energy_total'] += my_building.total_energy

            ## save results per building
            if self.save >= 2:
                my_building.save_csv()
            # my_building.save_load_duration_curve()

            ## save results per timestep
            if self.save == 3:
                if hot_water:
                    my_building.save_dhw_debug_csv()

            # Plot
            if self.plot == 2:
                my_building.plot_timeseries(space_heating=self._space_heating, Tb=True,
                                            hot_water=self._hot_water, total=True,
                                            xticks=('month', 3))

        r = building
        syn_b = [r.bid, r.footprint_area, r.use, r.free_walls,
                 r.lat, r.lon, r.dist_to_heat_source,
                 r.year_class, r.size_class, r.floors,
                 r.storey_area, r.area_corr_factor, r.heated_area, r.window_area,
                 r.dwellings, r.dwelling_size, r.occupants,
                 r.ref_level_roof, r.ref_level_wall, r.ref_level_floor, r.ref_level_window,
                 r.U, r.V, r.C, r.Tau,
                 r.Tset, r.act_start, r.act_end,
                 r.hw_demand_day, r.hw_tank_cap]
        result = [syn_b, sh, hw, total]

        # return flattened list	with one row
        result = [item for sublist in result for item in sublist]
        result = np.resize(result, (1, len(result)))

        return result

    #
    def create_building_object(self, building):
        """
        Creates instance of class Object
        """

        # create Building object
        my_building = Building(
            # Building data
            building, self.building_stock_stats,
            # Simulation
            [self.dt_vector, self.dt_vector_excel], self.resolution,
            self.number_of_typ_days, self.weights,
            self.Tamb, self.I,
            self._space_heating, self._hot_water, self._energy_only,
            # Space heating demand
            self.Tb0_str, self.dTset, self.dT_per_hour, self.eta, self.thermal_inertia,
            self._active_population, self._workday_weekend, self.sh_prob,
            self._solar_gains, self._internal_gains,
            self._night_set_back, self.schedule_nsb, self.T_nsb, self.power_reduction,
            # Hot water demand
            self.Tw, self.dhw_prob, self.hw_tank_limit, self.hw_flow,
            self.day_vector, self.seasonal_vector, self.min_vector,
            # Results
            self.result_dir, self.plot, self.save, self.debug)

        return my_building

    #
    def write_to_city_heat_demand(self, writerQueue, filename):
        """
        Writes building properties and heat demand to file
        """

        while True:
            try:
                # get from Queue
                (result, iii) = writerQueue.get(block=True, timeout=20)

                # add building to matrix to save
                if self.counter_heat < self.chunk_size:
                    self.b_to_save_heat[self.counter_heat, :] = result
                    self.counter_heat += 1

                # save results once chunk size is reached
                else:
                    # save
                    f = open(filename, 'ab')
                    fmt = ['%d', '%.2f', '%d', '%d',
                           '%.3f', '%.3f', '%.2f',
                           '%d', '%d', '%d',
                           '%.2f', '%.2f', '%.2f', '%.2f',
                           '%d', '%.2f', '%d',
                           '%d', '%d', '%d', '%d',
                           '%.2f', '%.2f', '%.2f', '%.2f',
                           '%d', '%d', '%d',
                           '%.2f', '%.2f',
                           '%.2f', '%.2f', '%.2f',
                           '%.2f', '%.2f', '%.2f',
                           '%.2f', '%.2f',
                           '%.2f', '%.2f']
                    np.savetxt(f, self.b_to_save_heat, delimiter=';', fmt=fmt)
                    f.close()

                    # restart counter
                    self.counter_heat = 0
                    self.b_to_save_heat = np.zeros([self.chunk_size, 40])

                    # Progress bar
                    # print('      bid: {}'.format(result[0][0]))
                    if self.debug >= 1:
                        sys.stdout.write('\r')
                        i = float((iii + 1) / len(self.buildings))
                        sys.stdout.write("\n            Saving:        [%-20s] %d%%" % ('|' * int(20 * i), 100 * i))
                        sys.stdout.flush()

            except Exception as e:
                # print('>>> Error: {}'.format(e))

                # save again to not loose data
                f = open(filename, 'ab')
                fmt = ['%d', '%.2f', '%d', '%d',
                       '%.3f', '%.3f', '%.2f',
                       '%d', '%d', '%d',
                       '%.2f', '%.2f', '%.2f', '%.2f',
                       '%d', '%.2f', '%d',
                       '%d', '%d', '%d', '%d',
                       '%.2f', '%.2f', '%.2f', '%.2f',
                       '%d', '%d', '%d',
                       '%.2f', '%.2f',
                       '%.2f', '%.2f', '%.2f',
                       '%.2f', '%.2f', '%.2f',
                       '%.2f', '%.2f',
                       '%.2f', '%.2f']

                # drop zero rows
                self.b_to_save_heat = self.b_to_save_heat[~np.all(self.b_to_save_heat == 0, axis=1)]
                np.savetxt(f, self.b_to_save_heat, delimiter=';', fmt=fmt)
                f.close()

                break

    # Domestic hot water demand
    # --------------------------------------------------------------------------------
    def initialize_dhw_probabilities(self):
        """
        Calculates dhw probabilities (daily consumption, event loads, flow rate and
        duration as interpolate objects.
        """

        # Daily specific dhw consumption [m3/m2 of living area]	as cdf
        x = [i[0] for i in self.building_stock_stats[1][0]]
        p = [i[1] for i in self.building_stock_stats[1][0]]
        dhw_cdf = UrbanHeatPro.create_interpolated_cdf(x, p)

        # Probability distribution of the DHW-load happening at specific time of the day
        #
        ### Probabilities in 0.2h-steps
        p_shower = [i[1] for i in self.building_stock_stats[1][1]]
        p_bath = [i[2] for i in self.building_stock_stats[1][1]]
        p_med_small = [i[3] for i in self.building_stock_stats[1][1]]
        #
        ### Interpolators
        x_min = 0.
        x_max = 23.8
        x_steps = 120
        x = np.linspace(x_min, x_max, num=x_steps, endpoint=False)
        shower_ptime = interpolate.interp1d(x, p_shower)
        bath_ptime = interpolate.interp1d(x, p_bath)
        medium_ptime = interpolate.interp1d(x, p_med_small)
        small_ptime = medium_ptime
        #
        prob_time = [shower_ptime, bath_ptime, medium_ptime, small_ptime]

        # Factor for probability distribution of the DHW-load happening at specific weekday
        shower_pwday = [i[1] for i in self.building_stock_stats[1][2]]
        bath_pwday = [i[2] for i in self.building_stock_stats[1][2]]
        medium_pwday = [i[3] for i in self.building_stock_stats[1][2]]
        small_pwday = medium_pwday
        prob_wday = [shower_pwday, bath_pwday, medium_pwday, small_pwday]

        # Load flow rate CDF
        x_min = 3.
        x_max = 20.
        x_steps = 1000
        x = np.linspace(x_min, x_max, num=x_steps, endpoint=False)
        ### Shower
        mean = self.building_stock_stats[1][3][0][0]
        sigma = self.building_stock_stats[1][3][1][0]
        norm_dist = UrbanHeatPro.create_normal_distribution(x, mean, sigma)
        shower_f_cdf = UrbanHeatPro.create_interpolated_cdf(x, norm_dist)
        #
        ### Bath
        mean = self.building_stock_stats[1][3][0][1]
        sigma = self.building_stock_stats[1][3][1][1]
        norm_dist = UrbanHeatPro.create_normal_distribution(x, mean, sigma)
        bath_f_cdf = UrbanHeatPro.create_interpolated_cdf(x, norm_dist)
        #
        x_min = 0.
        x_max = 10.
        x_steps = 1000
        x = np.linspace(x_min, x_max, num=x_steps, endpoint=False)
        ### Medium
        mean = self.building_stock_stats[1][3][0][2]
        sigma = self.building_stock_stats[1][3][1][2]
        norm_dist = UrbanHeatPro.create_normal_distribution(x, mean, sigma)
        medium_f_cdf = UrbanHeatPro.create_interpolated_cdf(x, norm_dist)
        #
        ### Small
        mean = self.building_stock_stats[1][3][0][3]
        sigma = self.building_stock_stats[1][3][1][3]
        norm_dist = UrbanHeatPro.create_normal_distribution(x, mean, sigma)
        small_f_cdf = UrbanHeatPro.create_interpolated_cdf(x, norm_dist)
        #
        flowRate_cdf = [shower_f_cdf, bath_f_cdf, medium_f_cdf, small_f_cdf]

        # Load duration CDF
        x_min = 3.
        x_max = 15.
        x_steps = 1000
        x = np.linspace(x_min, x_max, num=x_steps, endpoint=False)
        ### Shower
        mean = self.building_stock_stats[1][4][0][0]
        sigma = self.building_stock_stats[1][4][1][0]
        norm_dist = UrbanHeatPro.create_normal_distribution(x, mean, sigma)
        shower_d_cdf = UrbanHeatPro.create_interpolated_cdf(x, norm_dist)
        #
        ### Bath
        mean = self.building_stock_stats[1][4][0][1]
        sigma = self.building_stock_stats[1][4][1][1]
        norm_dist = UrbanHeatPro.create_normal_distribution(x, mean, sigma)
        bath_d_cdf = UrbanHeatPro.create_interpolated_cdf(x, norm_dist)
        #
        x_min = 0.
        x_max = 5.
        x_steps = 1000
        x = np.linspace(x_min, x_max, num=x_steps, endpoint=False)
        ### Medium
        mean = self.building_stock_stats[1][4][0][2]
        sigma = self.building_stock_stats[1][4][1][2]
        norm_dist = UrbanHeatPro.create_normal_distribution(x, mean, sigma)
        medium_d_cdf = UrbanHeatPro.create_interpolated_cdf(x, norm_dist)
        #
        ### Small
        mean = self.building_stock_stats[1][4][0][3]
        sigma = self.building_stock_stats[1][4][1][3]
        norm_dist = UrbanHeatPro.create_normal_distribution(x, mean, sigma)
        small_d_cdf = UrbanHeatPro.create_interpolated_cdf(x, norm_dist)
        #
        duration_cdf = [shower_d_cdf, bath_d_cdf, medium_d_cdf, small_d_cdf]

        dhw_prob = [dhw_cdf, prob_time, prob_wday, flowRate_cdf, duration_cdf]

        return dhw_prob

    #
    def calculate_seasonal_variation_vector(self, amplitude=0.1, max_point=45):
        """
        Creates a sine wave representing the change of the nominal consumption during
        the year due to the seasonal variation.

        Args:
            amplitude		<float>		Variation of consumption (% of nominal load)
            max_point		<int>		Day in year with the highest hot water consumption
                                        (lowest ambient temperature)

        Returns:
            seasonal_vector	<numpy array>

        """

        # initialize seasonal vector
        seasonal_vector = np.zeros(366)

        # create sine function
        for iii in range(1, 367):
            perc = float(iii) / 366.
            seasonal_vector[iii - 1] = np.sin(2. * np.pi * (perc + max_point / 366.)) * amplitude + 1.

        return seasonal_vector

    #
    def calculate_day_vector(self):
        """
        Calculates a vector of the days in the year included in the simulation
        time frame. Maximum length is 366.

        Returns:
            self.day_vector	<list>	list of day numbers in simulation time frame
                                    with start and end indices
        """

        num_days = self.number_of_typ_days

        # initialize day_vector
        day_vector = np.zeros([num_days, 3], dtype=int)

        # day_vector: [day in year, start time step]
        prev_day = 0
        day_count = 0
        for iii, date in enumerate(self.dt_vector):
            if not (date.timetuple().tm_yday == prev_day):
                day_vector[day_count][0] = date.timetuple().tm_yday
                day_vector[day_count][1] = iii  # index of dt_vector where a new day starts
                day_count += 1
            prev_day = date.timetuple().tm_yday

        # add end time step
        for iii, day in enumerate(day_vector):
            try:
                day_vector[iii][2] = day_vector[iii + 1][1] - 1
            except:  # last day
                day_vector[iii][2] = len(self.dt_vector) - 1

        return day_vector

    #
    def calculate_min_vector(self):
        """
        Calculates a vector of the simulation time steps in minutes of year. Maximum length is
        366*24*60.

        Returns:
            self.min_vector	<list>	list of time steps in minutes
        """

        # initialize min_vector
        min_vector = np.zeros([self.nts], dtype=int)

        # calculate minutes of year
        for iii, date in enumerate(self.dt_vector):
            day_in_year = date.timetuple().tm_yday
            hour_in_day = date.timetuple().tm_hour
            min_in_hour = date.timetuple().tm_min

            min_vector[iii] = (day_in_year - 1) * 24 * 60 + hour_in_day * 60 + min_in_hour

        return min_vector

    # Results
    # --------------------------------------------------------------------------------
    def plot_timeseries(self, space_heating=True, hot_water=True, total=True,
                        xticks=('month', 3)):
        """
        """

        if space_heating:
            fig_name = '{}/SpaceHeatingDemand_{}.png'.format(self.result_dir, self.rid)

            UrbanHeatPro.plot_timeseries(self.dt_vector,
                                         [self.space_heating_power], ['Space heating demand'],
                                         fig_name, xticks=xticks,
                                         ynumticks='auto', ylabel='Power [kW]', ylim0=True, yfactor=1e3)

        if hot_water:
            fig_name = '{}/HotWaterDemand_{}.png'.format(self.result_dir, self.rid)

            UrbanHeatPro.plot_timeseries(self.dt_vector,
                                         [self.hot_water_power], ['Hot water demand'],
                                         fig_name, xticks=xticks,
                                         ynumticks='auto', ylabel='Power [kW]', ylim0=True, yfactor=1e3)

        if total:
            fig_name = '{}/TotalHeatDemand_{}.png'.format(self.result_dir, self.rid)

            UrbanHeatPro.plot_stacked_timeseries(self.dt_vector,
                                                 [self.hot_water_power, self.space_heating_power],
                                                 ['Hot water', 'Space heating'],
                                                 fig_name, xticks=xticks,
                                                 ynumticks='auto', ylabel='Power [kW]', ylim0=True, yfactor=1e3)

    #
    def save_csv_syn_city_header(self, filename):
        """
        Saves key building parameters of every chunk.
        """

        f = open(filename, 'w')
        # Region
        f.write('{};{}'.format(self.region, self.name))
        # Building data
        f.write('\nbid;footprint_area;use;free_walls;' + \
                'lat;lon;dist_to_heat_source;' + \
                'year_class;size_class;floors;' + \
                'storey_area;area_corr_factor;heated_area;window_area;' + \
                'dwellings;dwelling_size;occupants;' + \
                'ref_level_roof;ref_level_wall;ref_level_floor;ref_level_window;' + \
                'U;V;C;Tau;' + \
                'Tset;act_start;act_end;' + \
                'hw_demand_day;hw_tank_cap\n')

        f.close()

    #
    def save_csv_energy_header(self, filename):
        """
        Saves key building parameters of every chunk.
        """

        f = open(filename, 'w')
        # Region
        f.write(self.region)
        # Building data
        f.write('\nbid;footprint_area;use;free_walls;' + \
                'lat;lon;dist_to_heat_source;' + \
                'year_class;size_class;floors;' + \
                'storey_area;area_corr_factor;heated_area;window_area;' + \
                'dwellings;dwelling_size;occupants;' + \
                'ref_level_roof;ref_level_wall;ref_level_floor;ref_level_window;' + \
                'U;V;C;Tau;' + \
                'Tset;act_start;act_end;' + \
                'hw_demand_day;hw_tank_cap;' + \
                'sh_demand;solar_gains;int_gains;' + \
                'sh_demand_per_h_area;solar_gains_per_h_area;int_gains_per_h_area;' + \
                'hw_demand;hw_demand_m3;' + \
                'total_demand;total_demand_per_h_area\n')

        f.close()

    #
    def save_csv_power(self):
        """
        Saves heat demand timeseries in csv file (space heating, hot water and total).
        """

        array_to_save = np.array([self.dt_vector_excel,
                                  self.Tamb,
                                  self.results['power_sh'],
                                  self.results['power_hw'],
                                  self.results['power_total'],
                                  self.results['power_total_delayed']]).transpose()

        # filename = '{}/CityHeatDemand-Profile_{}_{}_{}.csv'.format(self.result_dir, self.region, self.sce, self.rid)
        filename = '{}/CityHeatDemand-Profile_{}_{}.csv'.format(self.result_dir, self.name, self.rid)

        # Header
        with open(filename, 'w') as f:
            f.write('Run id;{};\n'.format(self.rid))
            f.write('Total buildings;{};\n'.format(len(self.buildings)))
            f.write('Commercial;{};\n'.format(np.sum(self.buildings.use == 0)))
            f.write('Industrial;{};\n'.format(np.sum(self.buildings.use == 1)))
            f.write('Public;{};\n'.format(np.sum(self.buildings.use == 2)))
            f.write('Residential;{};\n'.format(np.sum(self.buildings.use == 3)))
            f.write('SpaceHeatingDemand_GWh;{};\n'.format(self.results['energy_sh'] / 1e9))
            f.write('HotWaterDemand_GWh;{};\n'.format(self.results['energy_hw'] / 1e9))
            f.write('Base Load_GWh;{};\n'.format(np.sum(self.base_load) / 1e9))
            f.write('TotalHeatDemand_GWh;{};\n'.format(self.results['energy_total'] / 1e9))
            f.write('\ndatenum;Tamb_degC;SpaceHeatingDemand_MW;HotWaterDemand_MW;' +
                    'TotalHeatDemand_MW;TotalHeatDemand_Delayed_MW\n')

        # Building data
        with open(filename, 'ab') as f:
            np.savetxt(f, array_to_save, delimiter=';', fmt='%.3f')

    #
    def save_csv_energy(self):
        """
        Saves aggregated heat demand in csv file (space heating, hot water and total).
        """

        array_to_save = np.array([self.results['energy_sh'],
                                  self.results['energy_hw'],
                                  self.results['energy_total']])

        array_to_save = np.resize(array_to_save, (1, len(array_to_save)))

        # filename = '{}/CityHeatDemand-Aggregate_{}_{}_{}.csv'.format(self.result_dir, self.region, self.sce, self.rid)
        filename = '{}/CityHeatDemand-Aggregate_{}_{}.csv'.format(self.result_dir, self.name, self.rid)

        f = open(filename, 'a')
        # Region
        f.write(self.region)
        # Building data
        f.write('\nSHD_Wh;DHWD_Wh;TOTAL_Wh\n')
        f.close()

        # Building data
        with open(filename, 'ab') as f:
            np.savetxt(f, array_to_save, delimiter=';', fmt='%.3f')
