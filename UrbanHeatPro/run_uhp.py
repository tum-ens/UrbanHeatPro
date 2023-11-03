"""
runme.py
A. Molar-Cruz @ TUM ENS
"""

import multiprocessing
import os

import sorcery
import yaml

from UrbanHeatPro.Classes import Simulation
from UrbanHeatPro.Functions.uhp_utils import access_config_or_default, nested_get


# MAYBE convert to class?
def run_uhp(selected_region: str = None, simulation_name: str = None, buildings_use_filter="",
            settings_file="../settings/uhp_settings_example.yaml", result_dir=None):
    # CONTENT
    # ----------------------------------------------------------------------------------------------------
    # 1. SIMULATION
    #	 1.1 General
    #	 1.2 Scenarios
    #	 1.3 Multiprocessing
    # 2. CITY
    #	 2.1 Building data
    # 	 2.2 Connection factor
    # 	 2.3 City heat demand
    #	 2.4 Base load
    # 3. SPACE HEATING DEMAND
    #	 3.1 Flags
    #	 3.2 Refurbishment level
    #	 3.3 Initial temperature
    #	 3.4 Heating system
    #	 3.5 Demand Side Management
    # 4. HOT WATER DEMAND
    #    5.1 Hot water temperature
    #    5.2 Hot water tank
    # 5. REPORTING

    # INPUT DATA
    # ----------------------------------------------------------------------------------------------------
    default_config_path = "../settings/uhp_default_settings.yaml"
    with open(os.path.abspath(os.path.join(os.path.dirname(__file__), default_config_path)), 'r') as f:
        default_config: dict = yaml.load(f, Loader=yaml.FullLoader)
    with open(settings_file, 'r') as f:
        config: dict = yaml.load(f, Loader=yaml.FullLoader)

    # 1. SIMULATION
    #	 1.1 General
    #		 region 	<str>	name of region/city/urban area
    #	 	 N 			<int>	number of simulation runs. One run calculates the heat demand for the whole region.
    #	 	 resolution <int> 	temporal resolution [min]
    #	 	 offset 	<int> 	initial time step
    #	 	 length 	<int> 	number of time steps to simulate
    #		 timesteps	<int>	vector of time steps to simulate
    #		 number_of_typ_days <int> number of typical days to simulate
    # 	 1.2 Scenarios
    #		 sce_refurbishment	<str> Name of refurbishment scenario or None
    #		 filename_Tamb		<str> Name of ambient temperature scenario or None
    #	 1.3 Multiprocessing
    #		 processes 	<int> 	number of processes to use for multiprocessing (parallelization)
    #		 chunk_size <int>	number of buildings in chunk to save
    #

    region = selected_region if selected_region is not None else access_config_or_default(config, default_config,
                                                                                          ["simulation", "region"])
    N = access_config_or_default(config, default_config, ["simulation", "general", "N"])
    resolution = access_config_or_default(config, default_config, ["simulation", "general", "resolution"])
    offset = access_config_or_default(config, default_config, ["simulation", "general", "offset"])
    length = access_config_or_default(config, default_config, ["simulation", "general", "length"])
    timesteps = range(offset, offset + length)
    number_of_typ_days = access_config_or_default(config, default_config,
                                                  ["simulation", "general", "number_of_typ_days"])
    #
    sce_refurbishment = access_config_or_default(config, default_config,
                                                 ["simulation", "scenarios", "sce_refurbishment"])
    sce_Tamb = access_config_or_default(config, default_config, ["simulation", "scenarios", "sce_Tamb"])
    #
    tmp = nested_get(config, ["simulation", "multi_processing", "processes"], None)
    processes = multiprocessing.cpu_count() if tmp is None else tmp
    # number of lines in input file todo automate this parameter?
    chunk_size = access_config_or_default(config, default_config, ["simulation", "multi_processing", "chunk_size"])
    ###
    SIMULATION = [[region],
                  [N, resolution, timesteps, number_of_typ_days],
                  [sce_refurbishment, sce_Tamb],
                  [processes, chunk_size]]

    # to dict for yaml
    general = sorcery.dict_of(N, resolution, offset, length, number_of_typ_days)
    scenarios = sorcery.dict_of(sce_refurbishment, sce_Tamb)
    multi_processing = sorcery.dict_of(processes, chunk_size)
    simulation = sorcery.dict_of(region, general, scenarios, multi_processing)

    # 2. CITY
    # 	 2.1 Raw building data
    #		 filename_buildings	<str>		name of csv file with raw building data or None.
    #	 2.3 Synthetic city
    #	 	 filename_syn_city 	<str>		name of csv file with synthetic city or None.
    #	 2.2 Connection factor
    #		 connection_factor 	<float>		share of buildings connected to the network (as decimal)
    #	 2.3 City heat demand
    #		 _space_heating 	<boolean>	specifies if space heating demand is calculated.
    #								   		If False, heat losses and gains are also False.
    #		 _hot_water 		<boolean>	specifies if hot water demand is calculated.
    #		 _energy_only 		<boolean>	specifies if the focus is only on the aggregated heating energy
    # 								 		demand and not on the time series. If True, the hot water demand is
    #								 		calculated/added per day and not per time step.
    #	 2.4 Base load
    #		 base_load 			<float>		minimum load at every time step in W.
    #

    if buildings_use_filter == "":
        filename_buildings = 'buildings_{}.csv'.format(region)
    else:
        filename_buildings = 'buildings_{}_{}.csv'.format(region, buildings_use_filter)
    # filename_buildings  = None
    #
    filename_syn_city = access_config_or_default(config, default_config, ["city", "building_data", "filename_syn_city"])
    #
    connection_factor = access_config_or_default(config, default_config, ["city", "connection_factor"])
    #
    _space_heating = access_config_or_default(config, default_config, ["city", "city_heat_demand", "_space_heating"])
    _hot_water = access_config_or_default(config, default_config, ["city", "city_heat_demand", "_hot_water"])
    _energy_only = access_config_or_default(config, default_config, ["city", "city_heat_demand", "_energy_only"])
    #
    base_load = access_config_or_default(config, default_config, ["city", "base_load"])
    ###
    CITY = [[filename_buildings, filename_syn_city], [connection_factor],
            [_space_heating, _hot_water, _energy_only], [base_load]]

    # to dict for yaml
    building_data = sorcery.dict_of(filename_buildings, filename_syn_city)
    city_heat_demand = sorcery.dict_of(_space_heating, _hot_water, _energy_only)
    city = sorcery.dict_of(building_data, connection_factor, city_heat_demand, base_load)

    # 3. SPACE HEATING DEMAND
    #	 3.1 Flags
    #		 _trans_losses 			<boolean>	specifies if transmission losses (wall, windows, roof and floor)
    #										    are calculated.
    #		 _ventilation_losses 	<boolean>	specifies if ventilation and infiltration losses are considered.
    #		 _internal_gains 		<boolean>	specifies if internal gains are calculated.
    #		 _solar_gains 			<boolean>	specifies if solar losses are calculated.
    #		 _active_population 	<boolean>	specifies if statistics for active population are used to create
    #									   		synthetic population profiles (occupancy)
    #		 _workday_weekend		<boolean>   specifies if workdays and weekends are differentiated.
    #							  				Use only when full year is simulated.
    #		 _monthly_sh_prob		<boolean>	specifies if monthly probability of using heating is used
    #    3.2 Temperature
    # 		 Tb0_str 				<str>		Building initial temperature as string 'Tset' or 'Tamb'.
    #		 dTset 					<float>		Temperature difference to modify Tset_min, Tset_max in degC
    #	 3.3 Heating system
    #		 eta 					<float>		heating system efficiency
    #		 thermal_intertia 		<float>		weight of the delivered power from previous time step.
    #								   			i.e. how much can the output power change with respect to the previous time step?
    #	     dT_per_hour 			<float>		maximum temperature difference allowed in the building in degC / h.
    #	 3.4 Demand Side Management
    #		 _night_set_back 		<float>		share of buildings with night set-back
    #		 schedule_nsb 			<list>		[start, end] of night set-back in hours
    #		 T_nsb 					<float>		night set-back temperature in degC
    # 		 power_reduction 		<float>		reduced power as decimal. Input power = 1 - power_reduction
    # _trans_losses 		= True
    # _ventilation_losses = True
    _internal_gains = access_config_or_default(config, default_config,
                                               ["space_heating_demand", "flags", "_internal_gains"])
    _solar_gains = access_config_or_default(config, default_config, ["space_heating_demand", "flags", "_solar_gains"])
    _active_population = access_config_or_default(config, default_config,
                                                  ["space_heating_demand", "flags", "_active_population"])
    _workday_weekend = access_config_or_default(config, default_config,
                                                ["space_heating_demand", "flags", "_workday_weekend"])
    _monthly_sh_prob = access_config_or_default(config, default_config,
                                                ["space_heating_demand", "flags", "_monthly_sh_prob"])
    #
    Tb0_str = access_config_or_default(config, default_config, ["space_heating_demand", "temperature", "Tb0_str"])
    dTset = access_config_or_default(config, default_config, ["space_heating_demand", "temperature", "dTset"])
    #
    eta = access_config_or_default(config, default_config, ["space_heating_demand", "heating_system", "eta"])
    dT_per_hour = access_config_or_default(config, default_config,
                                           ["space_heating_demand", "heating_system", "dT_per_hour"])
    thermal_inertia = access_config_or_default(config, default_config,
                                               ["space_heating_demand", "heating_system", "thermal_inertia"])
    #
    _night_set_back = access_config_or_default(config, default_config,
                                               ["space_heating_demand", "demand_side_management", "_night_set_back"])
    schedule_nsb = access_config_or_default(config, default_config,
                                            ["space_heating_demand", "demand_side_management", "schedule_nsb"])
    T_nsb = access_config_or_default(config, default_config,
                                     ["space_heating_demand", "demand_side_management", "T_nsb"])
    power_reduction = access_config_or_default(config, default_config,
                                               ["space_heating_demand", "demand_side_management", "power_reduction"])
    #
    ###
    SPACE_HEATING = [[_internal_gains, _solar_gains, _active_population, _workday_weekend, _monthly_sh_prob],
                     [None],  # TODO !!! what is here?
                     [Tb0_str, dTset],
                     [eta, dT_per_hour, thermal_inertia],
                     [_night_set_back, schedule_nsb, T_nsb, power_reduction]]

    # to dict for yaml
    flags = sorcery.dict_of(_internal_gains, _solar_gains, _active_population, _workday_weekend, _monthly_sh_prob)
    empty = {"empty": None}
    temperature = sorcery.dict_of(Tb0_str, dTset)
    heating_system = sorcery.dict_of(eta, dT_per_hour, thermal_inertia)
    demand_side_management = sorcery.dict_of(_night_set_back, schedule_nsb, T_nsb, power_reduction)
    space_heating_demand = sorcery.dict_of(flags, empty, temperature, heating_system, demand_side_management)

    # 4. HOT WATER DEMAND
    #	 4.1 Hot water temperature
    #		 Tw 			<float> Hot water temperature in degC
    #    4.2 Hot water tank
    # 		 hw_tank_limit 	<float> Lower limit of hot water tank as decimal.
    #									Below this limit, hw tank is refilled.
    #		 hw_flow 		<float> Volume flow to refill hot water tank in L/min
    Tw = access_config_or_default(config, default_config, ["hot_water_demand", "Tw"])
    #
    hw_tank_limit = access_config_or_default(config, default_config,
                                             ["hot_water_demand", "hot_water_tank", "hw_tank_limit"])
    hw_flow = access_config_or_default(config, default_config, ["hot_water_demand", "hot_water_tank", "hw_flow"])
    ###
    HOT_WATER = [[Tw], [hw_tank_limit, hw_flow]]

    # to dict for yaml
    hot_water_tank = sorcery.dict_of(hw_tank_limit, hw_flow)
    hot_water_demand = sorcery.dict_of(Tw, hot_water_tank)

    # 5. REPORTING
    # 	0 No results saved or plotted
    #   1 Results per simulation
    #   2 Results per building --> as in the UHP_output_profile.csv example
    # 	3 Results per time step
    plot = access_config_or_default(config, default_config, ["reporting", "plot"])
    save = access_config_or_default(config, default_config, ["reporting", "save"])
    debug = access_config_or_default(config, default_config, ["reporting", "debug"])

    result_dir = access_config_or_default(config, default_config, ["reporting", "result_dir"]) if result_dir is None else result_dir

    ###
    REPORTING = [plot, save, debug, result_dir]

    # to dict for yaml
    reporting = sorcery.dict_of(plot, save, debug, result_dir)

    # MAIN
    # --------------------------------------------------------------------------------

    sim_dict = sorcery.dict_of(simulation, city, space_heating_demand, hot_water_demand, reporting)
    used_config_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../settings/uhp_settings_currently_used.yaml'))
    with open(used_config_path, 'w') as f:
        yaml.dump(sim_dict, f, sort_keys=False)
        print("Loaded and used configuration for the current simulation run saved at", used_config_path)

    # ----- Run the simulation with the given settings
    # Simulation name
    NAME = simulation_name if simulation_name not in [None, ''] else '{}_0'.format(region)
    multiprocessing.freeze_support()
    my_Simulation = Simulation(NAME, SIMULATION, CITY, SPACE_HEATING, HOT_WATER, REPORTING)
    my_Simulation.run()
