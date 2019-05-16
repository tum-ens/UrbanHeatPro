"""
runme.py
AMC @ TUM ENS
"""

import multiprocessing

import UrbanHeatPro.Functions as UrbanHeatPro
import UrbanHeatPro.Classes as UrbanHeatPro

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
region				= 'Unterhaching'
N 	 	   			= 1		  		
resolution 			= 60
(offset, length) 	= (0, 24*365)
timesteps 			= range(offset, offset + length)	  		
number_of_typ_days  = 365
#
sce_refurbishment   = None
sce_Tamb		    = None
#
processes			= 24
chunk_size			= 2439
###
SIMULATION			= [[region], 
					   [N, resolution, timesteps, number_of_typ_days],
					   [sce_refurbishment, sce_Tamb], 
					   [processes, chunk_size]]


# 2. CITY
#	 2.1 Raw building data
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
#filename_buildings  = None
filename_buildings  = 'buildings_Unterhaching.csv'
#
filename_syn_city	= None
#filename_syn_city	= 'SynCity_Unterhaching_0.csv'
#
connection_factor	= 1.
#
_space_heating 		= True
_hot_water			= True
_energy_only		= False
#
base_load			= 5. * 1e6
###
CITY				= [[filename_buildings, filename_syn_city], [connection_factor],
					   [_space_heating, _hot_water, _energy_only], [base_load]]

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
#_trans_losses 		= True
#_ventilation_losses = True
_internal_gains		= True
_solar_gains		= True	
_active_population  = True
_workday_weekend	= True
_monthly_sh_prob	= True
#
Tb0_str				= 'Tset'
dTset				= 0.			
#
eta					= 1.0
dT_per_hour			= 0.1	
thermal_inertia		= 0.4
#
_night_set_back		= 0.5
schedule_nsb	    = [23, 6]	
T_nsb				= 19
power_reduction		= 0
#
###
SPACE_HEATING		= [[_internal_gains, _solar_gains, _active_population, _workday_weekend, _monthly_sh_prob],
					   [None],
					   [Tb0_str, dTset],
					   [eta, dT_per_hour, thermal_inertia],
					   [_night_set_back, schedule_nsb, T_nsb, power_reduction]]

					   
# 4. HOT WATER DEMAND
#	 4.1 Hot water temperature
#		 Tw 			<float> Hot water temperature in degC
#    4.2 Hot water tank
# 		 hw_tank_limit 	<float> Lower limit of hot water tank as decimal.
#									Below this limit, hw tank is refilled.
#		 hw_flow 		<float> Volume flow to refill hot water tank in L/min
Tw				    = 60
#
hw_tank_limit		= 0.1
hw_flow				= 15
###
HOT_WATER			= [[Tw], [hw_tank_limit, hw_flow]]


# 5. REPORTING
# 	0 No results saved or plotted
#	1 Results per simulation
#   2 Results per building
# 	3 Results per time step
plot 				= 0
save 				= 1
debug 				= 1
###
REPORTING			= [plot, save, debug]

# MAIN
# --------------------------------------------------------------------------------
if __name__ == '__main__':
	
	# Simulation name
	NAME 	  = '{}_0'.format(region)
							  
	multiprocessing.freeze_support()
	my_Simulation = UrbanHeatPro.Simulation(NAME, SIMULATION, CITY, SPACE_HEATING, HOT_WATER, REPORTING)
	my_Simulation.run()