"""
Building.py
AMC @ TUM ENS
"""

import os
import numpy as np
from datetime import datetime, timedelta
from random import randint, uniform
import scipy

import UrbanHeatPro.Functions as UrbanHeatPro
from .SpaceHeatingDemand import SpaceHeatingDemand
from .HotWaterDemand 	 import HotWaterDemand

class Building():
# --------------------------------------------------------------------------------
	def __init__(self, dt_vectors, resolution, Tamb, I, _space_heating, _hot_water, \
					Tb0_str, dTset, dT_per_hour, eta, thermal_inertia, building_stock_stats, \
					bid, use, gfa, free_walls, lat, lon, distance2hp, building_typology, year_class_int, btype_int, \
					refurbishment_level, _active_population, _solar_gains, _internal_gains, \
					_night_set_back, schedule_nsb, T_nsb, power_reduction, \
					Tw, dhw_prob, hw_tank_limit, hw_flow, \
					day_vector, seasonal_vector, min_vector, \
					result_dir, plot, save, debug):
		
		# General data
		# --------------------------------------------------
		# Building data
		self.bid                  = bid							# Building id
		self.use                  = UrbanHeatPro.building_use_to_tuple(use) # Building use (int, str)
		self.free_walls           = free_walls					# Number of walls in contact with Tamb
		self.coords				  = (lat, lon)					# Coordinates of building centroid
		self.building_typology	  = building_typology
		self.year_class_int		  = year_class_int
		self.btype_int			  = btype_int
		self.year_class           = None 						# Construction year class (int, str)
		self.btype				  = None 						# Building type (int, str)
		self.env_areas			  = [0, 0, gfa, 0] 				# Envelope areas [Roof, Walls, Floor (gfa), Window]
		self.window_areas		  = [0, 0, 0, 0]				# Window area oriented to the [east, south, west, north]
		self.floors				  = 0							# number of floors
		self.f2f_height			  = 0				    		# floor to floor height [m]
		self.heated_area		  = 0							# total heated area (reference area) [m2]
		self.living_area		  = 0							# total living area [m2]
		self.dwellings			  = 0							# number of dwellings in building
		self.dwelling_size_cat	  = 0							# dwelling size category to calculate household size
		
		# Simulation
		self.dt_vector    		  = dt_vectors[0] 				# Vector of time steps as datetime objects
		self.dt_vector_excel	  = dt_vectors[1] 				# Vector of time steps as excel date
		self.nts          		  = len(self.dt_vector)			# Number of time steps
		self.resolution			  = resolution					# Temporal resolution in min
		self._space_heating		  = _space_heating				# Calculate space heating demand?
		self._hot_water			  = _hot_water					# Calculate hot water demand?
	
		# External factors	
		self.Tamb         		  = Tamb		 				# Ambient temperature vector in degC
		self.I	         	  	  = I					 		# Solar radiation vector in W/m2 [I_Gh, I_Dh, I_ex, hs]
		self.eta          		  = eta			 				# Heating process efficiency
		self.thermal_inertia	  = thermal_inertia				# Thermal inertia of the heating system
		
		# Result directory
		self.result_dir			  = result_dir					# Directory where results are stored
		
		
		# Space heating demand
		# --------------------------------------------------
		# Building stock statistics I
		self.GFA  				  = building_stock_stats[0][0]  # GFA of building typologies (TABULA)
		self.ARATIO 			  = building_stock_stats[0][1]  # Area ratio [Roof, Wall, Window]
		self.WRATIO_ORIENTATION	  = building_stock_stats[0][2]  # Window ratio [East, South, West, North]
		self.FLOORS				  = building_stock_stats[0][3]  # Number of floors
		self.U_RES				  = building_stock_stats[0][4]  # U-values for residential buildings
		self.V_RES				  = building_stock_stats[0][5]  # Air flow rate (ventilation losses) for residential
		self.C_RES				  = building_stock_stats[0][6]  # Thermal cap for residential buildings
		self.REFURBISHED_RES	  = building_stock_stats[0][7]  # Percentage of residential refurbished buildings
		self.SINGLE_DWELLING	  = building_stock_stats[0][8]  # Percentage of single dwellings for SFH and TH
		self.HOUSEHOLD_SIZE		  = building_stock_stats[0][9]  # Household size for dwelling size categories
		self.U_NRES				  = building_stock_stats[0][10] # U-values for non-residential buildings
		self.V_NRES				  = building_stock_stats[0][11]	# Air flow rate (ventilation losses) for non-residential
		self.C_NRES				  = building_stock_stats[0][12]	# Thermal cap for non-residential buildings
		self.TSET				  = building_stock_stats[0][13] # Target temperatures per building use
		self.SCHEDULE			  = building_stock_stats[0][14]	# Active hours per building use
		self.STOCK_RES			  = building_stock_stats[0][15] # Building stock statistics for residential
		self.STOCK_NRES			  = building_stock_stats[0][16] # Building stock statistics for non-residential
		
		# Building thermal properties
		self.U                    = 0.							# Transmission losses (U-value) [W/K]
		self.V					  = 0.							# Ventilation losees [W/K]
		self.C                    = 0.							# Thermal capacitance [J/K]
		self.Tau                  = 0.							# Time constant [s]
		self.Tb0_str			  = Tb0_str						# Initial building temperature as string: 'ambient' or 'Tset'
		self.Tb0                  = 0.							# Initial building temperature
		self.Tset                 = 0.							# [Target temperature, dT]  [degC]
		self.dTset				  = dTset						# Delta temperature (for Tset_min, Tset_max)
		self.dT_per_hour		  = dT_per_hour					# Maximum dT allowed in building per hour [degC]
		
		# Activity and occupancy in building
		self._active_population	  = _active_population			# Consider active population for occupancy vector
		self.active_hours         = [[0, 0], [0, 0]]			# Building active hours [(start0, end0), (start1, end1)]
		self.activity_vector	  = np.ones([self.nts]) * 1. 	# Building activity vector
		self.occupancy_vector	  = np.ones([self.nts]) * 1.	# Number of occupants in building per timestep
		self.household_vector	  = []							# Vector of households occupancy
		self.occupants			  = 0							# Number of occupants in building
		self.occupied			  = 1.							# Is building occupied?
		self.perc_occupied		  = 0.							# Percentage of occupied dwellings in building
		
		# Heat gains
		self._solar_gains	  	  = _solar_gains
		self._internal_gains  	  = _internal_gains
		
		# DSM
		self.refurbishment_level  = refurbishment_level
		self._night_set_back  	  = _night_set_back				# Share of buildings with nsb
		self.schedule_nsb		  = schedule_nsb				# [start, end] of nsb in h
		self.T_nsb		  		  = T_nsb						# Night set-back temperature in degC
		self.power_reduction 	  = power_reduction				# Percentage of power reduced (as decimal)
		
		# Results
		self.space_heating_power  = np.zeros([self.nts]) 		# Power required so that Tb >= Tset in W
		self.solar_gains  		  = np.zeros([self.nts]) 		# Power required so that Tb >= Tset in W
		self.internal_gains  	  = np.zeros([self.nts]) 		# Power required so that Tb >= Tset in W
		self.Tb  				  = np.zeros([self.nts]) 		# Building temp in degC
		self.space_heating_energy = 0.							# Aggregated heating demand in Wh
		self.space_heating_energy_per_area = 0.					# Aggregated heating demand in kWh/m2

		
		# Hot water demand
		# --------------------------------------------------
		# Domestic hot water consumption
		self.Tw					  = Tw							# Hot water temperature in [degC]
		self.daily_DHW			  = 0.							# Daily domestic hot water consumption in m3
		self.dhw_prob			  = dhw_prob					# Probabilities for dhw-loads
		self.hw_tank_capacity	  = 0.							# Hot water tank capacity in m3
		self.hw_tank_volume_t0	  = 0.							# Initial state of hot water tank in m3
		self.hw_tank_limit	   	  = hw_tank_limit 				# Hot water tank limit as perc (decimal)
		self.hw_flow		      = hw_flow				    	# Flow to refill hot water tank in L/min
		
		# Seasonality
		self.day_vector			  = day_vector					# Vector of days in simulation time frame
		self.seasonal_vector	  = seasonal_vector				# Sinusoidal function for seasonal variations of DHW consumption
		
		# Results
		self.hot_water_m3		  = np.zeros([self.nts])		# Hot water demand in m3 (instantaneous demand)
		self.hot_water_tank_m3	  = np.zeros([self.nts])		# Hot water demand in m3 (tank demand)
		self.hot_water_power      = np.zeros([self.nts])		# Hot water demand in W
		self.hot_water_energy	  = 0.							# Aggregated hot water demand in Wh
		self.dhw_debug	 		  = np.zeros([self.nts, 20])	# Hot water demand data per time step
		
		
		# Total heating energy demand
		# --------------------------------------------------
		# Instantaneous demand
		self.total_power   		  = np.zeros([self.nts])		# Total heat demand in W
		self.total_energy	  	  = 0.							# Aggregated total heat demand in Wh
		
		# Delayed demand
		self.min_vector			  = min_vector					# Vector of simulation time steps in minutes
		self.distance2hp		  = distance2hp					# Distance to heat plant in m
		self.delay				  = 0.							# Delay due to distance in min
		self.delayed_min_vector	  = np.zeros([self.nts])		# Vector of delayed time steps in min
		self.total_power_delayed  = np.zeros([self.nts])		# Delayed total heat demand in W
		self.total_energy_delayed = 0.							# Delayed total energy demand in Wh
		

		# Reporting
		# --------------------------------------------------
		self.plot				  = plot
		self.save				  = save
		self.debug				  = debug
				
# Heating energy demand
# --------------------------------------------------------------------------------		
	def calculate_space_heating_demand(self):
		"""
		Calculates building space heating demand as timeseries [W] and aggregated value [Wh]
		"""
		
		# self.U, self.C
		self.calculate_building_thermal_properties()
		
		# self. Tset
		self.calculate_building_Tset()
		
		# self.activity_vector, self.occupancy_vector
		self.calculate_building_activity_occupancy_vector()
		
		# Initial building temperature (Tb[0])
		if (self.Tb0_str   == 'Tamb'):
			self.Tb0        = self.Tamb[0]
		elif (self.Tb0_str == 'Tset'):
			self.Tb0        = self.Tset
		
		# Space heating demand
		# --------------------------------------------------
		self.my_space_heating_demand = SpaceHeatingDemand(self.dt_vector, self.resolution, self.heated_area,\
										self.Tamb, self.I, self.Tb0, self.dT_per_hour, self.eta, self.thermal_inertia, \
										self.U, self.V, self.C, self.Tset, self.dTset, \
										self.activity_vector, self.occupancy_vector, \
										self._solar_gains, self._internal_gains, \
										self._night_set_back, self.schedule_nsb, self.T_nsb, self.power_reduction,
										self.window_areas, self.coords, self.debug)
		self.my_space_heating_demand.calculate()
		
		# Power [W]
		self.space_heating_power	 = self.my_space_heating_demand.sh_power
		self.solar_gains			 = self.my_space_heating_demand.solar_gains
		self.internal_gains	 		 = self.my_space_heating_demand.internal_gains
		
		# Building temperature [degC]
		self.Tb	 					 = self.my_space_heating_demand.Tb
		
		# Energy [Wh]
		self.space_heating_energy    = (self.space_heating_power.sum() * self.resolution * 1 / 60)
		self.solar_gains_energy    	 = (self.solar_gains.sum() * self.resolution * 1 / 60)
		self.internal_gains_energy   = (self.internal_gains.sum() * self.resolution * 1 / 60)
		
		# Energy per (conditioned) area [kWh/m2]
		self.space_heating_energy_per_area = (self.space_heating_energy / self.heated_area) / 1e3
		self.solar_gains_per_area 		   = (self.solar_gains_energy / self.heated_area) / 1e3
		self.internal_gains_per_area 	   = (self.internal_gains_energy / self.heated_area) / 1e3
#
	def calculate_hot_water_demand(self, save_debug):
		"""
		Calculates building hot water demand as timeseries [W] and [m3] and aggregated value [Wh]
		Only for residential buildings
		
		args:
			save_debug	<boolean>	Is debug file saved?
		
		returns:
			self.hot_water_m3  
			self.hot_water_power
			self.hot_water_energy
		"""
		
		# Building areas
		self.calculate_building_envelope_areas()
		
		if (self.use[0] == 3):	# residential
		
			# self.activity_vector, self.occupancy_vector
			self.calculate_building_activity_occupancy_vector()
		
			# calculate daily nominal load
			self.calculate_daily_hot_water_demand()
						
			# define size and state of hot water tank
			self.define_hot_water_tank()
			
			# Hot water demand
			# --------------------------------------------------
			self.my_hot_water_demand = HotWaterDemand(self.dt_vector, self.resolution, \
											self.day_vector, self.seasonal_vector, self.activity_vector, \
											self.Tw, self.daily_DHW, self.dhw_prob, \
											self.hw_tank_capacity, self.hw_tank_limit, self.hw_tank_volume_t0, self.hw_flow, \
											self.result_dir, self.use, self.year_class, self.btype, self.bid, \
											self.debug, save_debug)
			self.my_hot_water_demand.calculate()
			
			if save_debug:
				self.dhw_debug		 = self.my_hot_water_demand.dhw_debug
			
			# Share of dhw per load
			self.dhw_load_share 	 = (self.my_hot_water_demand.agg_consumption / self.my_hot_water_demand.agg_consumption[4]) * 100
			
			# Power [W] 
			self.hot_water_power	 = self.my_hot_water_demand.dhw_power
			
			# Flow rate [m3]
			self.hot_water_m3		 = self.my_hot_water_demand.dhw_m3
			self.hot_water_tank_m3	 = self.my_hot_water_demand.dhw_tank_m3
			
			# Energy [Wh]
			self.hot_water_energy    = (self.hot_water_power.sum() * self.resolution * 1 / 60)
#
	def calculate_total_heat_demand(self):
		"""
		Agrregate the space heating and/or hot water demand time series.
		The total time series is delayed depending on the distance to the heat plant.
		"""

		# Space heating demand
		if self._space_heating:
			self.total_power  += self.space_heating_power
			self.total_energy += self.space_heating_energy
		
		# Hot water demand
		if self._hot_water:
			self.total_power  += self.hot_water_power
			self.total_energy += self.hot_water_energy
		
		# Delayed time series
		self.calculate_delayed_timeseries()
#
	def calculate_delayed_timeseries(self, flow_vel = 1.):
		"""
		Delays vector of heat demand depending on the distance of the building 
		centroid to the (geothermal) heat plant. A flow velocity of 1 m/s in the district heating network is 
		considered.
		"""
		
		# calculate flow velocity in minutes
		flow_vel_min     = flow_vel * 60. # m/min
		
		# calculate delay
		# Delay is two times the distance between the building and the power plant (round trip: flow
		# from the power plant to the building and back)
		self.delay 	   	 = 2. * self.distance2hp / flow_vel_min # min
		
		# add delay to datetime minute vector
		for iii, dt in enumerate(self.min_vector):
			self.delayed_min_vector[iii] = self.min_vector[iii] - self.delay
		
		# interpolate to calculate delayed time series in base datetime vector
		f 				 		  = scipy.interpolate.interp1d(self.delayed_min_vector, self.total_power, fill_value = "extrapolate")
		self.total_power_delayed  = f(self.min_vector)
		
		# calculate delayed energy demand
		self.total_energy_delayed = (self.total_power_delayed.sum() * self.resolution * 1 / 60)
					
# Building parameters
# --------------------------------------------------------------------------------
	def calculate_building_thermal_properties(self):
		"""
		Calculates equivalent U-value, thermal capacitance (C) and time constant (Tau)
		for the building. These properties are used in the first order thermal model.
		
		returns:
			self.U		<float>		Transmission losses (U-value) [W/K]
			self.V		<float>		Ventilation losses [W/K]
			self.C		<float>		Thermal capacitance [J/K]
			self.Tau	<float>		Time constant [s]
		"""
			
		# calculate building envelope areas in contact with ambient temperature
		self.calculate_building_envelope_areas()
		
		# assign refurbishment level from statistics (only in case it is not fixed)
		if np.mean(self.refurbishment_level) == 0:
			self.determine_building_refurbishment_level()
		
		# calculate building thermal properties per unit area
		u, v, c = self.get_building_thermal_properties_per_unit_area()
			
		# multiply u and c by envelope areas
		# consider adjustment factor for border situation acording to TABULA
		b_tr = 0.5	# adjustment factor for surface bordering on soil
		U_roof   = u[0] * self.env_areas[0]
		U_wall   = u[1] * self.env_areas[1]
		U_floor  = u[2] * self.gfa_h * b_tr
		U_window = u[3] * self.env_areas[3]
		
		C_roof   = c[0] * self.env_areas[0]
		C_wall   = c[1] * self.env_areas[1]
		C_floor  = c[2] * self.env_areas[2]
		
		# calculate ventilation losses, V [W/K]
		# air_specific_heat_capacity * air_density = 0.34 Wh/(m3 K)
		# air flow rate related to usage is considered as half the value in TABULA (0.4)
		self.V	 = (v[0] + v[1]) * self.heated_area * self.f2f_height * 0.34
		
		# Equivalent building thermal properties
		### Transmission losses (U-value), U [W/K]
		self.U 	 = U_roof + U_wall + U_floor + U_window
		"""
		print('          U_roof:   {}*{} = {}'.format(u[0],self.env_areas[0],U_roof))
		print('          U_wall:   {}*{} = {}'.format(u[1],self.env_areas[1],U_wall))
		print('          U_floor:  {}*{}*0.5 = {}'.format(u[2],self.gfa_h,U_floor))
		print('          U_window: {}*{} = {}'.format(u[3],self.env_areas[3],U_window))
		print('          U_eq:     {}'.format(self.U))
		raw_input()
		"""
		
		self.adjust_building_thermal_properties()
				
		### Thermal capacitance, C [J/K]
		self.C 	 = C_floor + C_wall + C_roof
		
		# Time constant, Tau [s]
		self.Tau = self.C / (self.U + self.V)		
#
	def adjust_building_thermal_properties(self):
		"""
		Empirical adjustment of U-values to match TABULA results
		"""

		if (np.mean(self.refurbishment_level) >= 1. and np.mean(self.refurbishment_level) < 2.):
			if   self.btype[0] == 0: # SFH
				self.U = self.U * 1.
			elif self.btype[0] == 1: # TH
				self.U = self.U * 1.2
			elif self.btype[0] == 2: # MFH
				self.U = self.U * 1.18
			elif self.btype[0] == 3: # AB
				self.U = self.U * 1.18
		elif (np.mean(self.refurbishment_level) >= 2. and np.mean(self.refurbishment_level) < 3.):
			if   self.btype[0] == 0: # SFH
				self.U = self.U * 1.35
			elif self.btype[0] == 1: # TH
				self.U = self.U * 1.31
			elif self.btype[0] == 2: # MFH
				self.U = self.U * 1.53
			elif self.btype[0] == 3: # AB
				self.U = self.U * 1.5
		elif (np.mean(self.refurbishment_level) >= 3.):
			if   self.btype[0] == 0: # SFH
				self.U = self.U * 1.6
			elif self.btype[0] == 1: # TH
				self.U = self.U * 1.65
			elif self.btype[0] == 2: # MFH
				self.U = self.U * 1.7
			elif self.btype[0] == 3: # AB
				self.U = self.U * 1.55
#
	def get_building_thermal_properties_per_unit_area(self):
		"""
		Gets building thermal properties from TABULA Web Tool data based on the
		building typology [year_class, btype]
		
		returns:
			u	<list>	[u_roof, u_wall, u_floor, u_window] in W/(K m2)
			v	<list>  [v_usage, v_infiltration] in 1/h
			c	<list>	[c_roof, c_wall, c_floor] in J/(K m2)
		"""	
			
		if (self.use[0] == 3): # RESIDENTIAL
			
			### ROOF
			if  self.refurbishment_level[0] == 1:
				u_roof   = self.U_RES[0][0][self.year_class[0]][self.btype[0]]
			elif self.refurbishment_level[0] == 2:
				u_roof   = self.U_RES[0][1][self.year_class[0]][self.btype[0]]
			elif self.refurbishment_level[0] == 3:
				u_roof   = self.U_RES[0][2][self.year_class[0]][self.btype[0]]
			
			### WALL
			if  self.refurbishment_level[1] == 1:
				u_wall   = self.U_RES[1][0][self.year_class[0]][self.btype[0]]
			elif self.refurbishment_level[1] == 2:	
				u_wall   = self.U_RES[1][1][self.year_class[0]][self.btype[0]]
			elif self.refurbishment_level[1] == 3:	
				u_wall   = self.U_RES[1][2][self.year_class[0]][self.btype[0]]
				
			### FLOOR
			if  self.refurbishment_level[2] == 1:
				u_floor  = self.U_RES[2][0][self.year_class[0]][self.btype[0]]
			elif self.refurbishment_level[2] == 2:
				u_floor  = self.U_RES[2][1][self.year_class[0]][self.btype[0]]	
			elif self.refurbishment_level[2] == 3:	
				u_floor  = self.U_RES[2][2][self.year_class[0]][self.btype[0]]	
			
			### WINDOW
			if  self.refurbishment_level[3] == 1:
				u_window = self.U_RES[3][0][self.year_class[0]][self.btype[0]]
			elif self.refurbishment_level[3] == 2:										
				u_window = self.U_RES[3][1][self.year_class[0]][self.btype[0]]
			elif self.refurbishment_level[3] == 3:											
				u_window = self.U_RES[3][2][self.year_class[0]][self.btype[0]]
				
			### Ventilation losses (air exchange rate) [1/h]
			overall_refurbishment_level = np.ceil(np.mean(self.refurbishment_level))
			if overall_refurbishment_level == 1:
				v_usage  = self.V_RES[0][0][self.year_class[0]][self.btype[0]]
				v_inf    = self.V_RES[1][0][self.year_class[0]][self.btype[0]]
			elif overall_refurbishment_level == 2:
				v_usage  = self.V_RES[0][1][self.year_class[0]][self.btype[0]]
				v_inf    = self.V_RES[1][1][self.year_class[0]][self.btype[0]]
			elif overall_refurbishment_level == 3:
				v_usage  = self.V_RES[0][2][self.year_class[0]][self.btype[0]]
				v_inf    = self.V_RES[1][2][self.year_class[0]][self.btype[0]]
			
			### Thermal capacitance [J/(K m2)]
			# >>> Consider refurbishment
			c_roof		 = self.C_RES[0][self.year_class[0]][self.btype[0]]
			c_wall		 = self.C_RES[1][self.year_class[0]][self.btype[0]]
			c_floor		 = self.C_RES[2][self.year_class[0]][self.btype[0]]
			
		else: # NON-RESIDENTIAL	
			
			### Transmission losses (U-values) [W/(K m2)]
			# >>> Consider refurbishment
			u_roof   = self.U_NRES[self.year_class[0]][0]
			u_wall   = self.U_NRES[self.year_class[0]][1]
			u_floor  = self.U_NRES[self.year_class[0]][2]
			u_window = self.U_NRES[self.year_class[0]][3]
			
			### Ventilation losses (air exchange rate) [1/h]
			# >>> Consider refurbishment
			v_usage   = self.V_NRES[self.year_class[0]][0]
			v_inf     = self.V_NRES[self.year_class[0]][1]
			
			### Thermal capacitance [J/(K m2)]
			# >>> Consider refurbishment
			c_roof   = self.C_NRES[self.year_class[0]][0]
			c_wall   = self.C_NRES[self.year_class[0]][1]
			c_floor  = self.C_NRES[self.year_class[0]][2]
			
		return [u_roof, u_wall, u_floor, u_window], [v_usage, v_inf], [c_roof, c_wall, c_floor]
#
	def calculate_building_envelope_areas(self):
		"""
		Calculates building envelope areas (wall, roof and window) depending on the 
		building use: 
		Residential: areas are calculated according to building typologies in TABULA.
					 To obtain the typology (construction year class and building type)
					 the gfa is matched to the typology GFA.
		Non-residential: number of floors and window-to-wall ratio are
					 derived from statistics and used to calculate the building areas.
					 Construction year class is also derived from statistics.
		
		An area_correction_factor is included since not the whole floor area is heated.
		From VDI 3807-2.
		
		
		returns:
			self.env_areas		<list>		[Roof, Walls, Floor (gfa), Window] in m2
			self.heated_area 	<float>		Total building heated area
		"""
			
		if (self.use[0] == 3): # RESIDENTIAL
		
			# Area correction factor for calculating heated areas
			area_correction_factor = 0.85
			
			# Conditioned (heated) GFA
			self.gfa_h = self.env_areas[2] * area_correction_factor
		
			# check if building typology is given
			if self.building_typology:			
				self.year_class           = UrbanHeatPro.year_class_to_tuple(self.use[0], self.year_class_int) # Construction year class (int, str)
				self.btype				  = UrbanHeatPro.building_type_to_tuple(self.use[0], self.btype_int)	    # Building type (int, str)
		
			else:
				# assign building typology from statistics	
				# check that the building data exists as C-values are not available for all buildings
				building_exists = False
				while building_exists == False:
					self.categorize_residential_building()
					a = self.C_RES[2][self.year_class[0]][self.btype[0]] # C_floor value
					b = self.GFA[self.year_class[0]][self.btype[0]]
					if (a > 0.0 and b > 0.0):
						building_exists = True
			
			# Floor-to-floor height
			# from http://www.ctbuh.org/HighRiseInfo/TallestDatabase/Criteria/HeightCalculator/tabid/1007/language/en-GB/Default.aspx
			# self.f2f_height = (randint(30, 32) / 10.0) # m
			# from TABULA
			self.f2f_height = 2.5 # m
			
			# Number of floors
			self.floors 	= self.FLOORS[self.year_class[0]][self.btype[0]] * (1 + (1 - area_correction_factor))

			# Building areas
			### Roof
			self.env_areas[0] 	= self.ARATIO[0][self.year_class[0]][self.btype[0]] * self.gfa_h
			### Wall
			self.env_areas[1]	= self.ARATIO[1][self.year_class[0]][self.btype[0]] * self.gfa_h * self.free_walls					
			### Window
			self.env_areas[3] 	= self.ARATIO[2][self.year_class[0]][self.btype[0]] * self.gfa_h
						
			for j in range(4):
				self.window_areas[j] = self.env_areas[3] * self.WRATIO_ORIENTATION[5 * j + self.btype[0]]

			
		else: # NON-RESIDENTIAL	
		
			# Area correction factor for calculating heated areas
			area_correction_factor = (randint(6, 8) / 10.0)
			
			# Conditioned (heated) GFA
			self.gfa_h = self.env_areas[2] * area_correction_factor
		
			# Construction year class
			self.categorize_non_residential_building()
			
			# Floor-to-floor height
			# from http://www.ctbuh.org/HighRiseInfo/TallestDatabase/Criteria/HeightCalculator/tabid/1007/language/en-GB/Default.aspx
			self.f2f_height = (randint(30, 39) / 10.0) # m
			
			# Number of floors
			# >>> Source missing
			lower_limit     = 1
			upper_limit     = 3
			self.calculate_number_of_floors(lower_limit, upper_limit)
			
			# Building areas
			### Roof
			self.env_areas[0]   = self.gfa_h
			### Wall
			width 				= np.sqrt(self.gfa_h) # [m] building is assumed to be a cube
			self.env_areas[1]   = self.free_walls * self.floors * self.f2f_height * width
			### Window
			# >>> Add Window-to-Wall ratio for the different building types
			# >>> Source missing
			self.env_areas[3]   = (randint(1, 4) / 10.0) * self.env_areas[1]

			# Window areas oriented to [east, south, west, north]
			for j in range(4):
				self.window_areas[j] = self.env_areas[3] * self.WRATIO_ORIENTATION[5 * j + 4]			
			
		# Total heated area (reference area in TABULA) in m2
		self.heated_area = self.gfa_h * self.floors
		
		# Total living space
		self.living_area = self.gfa_h * self.floors * (1 + (1 - area_correction_factor))
#
	def determine_building_refurbishment_level(self):
		"""
		Determines the refurbishment level for the different building components depending
		on building type and construction year. Based on national statistics.
		Refurbishment levels:
			1	National minimum requirement
			2	Improved standard
			3	Ambitious standard
		Refurbishment level in non residential buildings is assumed to be 1
		>>> Statistics on refurbishment in non-residential buildings
		
		returns:
			self.refurbishment_level	<list>		[roof, wall, floor, window]
		"""
		
		self.refurbishment_level = [1, 1, 1, 1]
		
		# RESIDENTIAL
		if (self.use[0] == 3):
			# check building type
			if (self.btype[0] == 0 or self.btype[0] == 1): # SFH or TH
				# check building construction year class
				if self.year_class[0] <= 5: # < 1978
					perc_refurbished = self.REFURBISHED_RES[0]
				elif self.year_class[0] > 5 and self.year_class[0] <= 7: # 1978 - 1994
					perc_refurbished = self.REFURBISHED_RES[1]
				elif self.year_class[0] > 7 and self.year_class[0] <= 9: # > 1994
					perc_refurbished = self.REFURBISHED_RES[2]
			else: # MFH or AB
				# check building construction year class
				if self.year_class[0] <= 5: # < 1978
					perc_refurbished = self.REFURBISHED_RES[3]
				elif self.year_class[0] > 5 and self.year_class[0] <= 7: # 1978 - 1994
					perc_refurbished = self.REFURBISHED_RES[4]
				elif self.year_class[0] > 7 and self.year_class[0] <= 9: # > 1994
					perc_refurbished = self.REFURBISHED_RES[5]
			
			# assign refurbishment level for every building component
			# if random number is smaller than percentage, then refurfishment level of
			# component is 2, otherwise is 1
			for component in range(4): # [roof, wall, floor, window]	
				rand_num = uniform(0, 1)
				if rand_num <= perc_refurbished[component]:
					self.refurbishment_level[component] = 2
		
		# NON-RESIDENTIAL
		else: 
			# all components are assumed to have refurbishment level 1.
			# No statistics
			self.refurbishment_level = [1, 1, 1, 1]
#	
	def categorize_residential_building(self):
		"""
		Categorizes building according to TABULA typologies. 
		Construction year class and building type are calculated by comparing the residential
		building ground floor area (gfa) with the GFA of typical buildings (from TABULA).
		Values are adapted to building stock statistics.
		
		returns:
			self.class_year	  <tuple>	(int, str)
			self.btype		  <tuple>	(int, str)
		"""
		
		gfa = self.env_areas[2]
		
		# get size of GFA matrix
		rows 			= len(self.GFA)		# construction year class
		cols 			= len(self.GFA[0])	# building type
		
		# initialize distance vectors
		distance        = np.zeros(rows * cols)
		distance_inv    = np.zeros(rows * cols)
		row_col 	    = [[] for _ in range(rows * cols)]
		t_distance_inv  = 0
		kkk 		    = 0
		
		# calculate inverse of distances
		# 1/d is used so that the minimum distance has a higher probability
		for iii in range(0, rows):
			for jjj in range(0, cols):
				distance[kkk]     = abs(self.GFA[iii][jjj] - gfa) + 0.01 # "+ 0.001" to avoid having x/0
				distance_inv[kkk] = (1.0 / distance[kkk]) * self.STOCK_RES[iii, jjj]
				row_col[kkk] 	  = [iii, jjj]
				t_distance_inv   += distance_inv[kkk] 
				kkk              += 1
		
		# normalize inverse distance
		norm_distance   = distance_inv / t_distance_inv

		# get cumulative density function
		cdf             = np.cumsum(norm_distance)

		# generate random number
		rnd       	    = np.random.uniform(0, 1, 1)

		# check random number with cdf
		index 			= np.argmax(rnd < cdf)
		
		# Construction year class
		self.year_class = UrbanHeatPro.year_class_to_tuple(self.use[0], row_col[index][0])
		
		# Building type
		self.btype 		= UrbanHeatPro.building_type_to_tuple(self.use[0], row_col[index][1])
#	
	def categorize_non_residential_building(self):
		"""
		int	construction year class
		0	< 1918
		1	1919 - 1976
		2	1977 - 1983
		3	1984 - 1994
		4	> 1995

		>>> Statistics on the non-residential buildings stock are missing
		"""
		
		# assign random construction year class
		year_class_int = randint(0, 4)
		
		self.year_class = UrbanHeatPro.year_class_to_tuple(self.use[0], year_class_int)
		self.btype		= UrbanHeatPro.building_type_to_tuple(self.use[0], self.btype_int)
#
	def calculate_number_of_floors(self, lower_limit, upper_limit):
		"""
		Calculates number of floors as random number between lower_limit and upper_limit.
		
		returns:
			self.floors		<int>
		"""
		
		self.floors = randint(lower_limit, upper_limit)
#
	def calculate_building_Tset(self):
		"""
		Derives a target temperature by choosing a random temperature from Tset_mean 
		+/- dT. Values differ for different building types.
		From http://tc76.org/spc100/docs/IBP%2018599/18599-10.pdf
				
		returns:
			self.Tset	<float>		Building target temperature in degC
		"""
		
		if   (self.use[0] == 0): # commercial
			Tset_mean  = self.TSET[0][0]
			dT		   = self.TSET[1][0]
		elif (self.use[0] == 1): # industrial
			Tset_mean  = self.TSET[0][1]
			dT		   = self.TSET[1][1]
		elif (self.use[0] == 2): # public
			Tset_mean  = self.TSET[0][2]
			dT		   = self.TSET[1][2]
		elif (self.use[0] == 3): # residential
			Tset_mean  = self.TSET[0][3]
			dT		   = self.TSET[1][3]
			
		self.Tset = float(randint(Tset_mean - dT, Tset_mean + dT))
#
	def calculate_daily_hot_water_demand(self):
		"""
		Returns the daily hot water demand by getting a random value 
		from the cdf function based on the statistics from VDI 3807-3 
		(specific dhw demand in m3/m2 of living area)
		
		returns:
			self.daily_DHW		[m3/day]
		"""
		
		# calculate specific dhw demand
		rand_num 	 = np.random.uniform(0, 1, 1)
		specific_dhw = self.dhw_prob[0](rand_num)[0]
		
		# calculate building dhw demand
		self.daily_DHW = specific_dhw * self.living_area 
#
	def define_hot_water_tank(self):
		"""
		Calculates size and initial state of hot water tank
		
		returns:
			self.hw_tank_capacity
			self.hw_tank_volume_t0
		"""
		
		# set tank capacity as function of daily DHW limit
		self.hw_tank_capacity = 1.5 * self.daily_DHW
		
		# define initial state of hot water tank as random percentage of capacity
		self.hw_tank_volume_t0 = np.random.uniform(0, 1, 1)[0] * self.hw_tank_capacity		
#
	def calculate_building_activity_occupancy_vector(self):
		"""
		Calculate vector of activity in building, i.e. percentage of occupied dwellings (for space heating)
		Active_hours (scheduled), building occupancy and weekends are considered
				
		returns:
			self.activity_vector	<numpy array> 
		"""
		
		# calculate number of occupants in building
		self.calculate_number_of_occupants()
		
		# calculate building active hours
		self.calculate_building_active_hours()
		
		# calculate occupant schedule
		if self._active_population:
			self.calculate_occupants_schedule()
		
		# calculate activity vector
		for iii in range(0, self.nts):
			hour 	= self.dt_vector[iii].hour
			wd		= self.dt_vector[iii].weekday() # 0 Monday ... 6 Sunday
			
			if   (self.use[0] == 0): # commercial
				if wd == 6: 	 	 # no activities on Sunday
					self.activity_vector[iii] = 0.
				else:
					if (hour < self.active_hours[0][0] or hour > self.active_hours[0][1]):
						self.activity_vector[iii] = 0.	
				
				# calculate occupancy vector
				self.occupancy_vector[iii] = self.activity_vector[iii] * self.occupants
			
			elif (self.use[0] == 1 or self.use[0] == 2): # industrial/public
				if (wd >= 5): 		 # no activities on Saturday/Sunday
					self.activity_vector[iii] = 0.
				else:
					if (hour < self.active_hours[0][0] or hour > self.active_hours[0][1]):
						self.activity_vector[iii] = 0.
				
				# calculate occupancy vector
				self.occupancy_vector[iii] = self.activity_vector[iii] * self.occupants
				
			elif (self.use[0] == 3): # residential
				# calculate occupancy vector
				self.occupancy_vector[iii] = self.activity_vector[iii] * self.occupants
				
				# recalculate building occupancy considering active population if the flag active_pop is True
				if self._active_population:
					# initialize counters
					unoccupied_dwellings      = 0
					occupants_not_in_building = 0
					if (wd < 5): 		 # weekends are considered as active
						# occupant loop
						for dwelling in range(self.dwellings):
							occupants_not_at_home  = 0
							for occupant in self.household_vector[dwelling]:	
								# check if occupant is at home
								if (hour > occupant[1][0][1] and hour < occupant[1][1][0]): # not at home								
									occupants_not_at_home     += 1
									occupants_not_in_building += 1
							# check if building is empty
							empty = (occupant[0] == occupants_not_at_home)
							if empty:
								unoccupied_dwellings += 1.
					self.occupancy_vector[iii] = self.occupants - occupants_not_in_building
					self.activity_vector[iii]  = (self.dwellings - unoccupied_dwellings) / self.dwellings
#
	def calculate_number_of_occupants(self):
		"""
		Calculates random number of occupants based on:
			- Household size and number of dwellings for residential buildings
			- Use area for non-residential buildings
				From: https://www.engineeringtoolbox.com/number-persons-buildings-d_118.html
		
		returns:
			self.occupants	<int>
		"""
		
		if   (self.use[0] == 0): # commercial
			area_per_person = randint(2, 7) 	# in m2
			# retail, supermarket, department stores
			self.occupants = self.heated_area * area_per_person
		elif (self.use[0] == 1): # industrial
			# light manufacturing, heavy manufacturing
			area_per_person = randint(10, 30)	# in m2
			self.occupants = self.heated_area * area_per_person
		elif (self.use[0] == 2): # public
			# municipal buildings, library, museum
			area_per_person = randint(2, 10)	# in m2
			self.occupants = self.heated_area * area_per_person
		elif (self.use[0] == 3): # residential
			self.calculate_number_of_dwellings()
			self.occupants = self.calculate_household_vector()
#		
	def calculate_building_active_hours(self):
		"""
		Assigns random start and end hours for building active hours.
		Values differ for different building types.
				
		returns:
			self.active_hours	<list>		[(start0, end0), (start1, end1)] in h
		"""
		
		### Residential
		if (self.use[0] == 3): 		
			self.active_hours = [[0, 23], [0, 0]]
	
		### Non-residential
		else:
			if   (self.use[0] == 0): # commercial
				start_mean = self.SCHEDULE[0][0][0]
				end_mean   = self.SCHEDULE[0][0][1]
				dt		   = self.SCHEDULE[1][0]
			elif (self.use[0] == 1): # industrial
				start_mean = self.SCHEDULE[0][1][0]
				end_mean   = self.SCHEDULE[0][1][1]
				dt		   = self.SCHEDULE[1][1]
			elif (self.use[0] == 2): # public
				start_mean = self.SCHEDULE[0][2][0]
				end_mean   = self.SCHEDULE[0][2][1]
				dt		   = self.SCHEDULE[1][2]
			
			start = randint(start_mean - dt, start_mean + dt)
			end   = randint(end_mean - dt, end_mean + dt)
			self.active_hours = [[start, end], [0, 0]]	
#
	def calculate_number_of_dwellings(self):
		"""
		Calculates number of dwellings based on the building living area and mean dwelling size
		SFH and TH have only 1 or 2 dwellings.
		
		returns:
			self.dwellings	<int>
		"""
		
		mean_dwelling_size = 82.7 # m2
		# from https://www.statistik.bayern.de/statistikkommunal/09184148.pdf
			
		### For SFH and TH
		# compare random number to single-dwelling buildings statistics
		rand_num = uniform(0, 1)
		if (self.btype[0]) == 0: # SFH
			perc_sdwelling = self.SINGLE_DWELLING[self.year_class[0]][self.btype[0]]
			if rand_num <= perc_sdwelling:
				self.dwellings = 1
			else:
				self.dwellings = 2
		elif (self.btype[0]) == 1: # TH
			perc_sdwelling = self.SINGLE_DWELLING[self.year_class[0]][self.btype[0]]
			if rand_num <= perc_sdwelling:
				self.dwellings = 1
			else:
				self.dwellings = 2
			
		### For MFH and AB
		# divide total living area by average dwelling size
		else:
			self.dwellings = int(np.ceil(self.living_area / mean_dwelling_size))
			
		# calculate dwelling size
		dwelling_size = self.living_area / self.dwellings
		self.determine_dwelling_size_category(dwelling_size)
#	
	def calculate_household_vector(self):
		"""
		Calculates vector of number of persons living in the dwellings of building
		based on the size in m2 of the dwelling.
		A schedule is assigned to every occupant based on studying/working schedule.
		
		returns:
			household_vector	<list>		[dwelling, [occupant, [schedule]]] for occupant and dwelling in building
		"""
		
		# calculate occupancy based on percentage of active population
		perc_working_pop   = 0.583 # percentage of population that works or studies
		# from https://ergebnisse.zensus2011.de/#StaticContent:091840148148,BEG_1_6_1,m,table
		
		# initialize vector with household size per building
		self.household_vector = [[] for dwelling in range(self.dwellings)]
		
		# calculate number of occupants
		occupant_counter = 0
		for dwelling in range(self.dwellings):		
			
			# Household size
			# x-values
			household_size   = range(1, 7)
			# CFD
			cfd = np.cumsum(self.HOUSEHOLD_SIZE[self.dwelling_size_cat])
			# get random_number and obtain x-value from cfd
			rand_num = np.random.uniform(0, 1, 1)
			index = np.argmax(rand_num < cfd)
			# Household size (number of occupants)
			occupants = household_size[index]
			
			# initialize vector with occupant schedule per dwelling per building
			#self.household_vector = [[[[],[]] for o in range(occupants)] for d in range(self.dwellings)]
			self.household_vector[dwelling] = [[[],[]] for occupant in range(occupants)]
			
			# add occupants to building
			occupant_counter += occupants
					
		return occupant_counter
#
	def calculate_occupants_schedule(self):
		"""
		A schedule is assigned to every occupant based on studying/working schedule.
		
		returns:
			household_vector	<list>		[dwelling, [occupant, [schedule]]] for occupant and dwelling in building
		"""
		
		# calculate occupancy based on percentage of active population
		perc_working_pop   = 0.583 # percentage of population that works or studies
		# from https://ergebnisse.zensus2011.de/#StaticContent:091840148148,BEG_1_6_1,m,table
		
		for dwelling in range(self.dwellings):
						
			# derive schedule for each occupant
			for iii, occupant in enumerate(self.household_vector[dwelling]):
				
				# add occupants in household vector
				occupant[0] = iii + 1
				
				# compare random number with percentage of working/student population
				# to determine if occupants are at home
				rand_num = np.random.uniform(0, 1, 1)
				if rand_num > perc_working_pop: # occupant leaves the building to work/study
					# calculate working/studying start/end time from "public" building usage
					start_mean = self.SCHEDULE[0][2][0]
					end_mean   = self.SCHEDULE[0][2][1]
					dt		   = self.SCHEDULE[1][2]
			
					start = randint(start_mean - dt, start_mean + dt)
					end   = randint(end_mean - dt, end_mean + dt)

					# define active_hours per occupant
					occupant[1] = [[0, start], [end, 23]]
				else:
					occupant[1] = [[0, 23], [0, 0]]
#
	def determine_dwelling_size_category(self, dwelling_size):
		"""
		Determine dwelling size category based on statistics 
		https://ergebnisse.zensus2011.de/#StaticContent:091840148148,GWZ_4_3_2,m,table
		"""
		
		if dwelling_size > 0 and dwelling_size <= 40:
			self.dwelling_size_cat = 0
		elif dwelling_size > 40 and dwelling_size <= 60:
			self.dwelling_size_cat = 1
		elif dwelling_size > 60 and dwelling_size <= 80:
			self.dwelling_size_cat = 2
		elif dwelling_size > 80 and dwelling_size <= 100:
			self.dwelling_size_cat = 3
		elif dwelling_size > 100 and dwelling_size <= 120:
			self.dwelling_size_cat = 4
		elif dwelling_size > 120 and dwelling_size <= 140:
			self.dwelling_size_cat = 5
		elif dwelling_size > 140 and dwelling_size <= 160:
			self.dwelling_size_cat = 6
		elif dwelling_size > 160 and dwelling_size <= 180:
			self.dwelling_size_cat = 7
		elif dwelling_size > 180 and dwelling_size <= 200:
			self.dwelling_size_cat = 8
		else:
			self.dwelling_size_cat = 9

# Results
# --------------------------------------------------------------------------------
	def plot_timeseries(self, \
				space_heating = True, Tb = False, \
				hot_water = True, \
				total = True):
		"""
		Plots heat demand timeseries
		"""
		
		if space_heating:
			fig_name = '{}/SpaceHeatingDemand_{}_{}_{}_{}.png'.format(self.result_dir, \
							self.use[0], self.year_class[0], self.btype[0], self.bid)
			
			UrbanHeatPro.plot_timeseries(self.dt_vector, \
					[self.space_heating_power], ['Space heating demand'], \
					fig_name, xticks = ('month', 3), \
					ynumticks = 'auto', ylabel = 'Power [kW]', ylim0 = True, yfactor = 1e3)
				
		if Tb:
			fig_name = '{}/BuildingTemperature_{}_{}_{}_{}.png'.format(self.result_dir, \
							self.use[0], self.year_class[0], self.btype[0], self.bid)
			
			UrbanHeatPro.plot_timeseries(self.dt_vector, \
			[self.Tamb, self.Tb], ['T_amb', 'T_b'], \
					fig_name, xticks = ('month', 3), \
					ynumticks = 'auto', ylabel = 'Temperature [degC]', ylim0 = False, yfactor = 1)
									
		if hot_water:
			fig_name = '{}/HotWaterDemand_{}_{}_{}_{}.png'.format(self.result_dir, \
							self.use[0], self.year_class[0], self.btype[0], self.bid)
			
			UrbanHeatPro.plot_timeseries(self.dt_vector, \
					[self.hot_water_power], ['Hot water demand'], \
					fig_name, xticks = ('month', 3), \
					ynumticks = 'auto', ylabel = 'Power [kW]', ylim0 = True, yfactor = 1e3)

		if total:
			fig_name = '{}/TotalHeatDemand_{}_{}_{}_{}.png'.format(self.result_dir, \
							self.use[0], self.year_class[0], self.btype[0], self.bid)
			
			UrbanHeatPro.plot_stacked_timeseries(self.dt_vector, \
					[self.hot_water_power, self.space_heating_power], \
					['Hot water', 'Space heating'], \
					fig_name, xticks = ('month', 3), \
					ynumticks = 'auto', ylabel = 'Power [kW]', ylim0 = True, yfactor = 1e3)	
#
	def save_csv(self):
		"""
		Saves key building parameters and heat demand (space heating, hot water and 
		total) as timeseries.
		"""
		
		array_to_save = np.array([self.dt_vector_excel, self.activity_vector, 
						 self.Tamb, 
						 self.Tb, 
						 self.occupancy_vector,
						 self.space_heating_power, self.internal_gains, self.solar_gains,
						 self.hot_water_power, self.hot_water_tank_m3,
						 self.total_power,
						 self.total_power_delayed]).transpose()
			
		filename = '{}/HeatDemand_{}_{}_{}_{}.csv'.format(self.result_dir, \
						self.use[0], self.year_class[0], self.btype[0], self.bid)
					
		# Building parameters
		with open(filename, 'w') as text_file:

			text_file.write('Building id;{};\n'.format(self.bid))
			text_file.write('Building use;{};{};\n'.format(self.use[0], self.use[1]))
			text_file.write('Year class;{};{};\n'.format(self.year_class[0], self.year_class[1]))
			text_file.write('Building type;{};{};\n'.format(self.btype[0], self.btype[1]))
			text_file.write('GFA [m2];{:.3f};\n'.format(self.env_areas[2]))
			text_file.write('Living area [m2];{:.3f};\n'.format(self.living_area))
			text_file.write('Heated area [m2];{:.3f};\n'.format(self.heated_area))
			text_file.write('No. floors;{:.3f};\n'.format(self.floors))
			text_file.write('No. free walls;{};\n'.format(self.free_walls))
			text_file.write('Refurbishment level;{};\n'.format(self.refurbishment_level))
			text_file.write('Tset [degC];{};\n'.format(self.Tset))
			text_file.write('Tb0 [degC];{};\n'.format(self.Tb[0]))
			text_file.write('Active hours [h];{};{};{};{};\n'.format(self.active_hours[0][0], self.active_hours[0][1], self.active_hours[1][0], self.active_hours[1][1]))
			text_file.write('Occupied?;{};\n'.format(self.occupied))
			text_file.write('U [kW/K];{:.5f};\n'.format(self.U / 1e3))
			text_file.write('V [kW/K];{:.5f};\n'.format(self.V / 1e3))
			text_file.write('C [MJ/K];{:.5f};\n'.format(self.C / 1e6))
			text_file.write('Tau [h];{:.3f};\n'.format(self.Tau / 3600.0))
			text_file.write('Space heating demand [kWh];{:.3f};\n'.format(self.space_heating_energy / 1e3))
			text_file.write('Space heating demand per liv_area [kWh/m2];{:.3f};\n'.format((self.space_heating_energy / self.living_area) / 1e3))
			text_file.write('Daily DHW demand [m3/day];{:.3f};\n'.format(self.daily_DHW))
			text_file.write('Hot water demand [m3];{:.3f};\n'.format(sum(self.hot_water_m3)))
			text_file.write('Hot water demand per liv_area [m3/m2];{:.3f};\n'.format(sum(self.hot_water_m3) / self.living_area))
			text_file.write('Hot water demand [kWh];{:.3f};\n'.format(self.hot_water_energy / 1e3))
			text_file.write('Hot water demand per liv_area [kWh/m2];{:.3f};\n'.format((self.hot_water_energy / self.living_area) / 1e3))
			text_file.write('Total heat demand [kWh];{:.3f};\n'.format(self.total_energy / 1e3))
			text_file.write('Total heat demand per liv_area [kWh/m2];{:.3f};\n'.format((self.total_energy / self.living_area) / 1e3))
			text_file.write('datenum;active?;Tamb_degC;Tb_degC;Occupants;SpaceHeatingD_W;int_gains_W;sol_gains_W;HotWaterD_W;HotWaterTank_m3;TotalD_W;DelayedTotalD_W\n')
		
		# Time series
		with open(filename, 'a') as f:
			np.savetxt(f, array_to_save, delimiter = ';', fmt='%.4f')
#
	def save_load_duration_curve(self):
		"""
		Save sorted demand
		"""
	
		filename = '{}/LoadDurationCurve_{}_{}_{}_{}.csv'.format(self.result_dir, \
						self.use[0], self.year_class[0], self.btype[0], self.bid)
		
		# sort demand to save
		array_to_save = np.sort(self.total_power_delayed)[::-1].transpose()
		
		with open(filename, 'a') as f:
			np.savetxt(f, array_to_save, delimiter = ';', fmt='%.4f')
#
	def save_dhw_debug_csv(self):
		"""
		Saves debug values for dhw demand.
		"""
		
		array_to_save = np.array(self.dhw_debug)
		
		filename = '{}/DHWDemand_{}_{}_{}_{}.csv'.format(self.result_dir, \
						self.use[0], self.year_class[0], self.btype[0], self.bid)
					
		# Building parameters
		with open(filename, 'w') as text_file:
			text_file.write('timstep;hour;day;day_in_year;daily_dhw;daily_limit;activity;shower?;bath?;medium?;small?;fr_shower;fr_bath;fr_medium;fr_small;d_shower;d_bath;d_medium;d_small;c_shower;c_bath;c_medium;c_small;daily_agg;daily_convergence\n')
		
		# Time series
		with open(filename, 'a') as f:
			np.savetxt(f, array_to_save, delimiter = ';', fmt='%.4f')