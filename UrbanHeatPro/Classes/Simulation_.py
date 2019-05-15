"""
Simulation.py
AMC @ TUM ENS
"""

import os
import numpy as np
import pandas as pd
from datetime import datetime
from random import randint

import UrbanHeatPro.Functions as UrbanHeatPro
from .City import City

class Simulation():
# --------------------------------------------------------------------------------	
	def __init__(self, NAME, SIMULATION, CITY, SPACE_HEATING, HOT_WATER, REPORTING):
		
		# SIMULATION
		self.name 				  = NAME					# Simulation name
		self.region				  = SIMULATION[0][0]		# Name of region/city/urban area
		self.N 	  				  = SIMULATION[1][0]		# Number of runs
		self.resolution			  = SIMULATION[1][1]		# Temporal resolution in min
		self.timesteps			  = SIMULATION[1][2]		# Simulation time steps as list
		self.nts          		  = len(self.timesteps)		# Number of time steps
		self.dt_vector    		  = None		 			# Vector of time steps as datetime objects
		self.dt_vector_excel	  = None					# Vector of time steps as excel date
		self.my_dir    			  = os.getcwd()				# Current directory
		self.processes			  = SIMULATION[2][0]		# Number of parallel processes
		
		# CITY
		### External factor
		self.Tamb				  = None					# Vector of ambient temperature in degC
		self.I				  	  = None					# Solar radiation vector in W/m2 [I_Gh, I_Dh, I_ex, hs]
		### Buildings
		self.buildings_filename	  = CITY[0][0]				# Name of shapefile containing city data
		self.building_typology	  = CITY[0][1]				# Is typology given?
		self.building_stock_stats = None					# Data from building stock statistics
		self.buildings			  = None					# GeoDataFrame with buildings
		self.connection_factor	  = CITY[1][0]				# Share of buildings connected to the network
		### Flags
		self._space_heating	      = CITY[2][0]				# calculate space heating demand?
		self._hot_water		      = CITY[2][1]				# calculate hot water demand?
		### Base load
		self.base_load			  = CITY[3][0]				# Base load in W (min load)
		
		# SPACE HEATING DEMAND
		self.SPACE_HEATING		  = SPACE_HEATING			

		# HOT WATER DEMAND
		self.HOT_WATER			  = HOT_WATER
		
		# RESULTS	
		### Directory
		self.result_dir			  = None					# Directory to save results
		
		### Space heating demand 
		self.space_heating_power  = None					# City Space heating demand in W
		self.space_heating_energy = np.zeros([self.N])		# Aggregated city heating demand in Wh
		
		### Hot water demand
		self.hot_water_power      = None					# City hot water demand in W
		self.hot_water_energy	  = np.zeros([self.N])		# Aggregated city hot water demand in Wh
		
		### Total heat demand
		self.total_power   		  = None					# City total heat demand in W
		self.total_energy	  	  = np.zeros([self.N])		# Aggregated city total heat demand in Wh
		self.total_power_delayed  = None					# City total delayed heat demand in W
		self.total_energy_delayed = np.zeros([self.N])		# Aggregated city total delayed heat demand in Wh
		
		### Buildings per typology
		self.bstock_res   		  = np.array([np.zeros([10, 4]) for run in range(self.N)])
		self.bstock_nres  		  = np.array([np.zeros([5]) for run in range(self.N)])
		
		### Thermal properties per building typology
		self.u_res   			  = [[[np.array([]) for col in range(4)] for row in range(10)] for run in range(self.N)]
		self.c_res   			  = [[[np.array([]) for col in range(4)] for row in range(10)] for run in range(self.N)]
		self.tau_res   			  = [[[np.array([]) for col in range(4)] for row in range(10)] for run in range(self.N)]
		self.u_nres  			  = [[np.array([]) for row in range(5)] for run in range(self.N)]
		self.c_nres  			  = [[np.array([]) for row in range(5)] for run in range(self.N)]
		self.tau_nres  			  = [[np.array([]) for row in range(5)] for run in range(self.N)]
		
		### Heat demand per unit area per building typology
		self.sh_per_area_res   	  = [[[np.array([]) for col in range(4)] for row in range(10)] for run in range(self.N)]
		self.sh_per_area_nres     = [[np.array([]) for row in range(5)] for run in range(self.N)]
		self.hw_per_area_res   	  = [[[np.array([]) for col in range(4)] for row in range(10)] for run in range(self.N)]
		self.hw_per_area_nres     = [[np.array([]) for row in range(5)] for run in range(self.N)]
		self.heat_per_area_res    = [[[np.array([]) for col in range(4)] for row in range(10)] for run in range(self.N)]
		self.heat_per_area_nres   = [[np.array([]) for row in range(5)] for run in range(self.N)]
	
		# REPORTING
		self.REPORTING			  = REPORTING
		self.plot				  = REPORTING[0]					# Plot level  [0, 1, 2]
		self.save 				  = REPORTING[1]					# Save level  [0, 1, 2]
		self.debug		   		  = REPORTING[2]					# Debug level [0, 1, 2]		
#	
	def run(self, analyze_building_stock = False, analyze_thermal_properties = False):
		"""
		Runs a complete simulation of N runs of the city heat demand.
		"""
	
		print('\n---------------------------------------------------')
		print('UrbanHeatPRO')
		print('Calculation of urban heating energy demand profiles')
		print('AMC @ TUM ENS')
		print('---------------------------------------------------\n')
		
		# Input data
		if self.debug  != 0:
			print('\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
			print('INPUT DATA')
			print('---------------------------------------------------')
		
		### Building data 			-> self.buildings
		self.read_building_data()
		
		### Statistical data 		-> self.building_stock_stats
		self.read_input_data_csv()

		### Simulation time steps 	-> self.dt_vector
		self.calculate_dt_vector()
		self.nts = len(self.dt_vector)
		
		### Weather data			-> self.Tamb, self.I
		self.read_Tamb()
		self.read_I()

		### Result directory		-> self.result_dir
		self.prepare_result_directory()
		
		# initialize results matrices
		### Space heating demand
		self.space_heating_power  = np.zeros([self.nts, self.N]) # City Space heating demand in W
		### Hot water demand
		self.hot_water_power      = np.zeros([self.nts, self.N]) # City hot water demand in W
		### Total heat demand
		self.total_power   		  = np.zeros([self.nts, self.N]) # City total heat demand in W
		### Total heat demand (delayed)
		self.total_power_delayed  = np.zeros([self.nts, self.N]) # City total heat demand in W
		
		# Heat demand calculations
		for run in range(self.N):
		
			# debug
			if not self.debug == 0:
				print('\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
				print('\nRun {}/{}'.format(run, self.N - 1))
				print('___________________________________________________')
				
			# Result directory
			result_dir_run = '{}/{}'.format(self.result_dir, run)
			if not os.path.exists(result_dir_run): 
				os.makedirs(result_dir_run)
				
			# create City object
			self.create_city_object(run, result_dir_run)
			
			# create synthetic city
			if self.debug  != 0:
				print('\nSYNTHETIC CITY')
				print('---------------------------------------------------')
			self.my_city.create_synthetic_city()
			
			# ========================
		
			# calculate heating energy demand
			if self.debug  != 0:
				print('HEAT DEMAND')
				print('---------------------------------------------------')
			self.my_city.calculate_city_heat_demand()			
			
			# add run results to simulation results
			if self._space_heating:
				self.space_heating_power[:, run] = self.my_city.results['power_sh']
				self.space_heating_energy[run]   = self.my_city.results['energy_sh']
			if self._hot_water:
				self.hot_water_power[:, run] 	 = self.my_city.results['power_hw']
				self.hot_water_energy[run] 	 	 = self.my_city.results['energy_hw']
			self.total_power[:, run] 		 	 = self.my_city.results['power_total']
			self.total_energy[run] 			 	 = self.my_city.results['energy_total']
			self.total_power_delayed[:, run] 	 = self.my_city.total_power_delayed
			self.total_energy_delayed[run] 		 = self.my_city.total_energy_delayed
			
			
			# Analysis of simulation vs statistics per run
			if analyze_building_stock:
				self.analyze_building_stock_per_run(run, result_dir_run)
			if analyze_thermal_properties:
				self.analyze_thermal_properties_per_run(run, result_dir_run)
					
		# plot simulation results
		if (self.plot >= 1):
			self.plot_power(space_heating = self._space_heating, hot_water = self._hot_water, total = True)
			self.plot_energy(space_heating = self._space_heating, hot_water = self._hot_water, total = True)
			
		# save simulation results
		if (self.save >= 1):
			self.save_csv_power()
			self.save_csv_energy()
			
		# Analysis of simulation vs statistics per simulation
		if analyze_building_stock:
			self.analyze_building_stock_per_sim()
		if analyze_thermal_properties:
			self.analyze_thermal_properties_per_sim()
#	
	def create_city_object(self, run, result_dir_run):
		"""
		Creates instance of class City
		"""
		
		SIMULATION 	 = [[self.region], [self.dt_vector, self.dt_vector_excel], [self.resolution], [self.processes], [run, result_dir_run]]
		CITY 		 = [[self.Tamb, self.I], [self.buildings, self.building_stock_stats, 
						self.building_typology], [self.connection_factor], [self._space_heating, self._hot_water],
						[self.base_load]]
		
		self.my_city = City(SIMULATION, CITY, self.SPACE_HEATING, self.HOT_WATER, self.REPORTING)

# Input data
# --------------------------------------------------------------------------------
	def read_input_data_csv(self):
		"""
		Returns input data in csv files as numpy arrays.
		"""
		
		input_dir = self.my_dir + '/input/'
		
		if self.debug  != 0:
			print('Statistical data (csv)')
			print('  ' + input_dir)
		
		# Input files
		# For residential buildings
		Area_csv  		= input_dir + 'Area.csv'
		Ratio_csv 		= input_dir + 'AreaRatio.csv'
		Floors_csv		= input_dir + 'Floors.csv'
		Windows_csv		= input_dir + 'WindowOrientationRatio.csv'
		U_L1_csv  		= input_dir + 'U_Level1.csv'
		U_L2_csv  		= input_dir + 'U_Level2.csv'
		U_L3_csv  		= input_dir + 'U_Level3.csv'
		V_L1_csv  		= input_dir + 'AirFlowRate_Level1.csv'
		V_L2_csv  		= input_dir + 'AirFlowRate_Level2.csv'
		V_L3_csv  		= input_dir + 'AirFlowRate_Level3.csv'
		C_csv     		= input_dir + 'C.csv'
		Sres_csv  		= input_dir + 'BuildingStock_Unterhaching.csv'
		refurbished_csv = input_dir + 'PercentageRefurbished.csv'
		dwelling_csv	= input_dir + 'SingleDwellingBuildings.csv'
		dsize_csv 		= input_dir + 'AverageDwellingSize.csv'
		household_csv	= input_dir + 'HouseholdSize.csv'
		dhw_demand_csv	= input_dir + 'dhw_Demand.csv'
		dhw_loads_csv	= input_dir + 'dhw_Loads.csv'
		dhw_pdayt_csv	= input_dir + 'dhw_ProbDaytime.csv'
		dhw_pwday_csv	= input_dir + 'dhw_ProbWeekday.csv'
		# For non-residential buildings
		U_NR_csv  		= input_dir + 'U_NonResidential.csv'
		V_NR_csv  		= input_dir + 'AirFlowRate_NonResidential.csv'
		C_NR_csv  		= input_dir + 'C_NonResidential.csv'
		Snres_csv 		= input_dir + 'BuildingStock_NonResidential.csv'
		# For both
		Tset_csv  		= input_dir + 'Tset.csv'
		Sched_csv 		= input_dir + 'ActiveHours.csv'
		
		#
		if self.debug == 2:
			print('    ' + 'Area.csv')
			print('    ' + 'Ratio.csv')
			print('    ' + 'Floors.csv')
			print('    ' + 'U_Level1.csv')
			print('    ' + 'U_Level2.csv')
			print('    ' + 'U_Level3.csv')
			print('    ' + 'C.csv')
			print('    ' + 'BuildingStock.csv')
			print('    ' + 'U_NonResidential.csv')
			print('    ' + 'C_NonResidential.csv')
			print('    ' + 'BuildingStock_NonResidential.csv')
			print('    ' + 'Tset.csv')
			print('    ' + 'ActiveHours.csv' + '\n')
			
		# For residential buildings
		# --------------------------------------------------
		###
		# AREA [m2]
		# From TABULA Web Tool
		FOOTPRINT = self.read_data_from_csv(Area_csv, usecols = range(1, 5))   # Footprint area [m2]
		WALL_R	  = self.read_data_from_csv(Ratio_csv, usecols = range(1, 5))  # Ratio (Wall area/4)/FOOTPRINT
		ROOF_R	  = self.read_data_from_csv(Ratio_csv, usecols = range(5, 9))  # Ratio Roof area/FOOTPRINT
		WINDOW_R  = self.read_data_from_csv(Ratio_csv, usecols = range(9, 13)) # Ratio Window area/FOOTPRINT
		
		# Window ratio per orientation
		# Derived from TABULA Web Tool for residential buildings
		# >> Source missing for non-residential buildings
		WRATIO_ORIENTATION = self.read_data_from_csv(Windows_csv, usecols = range(0, 20))[0] # Ratio Window area/FOOTPRINT
		
		# NUMBER OF FLOORS
		# Derived from TABULA Web Tool
		FLOORS 	  = self.read_data_from_csv(Floors_csv, usecols = range(1, 5)) # Number of Floors
		
		###
		# U-VALUE (thermal transmittance) [W/(K m2)]
		# From TABULA Web Tool
		# Refurbishment levels:
		#   1   National minimum requirement
		#	2	Improved standard
		# 	3	Ambitious standard/NZEB
		U_ROOF    = [self.read_data_from_csv(U_L1_csv, usecols = range(1, 5)),   # L1
					 self.read_data_from_csv(U_L2_csv, usecols = range(1, 5)),   # L2
					 self.read_data_from_csv(U_L3_csv, usecols = range(1, 5))]   # L3
		U_WALL    = [self.read_data_from_csv(U_L1_csv, usecols = range(5, 9)),   # L1
					 self.read_data_from_csv(U_L2_csv, usecols = range(5, 9)),   # L2
					 self.read_data_from_csv(U_L3_csv, usecols = range(5, 9))]   # L3
		U_FLOOR   = [self.read_data_from_csv(U_L1_csv, usecols = range(9, 13)),  # L1
					 self.read_data_from_csv(U_L2_csv, usecols = range(9, 13)),  # L2
					 self.read_data_from_csv(U_L3_csv, usecols = range(9, 13))]  # L3
		U_WINDOW  = [self.read_data_from_csv(U_L1_csv, usecols = range(13, 17)), # L1
					 self.read_data_from_csv(U_L2_csv, usecols = range(13, 17)), # L2
					 self.read_data_from_csv(U_L3_csv, usecols = range(13, 17))] # L3
		
		###
		# AIR FLOW RATE (for ventilation losses) [1/h]
		# From TABULA Web Tool
		# Refurbishment levels:
		#   1   National minimum requirement
		#	2	Improved standard
		# 	3	Ambitious standard/NZEB
		V_USAGE   = [self.read_data_from_csv(V_L1_csv, usecols = range(1, 5)),   # L1
					 self.read_data_from_csv(V_L2_csv, usecols = range(1, 5)),   # L2
					 self.read_data_from_csv(V_L3_csv, usecols = range(1, 5))]   # L3
		V_INF     = [self.read_data_from_csv(V_L1_csv, usecols = range(5, 9)),   # L1
					 self.read_data_from_csv(V_L2_csv, usecols = range(5, 9)),   # L2
					 self.read_data_from_csv(V_L3_csv, usecols = range(5, 9))]   # L3
		
		###
		# THERMAL MASS (heat capacity) [J/(K m2)]
		# From http://publications.lib.chalmers.se/records/fulltext/170378/170378.pdf
		# Refurbishment levels are not considered yet 
		C_ROOF    = self.read_data_from_csv(C_csv, usecols = range(1, 5))
		C_WALL    = self.read_data_from_csv(C_csv, usecols = range(5, 9))			 
		C_FLOOR   = self.read_data_from_csv(C_csv, usecols = range(9, 13))
		
		###
		# BUILDING STOCK
		# From http://episcope.eu/fileadmin/tabula/public/docs/scientific/DE_TABULA_ScientificReport_IWU.pdf
		STOCK_RES = self.read_data_from_csv(Sres_csv, usecols = range(1, 5))
				
		###
		# PERCENTAGE OF THERMALLY REFURFISHED ENVELOPE AREAS
		# National statistics (2010) from http://episcope.eu/building-typology/country/de/
		# Data only for two age classes: < 1978 and 1979 - 1994
		PERC_REFURBISHED = self.read_data_from_csv(refurbished_csv, usecols = range(1, 5))
		
		###
		# PERCENTAGE OF SFH AND TH WITH ONE DWELLING
		# Derived from https://ergebnisse.zensus2011.de/#StaticContent:091840148148,GWZ_1_2_2,m,table
		SINGLE_DWELLING = self.read_data_from_csv(dwelling_csv, usecols = range(1, 3))
		
		###
		# AVERAGE DWELLING SIZE [m2]
		# Extracted from https://ergebnisse.zensus2011.de/#StaticContent:09,GWZ_11_11,m,table
		AVG_DWELLING_SIZE = self.read_data_from_csv(dsize_csv, usecols = 1)
		
		###
		# HOUSEHOLD SIZE
		# From https://ergebnisse.zensus2011.de/#StaticContent:091840148148,GWZ_4_3_2,m,table
		HOUSEHOLD_SIZE = self.read_data_from_csv(household_csv, usecols = range(1, 7))
		
		###
		# SPECIFIC DAILY HOT WATER DEMAND [m3/m2 of living area]
		# From VDI 3807-3, section 6.2
		DHW_DEMAND = self.read_data_from_csv(dhw_demand_csv, usecols = range(0, 2))
		
		###
		# DHW-LOAD FLOW RATE [m3/min]
		# Adapted from http://sel.me.wisc.edu/trnsys/trnlib/iea-shc-task26/iea-shc-task26-load-profiles-description-jordan.pdf
		DHW_LOAD_FLOWRATE = self.read_data_from_csv(dhw_loads_csv, usecols = range(1, 5))[0:2]
		
		###
		# DHW-LOAD DURATION [min]
		# Adapted from http://sel.me.wisc.edu/trnsys/trnlib/iea-shc-task26/iea-shc-task26-load-profiles-description-jordan.pdf
		DHW_LOAD_DURATION = self.read_data_from_csv(dhw_loads_csv, usecols = range(1, 5))[2:4]
		
		###
		# PROBABILITY DISTRIBUTION OF DHW-LOADS DURING THE DAY
		# From http://sel.me.wisc.edu/trnsys/trnlib/iea-shc-task26/iea-shc-task26-load-profiles-description-jordan.pdf
		DHW_PDAYTIME = self.read_data_from_csv(dhw_pdayt_csv, usecols = range(0, 4))
		
		###
		# FACTORS FOR DHW-LOAD PROBABILITY DISTRIBUTION DURING THE WEEK
		# From http://sel.me.wisc.edu/trnsys/trnlib/iea-shc-task26/iea-shc-task26-load-profiles-description-jordan.pdf
		DHW_PWDAY = self.read_data_from_csv(dhw_pwday_csv, usecols = range(0, 4))
		# --------------------------------------------------
		
		# For non-residential buildings
		# --------------------------------------------------
		###
		# U-VALUE (thermal transmittance) [W/(K m2)]
		# From http://www.bbsr.bund.de/BBSR/DE/Veroeffentlichungen/BMVBS/Online/2011/DL_ON162011.pdf;jsessionid=E5A73B181D25945E0C3D9B5202D4A175.live21303?__blob=publicationFile&v=2
		# Refurbishment level is not considered
		U_NRES	   = self.read_data_from_csv(U_NR_csv, usecols = range(1, 5))
		
		###
		# AIR FLOW RATE (for ventilation losses) [1/h]
		# Adapted from TABULA Web Tool
		# Refurbishment level is not considered
		V_NRES	   = self.read_data_from_csv(V_NR_csv, usecols = range(1, 3))
		
		###
		# THERMAL MASS (heat capacity) [J/(K m2)]
		# Adapted from values from residential buildings
		# Refurbishment level is not considered
		# >>> Additional source missing
		C_NRES	   = self.read_data_from_csv(C_NR_csv, usecols = range(1, 4))
		
		###
		# BUILDING STOCK
		# From http://episcope.eu/fileadmin/tabula/public/docs/scientific/DE_TABULA_ScientificReport_IWU.pdf
		STOCK_NRES = self.read_data_from_csv(Snres_csv, usecols = (1))
		# --------------------------------------------------
		
		# For both residential and non-residential
		# --------------------------------------------------
		###
		# SET TEMPERATURE [*C]
		# From http://tc76.org/spc100/docs/IBP%2018599/18599-10.pdf
		# Target temperature and dT (how much Tset varies) are defined for every building use
		TSET      = self.read_data_from_csv(Tset_csv, usecols = (1))
		D_TSET    = self.read_data_from_csv(Tset_csv, usecols = (2))

		###
		# BUILDING SCHEDULE (Active hours) [h]
		# From http://tc76.org/spc100/docs/IBP%2018599/18599-10.pdf
		# Activity start and end hours as well as dt (how much the building schedule varies) 
		# are defined for every building use
		BSCHED    = self.read_data_from_csv(Sched_csv, usecols = (1, 2))
		D_BSCHED  = self.read_data_from_csv(Sched_csv, usecols = (3))
		# --------------------------------------------------
		
		for_space_heating_demand  = [FOOTPRINT, 							# 0
									 [WALL_R, ROOF_R, WINDOW_R], 			# 1
									 WRATIO_ORIENTATION, 					# 2
									 FLOORS, 								# 3
									 [U_ROOF, U_WALL, U_FLOOR, U_WINDOW],	# 4
									 [V_USAGE, V_INF],						# 5
									 [C_ROOF, C_WALL, C_FLOOR],				# 6
									 PERC_REFURBISHED, 						# 7
									 SINGLE_DWELLING, 						# 8
									 AVG_DWELLING_SIZE, 					# 9
									 HOUSEHOLD_SIZE,						# 10
									 U_NRES, 								# 11
									 V_NRES, 								# 12
									 C_NRES, 								# 13
									 [TSET, D_TSET], 						# 14
									 [BSCHED, D_BSCHED], 					# 15
									 STOCK_RES, 							# 16
									 STOCK_NRES]							# 17
		
		for_hot_water_demand	  =	[DHW_DEMAND, 							# 0
									 DHW_PDAYTIME, 							# 1
									 DHW_PWDAY,								# 2
									 DHW_LOAD_FLOWRATE, 					# 3
									 DHW_LOAD_DURATION]						# 4
		
		self.building_stock_stats = [for_space_heating_demand,
									 for_hot_water_demand]
#
	def read_building_data(self):
		"""
		Reads building data from csv file. Returns a pd.DataFrame.
		Columns are renamed to variables in UrbanHeatPro.
		"""
		
		# Filename
		input_dir 		= self.my_dir + '/input/Buildings/'
		filename 	  	= input_dir + self.buildings_filename
				
		# Building data
		self.buildings	= pd.read_csv(filename, delimiter = ";", skiprows = 1)
		if self.debug  != 0:
			print('\nBuilding data (csv)')
			print('  ' + filename)
			print('   {} buildings\n'.format(len(self.buildings)))
			
		# rename columns
		### Required columns
		new_columns = {'bid'  	: 'bid',
					   'area' 	: 'footprint_area',
					   'use'  	: 'use',
					   'walls'	: 'free_walls',
					   'lat'  	: 'lat',
					   'lon'  	: 'lon',
					   'dist' 	: 'dist_to_heat_source',
					   'NUTS_0'	: 'NUTS_0',
					   'NUTS_1'	: 'NUTS_1',
					   'NUTS_2'	: 'NUTS_2',
					   'NUTS_3'	: 'NUTS_3'}
		self.buildings.rename(index = str, columns = new_columns)
		
		### Optional columns
		try:
			new_columns = {'year' : 'year_class'}
			self.buildings.rename(index = str, columns = new_columns)
		except:
			pass
		#
		try:
			new_columns = {'size' : 'size_class'}
			self.buildings.rename(index = str, columns = new_columns)
		except:
			pass
		#
		try:
			new_columns = {'ref_roof'  	: 'ref_level_roof',
						   'ref_wall'  	: 'ref_level_wall',
						   'ref_floor'  : 'ref_level_floor',
						   'ref_window'	: 'ref_level_window',}
			self.buildings.rename(index = str, columns = new_columns)
		except:
			pass
		#
		try:
			new_columns = {'f' : 'floors'}
			self.buildings.rename(index = str, columns = new_columns)
		except:
			pass
		#
		try:
			new_columns = {'d' : 'dwellings'}
			self.buildings.rename(index = str, columns = new_columns)
		except:
			pass
		#
		try:
			new_columns = {'occ' : 'occupants'}
			self.buildings.rename(index = str, columns = new_columns)
		except:
			pass
#
	def read_Tamb(self):
		"""
		Reads Tamb data from csv file. The file contains the Tamb values for the whole year
		in simulation resolution. Only the simulation timesteps are extracted.
		"""
		
		filename  = '{}/input/{}'.format(self.my_dir, 'Tamb.csv') 
		self.Tamb = self.read_data_from_csv(filename, usecols = 0)
		self.Tamb = self.Tamb[self.timesteps]
#
	def read_I(self):
		"""
		Reads solar radiation data from csv file. The file contains I values [W/m2] for the whole year
		in simulation resolution in the form [I_Gh, I_Dh, I_ex, hs]. Only the simulation timesteps are extracted.
		"""
		
		filename  = '{}/input/{}'.format(self.my_dir, 'I.csv') 
		self.I = self.read_data_from_csv(filename, usecols = range(0, 4))
		self.I = self.I[self.timesteps]
#
	def prepare_result_directory(self):
		"""
		Creates a time stamped directory within the result folder.
		Returns path as string.
		"""

		result_dir_ = self.my_dir + '/results' 
		now = datetime.now().strftime('%d%m%Y %H%M')
		self.result_dir = os.path.join(result_dir_, '{} - {}'.format(now, self.name))
		if not os.path.exists(self.result_dir): 
			os.makedirs(self.result_dir)
			
		if self.debug  != 0:
			print('Results directory')
			print('  ' + '{}\n'.format(self.result_dir))
#
	def read_data_from_csv(self, my_file, usecols = None):
		"""
		Uses numpy to read csv file and returns content as numpy array.
		Two rows of header are always skipped.
		"""
		
		my_data = np.genfromtxt(my_file, dtype = 'float', delimiter = ';', \
							skip_header = 2, usecols = usecols)
		
		return my_data

# Datetime vectors
# --------------------------------------------------------------------------------
	def calculate_dt_vector(self):
		"""
		Calculates a vector of datetime objects based on the raw dt_matrix of the
		form [Y, M, D, h, m] and the simulation time steps.
				
		returns:
			self.dt_vector  <list>	List of datetime objects	
		"""
		
		# read csv file
		filename   	   = '{}/input/{}'.format(self.my_dir, 'YMdhm.csv') 
		dt_matrix_year = self.read_data_from_csv(filename, usecols = (0, 1, 2, 3, 4))
		dt_matrix	   = dt_matrix_year[self.timesteps, :]
		
		# initialize vector
		self.dt_vector 		 = [[] for dt in range(len(dt_matrix))]
		self.dt_vector_excel = [[] for dt in range(len(dt_matrix))]
		
		# calculate datetime objects
		for dt, iii in enumerate(dt_matrix):
			year       = int(round(iii[0], 1))
			month      = int(round(iii[1], 1))
			day        = int(round(iii[2], 1))
			hour       = int(round(iii[3], 1))
			minute     = int(round(iii[4], 1))
			self.dt_vector[dt] = datetime(year, month, day, hour, minute)
			self.dt_vector_excel[dt] = self.convert_datetime_to_excel_date(self.dt_vector[dt])
			
		if self.debug  != 0:
			print('\nSimulation data')
			print('  ' + 'Simulation steps: {}'.format(len(self.dt_vector)))
			print('  ' + 'Resolution [min]: {}'.format(self.resolution))
			print('  ' + 'Simulation runs : {}'.format(self.N) + '\n')
#
	def convert_datetime_to_excel_date(self, dt):
		"""
		"""
		
		temp = datetime(1899, 12, 30)
		delta = dt - temp
		dt_excel = float(delta.days) + (float(delta.seconds) / 86400)
		return dt_excel
		
# Results		
# --------------------------------------------------------------------------------
	def plot_power(self, space_heating = True, hot_water = True, total = True):
		"""
		Plot min, max, and mean power values for each time step.
		"""
		if space_heating:
			min  = self.space_heating_power.min(axis = 1)
			max  = self.space_heating_power.max(axis = 1)
			mean = self.space_heating_power.mean(axis = 1)
			
			fig_name = '{}/SpaceHeatingDemand_ts.png'.format(self.result_dir)
			
			UrbanHeatPro.plot_timeseries(self.dt_vector, \
					[min, max, mean], ['min', 'max', 'mean'], \
					fig_name, xticks = ('month', 3), \
					ynumticks = 'auto', ylabel = 'Space Heating Demand [MW]', ylim0 = True, yfactor = 1e6)
					
		if hot_water:
			min  = self.hot_water_power.min(axis = 1)
			max  = self.hot_water_power.max(axis = 1)
			mean = self.hot_water_power.mean(axis = 1)
			
			fig_name = '{}/HotWaterDemand_ts.png'.format(self.result_dir)
			
			UrbanHeatPro.plot_timeseries(self.dt_vector, \
					[min, max, mean], ['min', 'max', 'mean'], \
					fig_name, xticks = ('month', 3), \
					ynumticks = 'auto', ylabel = 'Hot Water Demand [MW]', ylim0 = True, yfactor = 1e6)			

		if total:
			min  = self.total_power.min(axis = 1)
			max  = self.total_power.max(axis = 1)
			mean = self.total_power.mean(axis = 1)
			
			fig_name = '{}/TotalHeatDemand_ts.png'.format(self.result_dir)
			
			UrbanHeatPro.plot_timeseries(self.dt_vector, \
					[min, max, mean], ['min', 'max', 'mean'], \
					fig_name, xticks = ('month', 3), \
					ynumticks = 'auto', ylabel = 'Heat Demand [MW]', ylim0 = True, yfactor = 1e6)
#
	def plot_energy(self, space_heating = True, hot_water = True, total = True):
		"""
		Plots histogram of aggregated heat demand for all simulations
		"""
		
		if space_heating:
			ylabel = 'Space Heating Demand [GWh]'
			fig_name = '{}/SpaceHeatingDemand_hist.png'.format(self.result_dir)
			UrbanHeatPro.plot_histogram(self.space_heating_energy, ylabel, fig_name, factor = 1e9)
			
		if hot_water:
			ylabel = 'Hot Water Demand [GWh]'
			fig_name = '{}/HotWaterDemand_hist.png'.format(self.result_dir)
			UrbanHeatPro.plot_histogram(self.hot_water_energy, ylabel, fig_name, factor = 1e9)
			
		if total:
			ylabel = 'Total Heat Demand [GWh]'
			fig_name = '{}/TotalHeatDemand_hist.png'.format(self.result_dir)
			UrbanHeatPro.plot_histogram(self.total_energy, ylabel, fig_name, factor = 1e9)		
#
	def save_csv_power(self):
		"""
		Saves heat demand timeseries in csv files (space heating, hot water and total).
		"""
		
		### Space heating demand
		filename = '{}/SpaceHeatingDemand.csv'.format(self.result_dir)
		with open(filename, 'w') as f:
			np.savetxt(f, self.space_heating_power, delimiter = ';', fmt='%.3f')
		
		### Hot water demand
		filename = '{}/HotWaterDemand.csv'.format(self.result_dir)
		with open(filename, 'w') as f:
			np.savetxt(f, self.hot_water_power, delimiter = ';', fmt='%.3f')
			
		### Total heat demand
		filename = '{}/TotalHeatDemand.csv'.format(self.result_dir)
		with open(filename, 'w') as f:
			np.savetxt(f, self.total_power, delimiter = ';', fmt='%.3f')
			
		### Total heat demand (delayed)
		filename = '{}/TotalHeatDemand_delayed.csv'.format(self.result_dir)
		with open(filename, 'w') as f:
			np.savetxt(f, self.total_power_delayed, delimiter = ';', fmt='%.3f')
#
	def save_csv_energy(self):
		"""
		Saves key building parameters and heat energy demand (space heating, hot water and 
		total).
		"""
		
		array_to_save = [self.space_heating_energy, 
						 self.hot_water_energy, 
						 self.total_energy]
		
		filename = '{}/EnergyPerRun.csv'.format(self.result_dir)
		
		# Header
		with open(filename, 'w') as text_file:
			text_file.write('SpaceHeatingDemand_Wh;HotWaterDemand_Wh;' +
							'TotalHeatDemand_Wh\n')
				
		# Energy per run
		with open(filename, 'a') as f:
			np.savetxt(f, array_to_save, delimiter = ';', fmt='%.3f')

# Analysis
# --------------------------------------------------------------------------------
	def analyze_thermal_properties_per_run(self, run, result_dir_run):
		"""
		Extracts data from object City and plots a histogram showing the values of the 
		thermal properties for all buildings in the city. Buildings are classified according 
		to year construction class and building type. Analyzed thermal parameters are: 
			U-value [W/K]
			Thermal mass [J/K] 
			Tau [s]
			Space heat demand per unit area [Wh/m2] 
			Hot water demand per unit area [Wh/m2]
			Total heat demand per unit area [Wh/m2]
		"""
	
		# extract data from City object
		for year_class in range(10):
			for btype in range(4):
			
				# extract all buildings B with specific year_class, btype and use
				c0  = self.my_city.energy_per_building[:, 1] == 3  # RESIDENTIAL
				c1  = self.my_city.energy_per_building[:, 2] == year_class
				c2  = self.my_city.energy_per_building[:, 3] == btype
				B   = self.my_city.energy_per_building[c0 * c1 * c2]
			
				# U-value [W/K]
				self.u_res[run][year_class][btype] = np.append(self.u_res[run][year_class][btype], \
														B[:, 13])
				# Thermal mass [J/K]
				self.c_res[run][year_class][btype] = np.append(self.c_res[run][year_class][btype], \
														B[:, 14])							
				# Tau [s]
				self.tau_res[run][year_class][btype] = np.append(self.tau_res[run][year_class][btype], \
														B[:, 15])	
				# Space heat demand per unit area [Wh/m2]
				self.sh_per_area_res[run][year_class][btype] = np.append(self.sh_per_area_res[run][year_class][btype], \
														(B[:, 16] / B[:, 5]))
				# Hot water demand per unit area [Wh/m2]
				self.hw_per_area_res[run][year_class][btype] = np.append(self.hw_per_area_res[run][year_class][btype], \
														(B[:, 17] / B[:, 5]))															 
				# Total heat demand per unit area [Wh/m2]
				self.heat_per_area_res[run][year_class][btype] = np.append(self.heat_per_area_res[run][year_class][btype], \
														(B[:, 18] / B[:, 5]))
		
		for year_class in range(5):
		
			# extract all buildings B with specific year_class, btype and use
			C0  = self.my_city.energy_per_building[:, 1] != 3  # NON-RESIDENTIAL
			c1  = self.my_city.energy_per_building[:, 2] == year_class
			B   = self.my_city.energy_per_building[c0 * c1]
														
			# U-value [W/K]
			self.u_nres[run][year_class] = np.append(self.u_nres[run][year_class], B[:, 13])
			# Thermal mass [J/K]
			self.c_nres[run][year_class] = np.append(self.c_nres[run][year_class], B[:, 14])							
			# Tau [s]
			self.tau_nres[run][year_class] = np.append(self.tau_nres[run][year_class], B[:, 15])
			# Space heat demand per unit area [Wh/m2]
			self.sh_per_area_nres[run][year_class] = np.append(self.sh_per_area_nres[run][year_class], \
														(B[:, 16] / B[:, 5]))
			# Hot water demand per unit area [Wh/m2]
			self.hw_per_area_nres[run][year_class] = np.append(self.hw_per_area_nres[run][year_class], \
														(B[:, 17] / B[:, 5]))															 
			# Total heat demand per unit area [Wh/m2]
			self.heat_per_area_nres[run][year_class] = np.append(self.heat_per_area_nres[run][year_class], \
														(B[:, 18] / B[:, 5]))													
		
		# plot histograms per run 
		### RESIDENTIAL
		###### U-value [MW/K]
		figname = '{}/histU_res_{}.png'.format(result_dir_run, run)
		UrbanHeatPro.plot_histogram_table(3, self.u_res[run], 'U-value [MW/K]', figname,
								factor = 1e6, figsize = (30, 25))
		###### U-value [MW/K]
		figname = '{}/histC_res_{}.png'.format(result_dir_run, run)
		UrbanHeatPro.plot_histogram_table(3, self.c_res[run], 'Thermal mass [GJ/K]', figname,
								factor = 1e9, figsize = (30, 25))
		###### Tau [h]
		figname = '{}/histTau_res_{}.png'.format(result_dir_run, run)
		UrbanHeatPro.plot_histogram_table(3, self.tau_res[run], 'Tau [h]', figname,
								factor = 3600, figsize = (30, 25))
		###### Space heating demand per unit area [kWh/m2]
		figname = '{}/histSHD_perUnitArea_res_{}.png'.format(result_dir_run, run)
		UrbanHeatPro.plot_histogram_table(3, self.sh_per_area_res[run], 'Energy [kWh/m2]', figname,
								factor = 1e3, figsize = (30, 25))
		###### Hot water demand per unit area [kWh/m2]
		figname = '{}/histHWD_perUnitArea_res_{}.png'.format(result_dir_run, run)
		UrbanHeatPro.plot_histogram_table(3, self.hw_per_area_res[run], 'Energy [kWh/m2]', figname,
								factor = 1e3, figsize = (30, 25))
		###### Total heat demand per unit area [kWh/m2]
		figname = '{}/histHeat_perUnitArea_res_{}.png'.format(result_dir_run, run)
		UrbanHeatPro.plot_histogram_table(3, self.heat_per_area_res[run], 'Energy [kWh/m2]', figname,
								factor = 1e3, figsize = (30, 25))
								
		### NON-RESIDENTIAL
		###### U-value [MW/K]
		figname = '{}/histU_nres_{}.png'.format(result_dir_run, run)
		UrbanHeatPro.plot_histogram_table(0, self.u_nres[run], 'U-value [MW/K]', figname,
								factor = 1e6, figsize = (30, 25))
		###### U-value [MW/K]
		figname = '{}/histC_nres_{}.png'.format(result_dir_run, run)
		UrbanHeatPro.plot_histogram_table(0, self.c_nres[run], 'Thermal mass [GJ/K]', figname,
								factor = 1e9, figsize = (30, 25))
		###### Tau [h]
		figname = '{}/histTau_nres_{}.png'.format(result_dir_run, run)
		UrbanHeatPro.plot_histogram_table(0, self.tau_nres[run], 'Tau [h]', figname,
								factor = 3600, figsize = (30, 25))
		###### Space heating demand per unit area [kWh/m2]
		figname = '{}/histSHD_perUnitArea_nres_{}.png'.format(result_dir_run, run)
		UrbanHeatPro.plot_histogram_table(0, self.sh_per_area_nres[run], 'Energy [kWh/m2]', figname,
								factor = 1e3, figsize = (30, 25))
		###### Hot water demand per unit area [kWh/m2]
		figname = '{}/histHWD_perUnitArea_nres_{}.png'.format(result_dir_run, run)
		UrbanHeatPro.plot_histogram_table(0, self.hw_per_area_nres[run], 'Energy [kWh/m2]', figname,
								factor = 1e3, figsize = (30, 25))
		###### Total heat demand per unit area [kWh/m2]
		figname = '{}/histHeat_perUnitArea_nres_{}.png'.format(result_dir_run, run)
		UrbanHeatPro.plot_histogram_table(0, self.heat_per_area_nres[run], 'Energy [kWh/m2]', figname,
								factor = 1e3, figsize = (30, 25))
#
	def analyze_thermal_properties_per_sim(self):
		"""
		Appends thermal properties for all runs and plot histograms
		"""
		
		# initialize temporal matrices
		U_res     = [[np.array([]) for col in range(4)] for row in range(10)]
		C_res     = [[np.array([]) for col in range(4)] for row in range(10)]
		TAU_res   = [[np.array([]) for col in range(4)] for row in range(10)]
		SH_res	  = [[np.array([]) for col in range(4)] for row in range(10)]
		HW_res	  = [[np.array([]) for col in range(4)] for row in range(10)]
		HEAT_res  = [[np.array([]) for col in range(4)] for row in range(10)]
		U_nres    = [np.array([]) for col in range(5)]
		C_nres    = [np.array([]) for col in range(5)]
		TAU_nres  = [np.array([]) for col in range(5)]
		SH_nres	  = [np.array([]) for col in range(5)]
		HW_nres	  = [np.array([]) for col in range(5)]
		HEAT_nres = [np.array([]) for col in range(5)]
		
		# append simulation values
		for run in range(self.N):
			
			# RESIDENTIAL
			for year_class in range(10):
				for btype in range(4):
				
					U_res[year_class][btype]    = np.append(U_res[year_class][btype], self.u_res[run][year_class][btype])
					C_res[year_class][btype]    = np.append(C_res[year_class][btype], self.c_res[run][year_class][btype])
					TAU_res[year_class][btype]  = np.append(TAU_res[year_class][btype], self.tau_res[run][year_class][btype])
					SH_res[year_class][btype]   = np.append(SH_res[year_class][btype], self.sh_per_area_res[run][year_class][btype])
					HW_res[year_class][btype]   = np.append(HW_res[year_class][btype], self.hw_per_area_res[run][year_class][btype])
					HEAT_res[year_class][btype] = np.append(HEAT_res[year_class][btype], self.heat_per_area_res[run][year_class][btype])
			
			# NON-RESIDENTIAL
			for year_class in range(5):
				
				U_nres[year_class] 		 = np.append(U_nres[year_class], self.u_nres[run][year_class])
				C_nres[year_class]		 = np.append(C_nres[year_class], self.c_nres[run][year_class])
				TAU_nres[year_class]	 = np.append(TAU_nres[year_class], self.tau_nres[run][year_class])
				SH_nres[year_class]		 = np.append(SH_nres[year_class], self.sh_per_area_nres[run][year_class])
				HW_nres[year_class]		 = np.append(HW_nres[year_class], self.hw_per_area_nres[run][year_class])
				HEAT_nres[year_class]	 = np.append(HEAT_nres[year_class], self.heat_per_area_nres[run][year_class])
		
		# plot histograms
		### RESIDENTIAL
		###### U-value [MW/K]
		figname = '{}/histU_res.png'.format(self.result_dir)
		UrbanHeatPro.plot_histogram_table(3, U_res, 'U-value [MW/K]', figname,
								factor = 1e6, figsize = (30, 25))
		###### U-value [MW/K]
		figname = '{}/histC_res.png'.format(self.result_dir)
		UrbanHeatPro.plot_histogram_table(3, C_res, 'Thermal mass [GJ/K]', figname,
								factor = 1e9, figsize = (30, 25))
		###### Tau [h]
		figname = '{}/histTau_res.png'.format(self.result_dir)
		UrbanHeatPro.plot_histogram_table(3, TAU_res, 'Tau [h]', figname,
								factor = 3600, figsize = (30, 25))
		###### Space heating demand per unit area [kWh/m2]
		figname = '{}/histSHD_perUnitArea_res.png'.format(self.result_dir)
		UrbanHeatPro.plot_histogram_table(3, SH_res, 'Energy [kWh/m2]', figname,
								factor = 1e3, figsize = (30, 25))
		###### Hot water demand per unit area [kWh/m2]
		figname = '{}/histHWD_perUnitArea_res.png'.format(self.result_dir)
		UrbanHeatPro.plot_histogram_table(3, HW_res, 'Energy [kWh/m2]', figname,
								factor = 1e3, figsize = (30, 25))
		###### Total heat demand per unit area [kWh/m2]
		figname = '{}/histHeat_perUnitArea_res.png'.format(self.result_dir)
		UrbanHeatPro.plot_histogram_table(3, HEAT_res, 'Energy [kWh/m2]', figname,
								factor = 1e3, figsize = (30, 25))
								
		### NON-RESIDENTIAL
		###### U-value [MW/K]
		figname = '{}/histU_nres.png'.format(self.result_dir)
		UrbanHeatPro.plot_histogram_table(0, U_nres, 'U-value [MW/K]', figname,
								factor = 1e6, figsize = (30, 25))
		###### U-value [MW/K]
		figname = '{}/histC_nres.png'.format(self.result_dir)
		UrbanHeatPro.plot_histogram_table(0, C_nres, 'Thermal mass [GJ/K]', figname,
								factor = 1e9, figsize = (30, 25))
		###### Tau [h]
		figname = '{}/histTau_nres.png'.format(self.result_dir)
		UrbanHeatPro.plot_histogram_table(0, TAU_nres, 'Tau [h]', figname,
								factor = 3600, figsize = (30, 25))
		###### Space heating demand per unit area [kWh/m2]
		figname = '{}/histSHD_perUnitArea_nres.png'.format(self.result_dir)
		UrbanHeatPro.plot_histogram_table(0, SH_nres, 'Energy [kWh/m2]', figname,
								factor = 1e3, figsize = (30, 25))
		###### Hot water demand per unit area [kWh/m2]
		figname = '{}/histHWD_perUnitArea_nres.png'.format(self.result_dir)
		UrbanHeatPro.plot_histogram_table(0, HW_nres, 'Energy [kWh/m2]', figname,
								factor = 1e3, figsize = (30, 25))
		###### Total heat demand per unit area [kWh/m2]
		figname = '{}/histHeat_perUnitArea_nres.png'.format(self.result_dir)
		UrbanHeatPro.plot_histogram_table(0, HEAT_nres, 'Energy [kWh/m2]', figname,
								factor = 1e3, figsize = (30, 25))
#
	def analyze_building_stock_per_run(self, run, result_dir_run):
		"""
		Compares the distribution of the synthetic building stock generated per run
		with the statistics.
		"""
			
		# RESIDENTIAL
		# total buildings
		count_res  = len(self.my_city.energy_per_building[self.my_city.energy_per_building[:, 1] == 3])
		
		for year_class in range(10):
			for btype in range(4):
		
				# extract all buildings B with specific year_class, btype and use
				c0  = self.my_city.energy_per_building[:, 1] == 3  # RESIDENTIAL
				c1  = self.my_city.energy_per_building[:, 2] == year_class
				c2  = self.my_city.energy_per_building[:, 3] == btype
				B   = self.my_city.energy_per_building[c0 * c1 * c2]
			
				# add building distribution to category
				self.bstock_res[run][year_class][btype]   = len(B)/float(count_res)
		
		# NON-RESIDENTIAL
		count_nres = len(self.my_city.energy_per_building[self.my_city.energy_per_building[:, 1] != 3])
		
		for year_class in range(5):
			
			# extract all buildings B with specific year_class, btype and use
			c0  = self.my_city.energy_per_building[:, 1] != 3  # NON-RESIDENTIAL
			c1  = self.my_city.energy_per_building[:, 2] == year_class
			B   = self.my_city.energy_per_building[c0 * c1]
			
			# add building distribution to category
			self.bstock_nres[run][year_class] 			  = len(B)/float(count_nres)
			count_nres += 1
					
		# Plots
		# Building distribution per category from statistics [%]
		stat_res   = self.building_stock_stats[8]
		stat_nres  = self.building_stock_stats[9]
		
		# plot imshow to compare with statistics
		fig_name    = '{}/BuildingStock_res_{}.png'.format(result_dir_run, run)
		UrbanHeatPro.plot_imshow_comparison(3, self.bstock_res[run], stat_res, fig_name, cmap = 'RdBu')
		#
		fig_name    = '{}/BuildingStock_nres_{}.png'.format(result_dir_run, run)
		UrbanHeatPro.plot_imshow_comparison(0, self.bstock_nres[run], stat_nres, fig_name, cmap = 'RdBu')
#
	def analyze_building_stock_per_sim(self):
		"""
		Compares the MEAN distribution of the synthetic building stock generated per run
		with the statistics.
		"""
		
		# Building distribution per category from statistics [%]
		stat_res   = self.building_stock_stats[8]
		stat_nres  = self.building_stock_stats[9]
		
		# Simulation building stock results
		mean_res   = self.bstock_res.mean(axis = 0)
		mean_nres  = self.bstock_nres.mean(axis = 0)
		
		# plot imshow to compare with statistics
		fig_name    = '{}/BuildingStock_res.png'.format(self.result_dir)
		UrbanHeatPro.plot_imshow_comparison(3, mean_res, stat_res, fig_name, cmap = 'RdBu')
		#
		fig_name    = '{}/BuildingStock_nres.png'.format(self.result_dir)
		UrbanHeatPro.plot_imshow_comparison(0, mean_nres, stat_nres, fig_name, cmap = 'RdBu')