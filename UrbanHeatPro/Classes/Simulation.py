"""
Simulation.py
AMC @ TUM ENS
"""

import os
import numpy as np
import pandas as pd
from datetime import datetime
from random import randint
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import squareform, pdist
import copy
#import ipdb

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
		self.timesteps			  = SIMULATION[1][2]		# Vector of timesteps
		self.number_of_typ_days   = SIMULATION[1][3]		# Number of typical days to simulate
		self.weights			  = np.ones(self.number_of_typ_days) # Weight of typical days
		self.dt_vector    		  = None		 			# Vector of time steps as datetime objects
		self.dt_vector_excel	  = None					# Vector of time steps as excel date
		self.my_dir    			  = os.getcwd()				# Current directory
		self.sce_refurbishment	  = SIMULATION[2][0]		# Name of scenario for refurbishment stats
		self.sce_Tamb			  = SIMULATION[2][1]		# Name of scenario for ambient temperature
		self.processes			  = SIMULATION[3][0]		# Number of parallel processes
		self.chunk_size			  = SIMULATION[3][1]		# Number of buildings in chunk to save
		self.sce				  = None
		
		# CITY
		### External factor
		self.Tamb				  = None					# Vector of ambient temperature in degC
		self.I				  	  = None					# Solar radiation vector in W/m2 [I_Gh, I_Dh, I_ex, hs]
		### Buildings
		self.buildings_filename	  = CITY[0][0]				# Name of file containing raw building data
		self.syncity_filename	  = CITY[0][1]				# Name of file containing syncity data
		self.building_stock_stats = None					# Data from building stock statistics
		self.buildings			  = None					# DataFrame with buildings
		self.connection_factor	  = CITY[1][0]				# Share of buildings connected to the network
		### Flags
		self._space_heating	      = CITY[2][0]				# calculate space heating demand?
		self._hot_water		      = CITY[2][1]				# calculate hot water demand?
		self._energy_only		  = CITY[2][2]				# calculate only aggregated demand?
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
	def run(self, include_date = True):
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
		
		## Building data 			-> self.buildings
		### Raw buildings
		if self.buildings_filename is not None: 
			self.read_raw_building_data()
		### Synthetic city	
		if self.syncity_filename is not None:
			input_dir 		= self.my_dir + '/input/SynCity/'
			filename 	  	= input_dir + self.syncity_filename
			self.buildings  = self.read_syn_city(filename)
		
		### Statistical data 		-> self.building_stock_stats
		self.read_input_data_csv()
		
		### Result directory		-> self.result_dir
		self.prepare_result_directory(include_date = include_date)
		
		### Weather data			-> self.Tamb, self.I
		self.read_Tamb()
		self.read_I()

		### Simulation time steps 	-> self.dt_vector
		# calculate days
		days = len(self.timesteps) / 24
		if (self.number_of_typ_days < days):
			self.timesteps, self.weights = self.calculate_typical_days()
		self.calculate_dt_vector()
		self.nts = len(self.dt_vector)
		self.filter_weather_data()
		
		
		# initialize results matrices
		### Space heating demand
		self.space_heating_power  = np.zeros([self.nts, self.N]) # City Space heating demand in W
		### Hot water demand
		self.hot_water_power      = np.zeros([self.nts, self.N]) # City hot water demand in W
		### Total heat demand
		self.total_power   		  = np.zeros([self.nts, self.N]) # City total heat demand in W
		### Total heat demand (delayed)
		self.total_power_delayed  = np.zeros([self.nts, self.N]) # City total heat demand in W
		
		# MAIN
		# ----------------------------------------------------------------
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
				
			# modify chunk_size
			#print('chunk_size old: {}'.format(self.chunk_size))
			#self.chunk_size = int(np.floor(len(self.buildings) // np.floor(len(self.buildings) // self.chunk_size)))
			#print('chunk_size new: {}'.format(self.chunk_size))
			
			# create City object
			self.create_city_object(run, result_dir_run)
			
			# create synthetic city
			if self.syncity_filename is None:  
				if self.debug  != 0:
					print('\nSYNTHETIC CITY')
					print('---------------------------------------------------')
				self.my_city.create_synthetic_city()
				if self.debug  != 0:
					print('\nSynCity created! ***')
				
				# read created syn city to update building dataframe
				#filename = '{}/{}/SynCity_{}_{}.csv'.format(self.result_dir, run, self.region, run)
				filename = '{}/{}/SynCity_{}_{}.csv'.format(self.result_dir, run, self.name, run)
				self.buildings = self.read_syn_city(filename)
				
				# update city object
				self.my_city.buildings = self.buildings
				
			# update synthetic city with scenario data
			## update refurbishment levels
			if self.debug  != 0:
				if self.sce_refurbishment or self.sce_Tamb:
					print('\nSCENARIOS')
					print('---------------------------------------------------')
			if self.sce_refurbishment:
				if self.debug  != 0:
					print('Refurbishment scenario: {}'.format(self.sce_refurbishment))
				self.read_refurbishment_matrices()
				self.my_city.sce = self.sce # update sce name
				self.my_city.update_synthetic_city(self.ref_matrix_res, self.ref_matrix_nres)
				if self.debug  != 0:
					print('\nSynCity updated: Refurbishment scenario ***')
			## update temperature vector
			if self.sce_Tamb:
				if self.debug  != 0:
					print('Tamb scenario: {}'.format(self.sce_Tamb))
				self.update_Tamb()
				if self.debug  != 0:
					print('\nSynCity updated: Temperature scenario ***')
				
			# calculate heating energy demand
			if self.debug  != 0:
				print('\nHEAT DEMAND')
				print('---------------------------------------------------')
			self.my_city.calculate_city_heat_demand()			
#	
	def create_city_object(self, run, result_dir_run):
		"""
		Creates instance of class City
		"""
		
		SIMULATION 	 = [[self.region, self.sce], [self.dt_vector, self.dt_vector_excel], [self.resolution],
						[self.processes, self.chunk_size], [self.number_of_typ_days, self.weights], 
						[run, result_dir_run]]
		CITY 		 = [[self.Tamb, self.I], [self.buildings, self.building_stock_stats], 
						[self.connection_factor], [self._space_heating, self._hot_water, self._energy_only],
						[self.base_load]]
		
		self.my_city = City(self.name, SIMULATION, CITY, self.SPACE_HEATING, self.HOT_WATER, self.REPORTING)

		
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
		
		#
		# File names
		#
		# Building typology
		# --------------------------------------------------
		input_dir_ = input_dir + '/Building Typology/'
		
		## RESIDENTIAL
		### Areas
		Area_csv  		= input_dir_ + 'EnvelopeArea_Residential.csv'
		Ratio_csv 		= input_dir_ + 'AreaRatio_Residential.csv'
		Windows_csv		= input_dir_ + 'WindowOrientationRatio_Residential.csv'
		Floors_csv		= input_dir_ + 'Floors_Residential.csv'
		### Thermal properties
		U1_res_csv  	= input_dir_ + 'U-1_Residential.csv'
		U2_res_csv  	= input_dir_ + 'U-2_Residential.csv'
		U3_res_csv  	= input_dir_ + 'U-3_Residential.csv'
		C_res_csv     	= input_dir_ + 'C_Residential.csv'
		V1_res_csv  	= input_dir_ + 'AirFlowRate-1_Residential.csv'
		V2_res_csv  	= input_dir_ + 'AirFlowRate-2_Residential.csv'
		V3_res_csv  	= input_dir_ + 'AirFlowRate-3_Residential.csv'
		
		## NON-RESIDENTIAL
		### Thermal properties
		U1_nres_csv  	= input_dir_ + 'U-1_NonResidential.csv'
		U2_nres_csv  	= input_dir_ + 'U-2_NonResidential.csv'
		U3_nres_csv  	= input_dir_ + 'U-3_NonResidential.csv'
		C_nres_csv     	= input_dir_ + 'C_NonResidential.csv'
		V_nres_csv  	= input_dir_ + 'AirFlowRate_NonResidential.csv'
		
		## BOTH
		Tset_csv  		= input_dir_ + 'Tset.csv'
		sh_prob_csv		= input_dir_ + 'MonthlySpaceHeatingProbability.csv'
		
		# Domestic Hot Water
		# --------------------------------------------------
		input_dir_ = input_dir + '/Domestic Hot Water/'
		
		dhw_demand_csv	= input_dir_ + 'dhw_Demand.csv'
		dhw_loads_csv	= input_dir_ + 'dhw_Loads.csv'
		dhw_pdayt_csv	= input_dir_ + 'dhw_ProbDaytime.csv'
		dhw_pwday_csv	= input_dir_ + 'dhw_ProbWeekday.csv'
		
		# Regional data
		# --------------------------------------------------
		input_dir_ = input_dir + '/Regional Data/' + self.region + '/'
		
		## RESIDENTIAL
		stock_res_csv  	= input_dir_ + 'BuildingStock_Residential_' + self.region + '.csv'
		dwelling_csv	= input_dir_ + 'SingleDwellingBuildings_' + self.region + '.csv'
		dsize_csv 		= input_dir_ + 'AverageDwellingSize_' + self.region + '.csv'
		household_csv	= input_dir_ + 'HouseholdSize_' + self.region + '.csv'	
		ref_res_csv 	= input_dir_ + 'CurrentRefurbished_Residential_' + self.region + '.csv'
		max_ref_res_csv	= input_dir_ + 'MaxRefurbished_Residential_' + self.region + '.csv'

		## NON-RESIDENTIAL
		stock_nres_csv	= input_dir_ + 'BuildingStock_NonResidential_' + self.region + '.csv'
		ref_nres_csv 	= input_dir_ + 'CurrentRefurbished_NonResidential_' + self.region + '.csv'
		max_ref_nres_csv= input_dir_ + 'MaxRefurbished_NonResidential_' + self.region + '.csv'
		
		## BOTH
		
		Sched_csv 		= input_dir_ + 'ActiveHours_' + self.region + '.csv'
		
		#
		# read files
		#
		# Building typology
		# --------------------------------------------------
		
		## RESIDENTIAL
		
		### ENVELOPE AREA [m2]
		### From TABULA Web Tool
		FOOTPRINT = self.read_data_from_csv(Area_csv, usecols = range(1, 5))   # Footprint area [m2]
		WALL_R	  = self.read_data_from_csv(Ratio_csv, usecols = range(1, 5))  # Ratio (Wall area/4)/FOOTPRINT
		ROOF_R	  = self.read_data_from_csv(Ratio_csv, usecols = range(5, 9))  # Ratio Roof area/FOOTPRINT
		WINDOW_R  = self.read_data_from_csv(Ratio_csv, usecols = range(9, 13)) # Ratio Window area/FOOTPRINT
		
		### WINDOW RATIO PER ORIENTATION
		### Derived from TABULA Web Tool for residential buildings
		### >> Source missing for non-residential buildings
		WRATIO_ORIENTATION = self.read_data_from_csv(Windows_csv, usecols = range(0, 20))[0] # Ratio Window area/FOOTPRINT
		
		### NUMBER OF FLOORS
		### Derived from TABULA Web Tool
		FLOORS 	  = self.read_data_from_csv(Floors_csv, usecols = range(1, 5)) # Number of Floors
		
		### U-VALUE (thermal transmittance) [W/(K m2)]
		### From TABULA Web Tool
		### Refurbishment levels:
		###   	1   National minimum requirement
		###		2	Improved standard
		### 	3	Ambitious standard/NZEB
		U_ROOF_R    = [self.read_data_from_csv(U1_res_csv, usecols = range(1, 5)),  # L1
					 self.read_data_from_csv(U2_res_csv, usecols = range(1, 5)),   	# L2
					 self.read_data_from_csv(U3_res_csv, usecols = range(1, 5))]   	# L3
		U_WALL_R    = [self.read_data_from_csv(U1_res_csv, usecols = range(5, 9)),  # L1
					 self.read_data_from_csv(U2_res_csv, usecols = range(5, 9)),    # L2
					 self.read_data_from_csv(U3_res_csv, usecols = range(5, 9))]    # L3
		U_FLOOR_R   = [self.read_data_from_csv(U1_res_csv, usecols = range(9, 13)), # L1
					 self.read_data_from_csv(U2_res_csv, usecols = range(9, 13)),   # L2
					 self.read_data_from_csv(U3_res_csv, usecols = range(9, 13))]   # L3
		U_WINDOW_R  = [self.read_data_from_csv(U1_res_csv, usecols = range(13, 17)),# L1
					 self.read_data_from_csv(U2_res_csv, usecols = range(13, 17)),  # L2
					 self.read_data_from_csv(U3_res_csv, usecols = range(13, 17))]  # L3
		
		### AIR FLOW RATE (for ventilation losses) [1/h]
		### From TABULA Web Tool
		### Refurbishment levels:
		###   	1   National minimum requirement
		###		2	Improved standard
		# ##	3	Ambitious standard/NZEB
		V_USAGE_R   = [self.read_data_from_csv(V1_res_csv, usecols = range(1, 5)),  # L1
					 self.read_data_from_csv(V2_res_csv, usecols = range(1, 5)),   	# L2
					 self.read_data_from_csv(V3_res_csv, usecols = range(1, 5))]   	# L3
		V_INF_R     = [self.read_data_from_csv(V1_res_csv, usecols = range(5, 9)),  # L1
					 self.read_data_from_csv(V2_res_csv, usecols = range(5, 9)),    # L2
					 self.read_data_from_csv(V3_res_csv, usecols = range(5, 9))]    # L3
		
		### THERMAL MASS (heat capacity) [J/(K m2)]
		### From http://publications.lib.chalmers.se/records/fulltext/170378/170378.pdf
		### Refurbishment levels are not considered yet 
		C_ROOF_R    = self.read_data_from_csv(C_res_csv, usecols = range(1, 5))
		C_WALL_R    = self.read_data_from_csv(C_res_csv, usecols = range(5, 9))			 
		C_FLOOR_R   = self.read_data_from_csv(C_res_csv, usecols = range(9, 13))
		
		
		## NON-RESIDENTIAL
		
		### U-VALUE (thermal transmittance) [W/(K m2)]
		### From http://www.bbsr.bund.de/BBSR/DE/Veroeffentlichungen/BMVBS/Online/2011/DL_ON162011.pdf;jsessionid=E5A73B181D25945E0C3D9B5202D4A175.live21303?__blob=publicationFile&v=2
		### Adapted to TABULA refurbishment levels
		U_ROOF_NR     = [self.read_data_from_csv(U1_nres_csv, usecols = 1), # L1
					 self.read_data_from_csv(U2_nres_csv, usecols = 1),   	# L2
					 self.read_data_from_csv(U3_nres_csv, usecols = 1)]   	# L3
		U_WALL_NR     = [self.read_data_from_csv(U1_nres_csv, usecols = 2), # L1
					 self.read_data_from_csv(U2_nres_csv, usecols = 2),   	# L2
					 self.read_data_from_csv(U3_nres_csv, usecols = 2)]   	# L3
		U_FLOOR_NR    = [self.read_data_from_csv(U1_nres_csv, usecols = 3), # L1
					 self.read_data_from_csv(U2_nres_csv, usecols = 3),  	# L2
					 self.read_data_from_csv(U3_nres_csv, usecols = 3)]  	# L3
		U_WINDOW_NR   = [self.read_data_from_csv(U1_nres_csv, usecols = 4), # L1
					 self.read_data_from_csv(U2_nres_csv, usecols = 4), 	# L2
					 self.read_data_from_csv(U3_nres_csv, usecols = 4)] 	# L3
					
		### AIR FLOW RATE (for ventilation losses) [1/h]
		### Adapted from TABULA refurbishment levels
		V_USAGE_NR    = self.read_data_from_csv(V_nres_csv, usecols = 1)
		V_INF_NR      = self.read_data_from_csv(V_nres_csv, usecols = 2)
		
		### THERMAL MASS (heat capacity) [J/(K m2)]
		### Adapted from values from residential buildings
		### Refurbishment levels are not considered yet 
		C_ROOF_NR     = self.read_data_from_csv(C_nres_csv, usecols = 1)
		C_WALL_NR     = self.read_data_from_csv(C_nres_csv, usecols = 2)			 
		C_FLOOR_NR    = self.read_data_from_csv(C_nres_csv, usecols = 3)
		
		# --------------------------------------------------
		
		
		# Domestic Hot Water
		# --------------------------------------------------
		
		### SPECIFIC DAILY HOT WATER DEMAND [m3/m2 of living area]
		### From VDI 3807-3, section 6.2
		DHW_DEMAND = self.read_data_from_csv(dhw_demand_csv, usecols = range(0, 2))
		
		### DHW-LOAD FLOW RATE [m3/min]
		### Adapted from http://sel.me.wisc.edu/trnsys/trnlib/iea-shc-task26/iea-shc-task26-load-profiles-description-jordan.pdf
		DHW_LOAD_FLOWRATE = self.read_data_from_csv(dhw_loads_csv, usecols = range(1, 5))[0:2]
		
		### DHW-LOAD DURATION [min]
		### Adapted from http://sel.me.wisc.edu/trnsys/trnlib/iea-shc-task26/iea-shc-task26-load-profiles-description-jordan.pdf
		DHW_LOAD_DURATION = self.read_data_from_csv(dhw_loads_csv, usecols = range(1, 5))[2:4]
		
		### PROBABILITY DISTRIBUTION OF DHW-LOADS DURING THE DAY
		### From http://sel.me.wisc.edu/trnsys/trnlib/iea-shc-task26/iea-shc-task26-load-profiles-description-jordan.pdf
		DHW_PDAYTIME = self.read_data_from_csv(dhw_pdayt_csv, usecols = range(0, 4))
		
		### FACTORS FOR DHW-LOAD PROBABILITY DISTRIBUTION DURING THE WEEK
		### From http://sel.me.wisc.edu/trnsys/trnlib/iea-shc-task26/iea-shc-task26-load-profiles-description-jordan.pdf
		DHW_PWDAY = self.read_data_from_csv(dhw_pwday_csv, usecols = range(0, 4))
		
		# --------------------------------------------------
		
		
		# Regional data
		# --------------------------------------------------
		
		## RESIDENTIAL
		
		### BUILDING STOCK
		### From http://episcope.eu/fileadmin/tabula/public/docs/scientific/DE_TABULA_ScientificReport_IWU.pdf
		STOCK_RES = self.read_data_from_csv(stock_res_csv, usecols = range(1, 5))
				
		### PERCENTAGE OF SFH AND TH WITH ONE DWELLING
		### Derived from https://ergebnisse.zensus2011.de/#StaticContent:091840148148,GWZ_1_2_2,m,table
		SINGLE_DWELLING = self.read_data_from_csv(dwelling_csv, usecols = range(1, 3))
		
		### AVERAGE DWELLING SIZE [m2]
		### Extracted from https://ergebnisse.zensus2011.de/#StaticContent:09,GWZ_11_11,m,table
		AVG_DWELLING_SIZE = self.read_data_from_csv(dsize_csv, usecols = 1)
		
		### HOUSEHOLD SIZE
		### From https://ergebnisse.zensus2011.de/#StaticContent:091840148148,GWZ_4_3_2,m,table
		HOUSEHOLD_SIZE = self.read_data_from_csv(household_csv, usecols = range(1, 7))
		
		### CURRENT PERCENTAGE OF THERMALLY REFURFISHED ENVELOPE AREAS
		### Adapted from national statistics (2010) from http://episcope.eu/building-typology/country/de/
		### Additional sources missing
		data = self.read_data_from_csv(ref_res_csv, usecols = range(2, 6))
		CURRENT_REF_RES = [data[0:11,:], data[11:22,:], data[22:33,:], data[33:44,:]]

		### MAXIMUM PERCENTAGE OF POTENTIAL REFURFISHED ENVELOPE AREAS
		### Adapted from https://www.umweltbundesamt.de/sites/default/files/medien/378/publikationen/climate_change_06_2016_klimaneutraler_gebaeudebestand_2050.pdf
		# >>> Read with geopandas
		data = self.read_data_from_csv(max_ref_res_csv, usecols = range(2, 6))
		MAX_REF_RES = [data[0:11,:], data[11:22,:], data[22:33,:], data[33:44,:]]
		
		## NON-RESIDENTIAL
		
		### BUILDING STOCK
		### From >>> Source missing
		STOCK_NRES = self.read_data_from_csv(stock_nres_csv, usecols = 1)
		
		### CURRENT PERCENTAGE OF THERMALLY REFURFISHED ENVELOPE AREAS
		###  From >>> source missing
		# >>> Read with geopandas
		CURRENT_REF_NRES = self.read_data_from_csv(ref_nres_csv, usecols = range(1, 5))
		
		### MAXIMUM PERCENTAGE OF POTENTIAL REFURFISHED ENVELOPE AREAS
		### Adapted from https://www.umweltbundesamt.de/sites/default/files/medien/378/publikationen/climate_change_06_2016_klimaneutraler_gebaeudebestand_2050.pdf
		# >>> Read with geopandas
		MAX_REF_NRES = self.read_data_from_csv(max_ref_nres_csv, usecols = range(1, 5))
		
		## BOTH

		### SET TEMPERATURE [*C]
		### From http://tc76.org/spc100/docs/IBP%2018599/18599-10.pdf
		### Target temperature and dT (how much Tset varies) are defined for every building use
		TSET      = self.read_data_from_csv(Tset_csv, usecols = (1))
		D_TSET    = self.read_data_from_csv(Tset_csv, usecols = (2))

		### BUILDING SCHEDULE (Active hours) [h]
		### From http://tc76.org/spc100/docs/IBP%2018599/18599-10.pdf
		### Activity start and end hours as well as dt (how much the building schedule varies) 
		### are defined for every building use
		BSCHED    = self.read_data_from_csv(Sched_csv, usecols = (1, 2))
		D_BSCHED  = self.read_data_from_csv(Sched_csv, usecols = (3))
		
		###	MONTHLY PROBABILITY OF USING SPACE HEATING DEMAND
		###	Heating period in Germany is from 01.10 to 30.04
		### Probabilities of using heat are adjusted accordingly
		### >>> Source missing for monthly probabilities
		SH_PROB   = self.read_data_from_csv(sh_prob_csv, usecols = (1))
		
		# --------------------------------------------------
		
		#
		# organize data
		#
		for_space_heating_demand  = [FOOTPRINT, 										# 0
									 [WALL_R, ROOF_R, WINDOW_R], 						# 1
									 WRATIO_ORIENTATION, 								# 2
									 FLOORS, 											# 3
									 [U_ROOF_R, U_WALL_R, U_FLOOR_R, U_WINDOW_R],		# 4
									 [V_USAGE_R, V_INF_R],								# 5
									 [C_ROOF_R, C_WALL_R, C_FLOOR_R],					# 6
									 [CURRENT_REF_RES, MAX_REF_RES], 					# 7
									 SINGLE_DWELLING, 									# 8
									 AVG_DWELLING_SIZE, 								# 9
									 HOUSEHOLD_SIZE,									# 10
									 [U_ROOF_NR, U_WALL_NR, U_FLOOR_NR, U_WINDOW_NR],	# 11
									 [V_USAGE_NR, V_INF_NR],							# 12
									 [C_ROOF_NR, C_WALL_NR, C_FLOOR_NR],				# 13 								
									 [CURRENT_REF_NRES, MAX_REF_NRES],					# 14
									 [TSET, D_TSET], 									# 15
									 [BSCHED, D_BSCHED], 								# 16
									 STOCK_RES, 										# 17
									 STOCK_NRES,										# 18
									 SH_PROB]											# 19
		
		for_hot_water_demand	  =	[DHW_DEMAND, 										# 0
									 DHW_PDAYTIME, 										# 1
									 DHW_PWDAY,											# 2
									 DHW_LOAD_FLOWRATE, 								# 3
									 DHW_LOAD_DURATION]									# 4
		
		self.building_stock_stats = [for_space_heating_demand,
									 for_hot_water_demand]
#
	def read_raw_building_data(self):
		"""
		Reads building data from csv file. Returns a pd.DataFrame.
		Columns are renamed to variables in UrbanHeatPro.
		"""
		
		# Filename
		input_dir 		= self.my_dir + '/input/Buildings/'
		filename 	  	= input_dir + self.buildings_filename
				
		# Building data
		self.buildings	= pd.read_csv(filename, delimiter = ";", skiprows = 0)
		if self.debug  != 0:
			print('\nRaw building data (csv)')
			print('  ' + filename)
			print('   {} buildings\n'.format(len(self.buildings)))	
		
		# rename columns
		### Required columns
		new_columns = {"bid"  		: "bid",
					   "area" 		: "footprint_area",
					   "use"  		: "use",
					   "free_walls"	: "free_walls",
					   "lat"  		: "lat",
					   "lon"  		: "lon",
					   "dist2hp"	: "dist_to_heat_source"}
		self.buildings = self.buildings.rename(columns = new_columns)
		
		### Optional columns
		try:
			new_columns = {"year" : "year_class"}
			self.buildings = self.buildings.rename(columns = new_columns)
		except:
			pass
		#
		try:
			new_columns = {"btype" : "size_class"}
			self.buildings = self.buildings.rename(columns = new_columns)
		except:
			pass
		#
		try:
			new_columns = {"ref_roof"  	: "ref_level_roof",
						   "ref_wall"  	: "ref_level_wall",
						   "ref_floor"  : "ref_level_floor",
						   "ref_window"	: "ref_level_window",}
			self.buildings = self.buildings.rename(columns = new_columns)
		except:
			pass
		#
		try:
			new_columns = {"f" : "floors"}
			self.buildings = self.buildings.rename(columns = new_columns)
		except:
			pass
		#
		try:
			new_columns = {"d" : "dwellings"}
			self.buildings = self.buildings.rename(columns = new_columns)
		except:
			pass
		#
		try:
			new_columns = {"occ" : "occupants"}
			self.buildings = self.buildings.rename(columns = new_columns)
		except:
			pass
#
	def read_syn_city(self, filename):
		"""
		Reads existing syn city file
		"""
				
		# Building data
		buildings	= pd.read_csv(filename, delimiter = ";", skiprows = 1)
		if self.debug  != 0:
			print('\nSynthetic city (csv)')
			print('  ' + filename)
			print('   {} buildings\n'.format(len(buildings)))
			
		return buildings
#
	def read_Tamb(self):
		"""
		Reads Tamb data from csv file. The file contains the Tamb values for the whole year
		in simulation resolution. Only the simulation timesteps are extracted.
		"""
		
		filename  = '{}/input/Regional Data/{}/Tamb_{}.csv'.format(self.my_dir, self.region, self.region) 
		self.Tamb = self.read_data_from_csv(filename, usecols = 0)		
#
	def read_I(self):
		"""
		Reads solar radiation data from csv file. The file contains I values [W/m2] for the whole year
		in simulation resolution in the form [I_Gh, I_Dh, I_ex, hs]. Only the simulation timesteps are extracted.
		"""
		
		filename  = '{}/input/Regional Data/{}/I_{}.csv'.format(self.my_dir, self.region, self.region) 
		self.I = self.read_data_from_csv(filename, usecols = range(0, 4))		
#
	def filter_weather_data(self):
		"""
		Filter weather data with timesteps vector with typical days
		"""
		self.Tamb = self.Tamb[self.timesteps]
		self.I = self.I[self.timesteps]
#
	def update_Tamb(self):
		"""
		Updates Tamb of City object according to the scenario to simulate
		"""

		if self.sce_refurbishment:
			self.sce = '{}_{}'.format(self.sce_refurbishment, self.sce_Tamb)
		else:
			self.sce = '{}'.format(self.sce_Tamb)
		filename  = '{}/input/Scenarios/{}/Tamb_{}_{}.csv'.format(self.my_dir, self.sce, self.sce_Tamb, self.region) 
		Tamb = self.read_data_from_csv(filename, usecols = 0)
		self.my_city.Tamb = Tamb[self.timesteps]
#
	def read_refurbishment_matrices(self):
		"""
		Reads refurbishment matrix for residential and non residential buildings
		for scenario simulated.
		"""
		
		self.sce = '{}'.format(self.sce_refurbishment, self.sce_Tamb)
		if self.sce_Tamb:
			self.sce += '_{}'.format(self.sce_Tamb)
		
		# Residential
		size_class = ['sfh', 'th', 'mfh', 'ab']
		filename  = '{}/input/Scenarios/{}/Refurbishment_{}_Residential.csv'.format(self.my_dir, self.sce, self.sce_refurbishment) 
		ref_matrix_res = pd.read_csv(filename, delimiter =';') # df
		ref_matrix_res = ref_matrix_res.set_index(['size_class', 'year_class']) # df
		self.ref_matrix_res = [ref_matrix_res.loc[sc, :].values for sc in size_class] # np.array
		
		# Non residential
		filename  = '{}/input/Scenarios/{}/Refurbishment_{}_NonResidential.csv'.format(self.my_dir, self.sce, self.sce_refurbishment)
		ref_matrix_nres = pd.read_csv(filename, delimiter =';') # df
		ref_matrix_nres = ref_matrix_nres.set_index('year_class') # df
		self.ref_matrix_nres = ref_matrix_nres.values # np.array
#
	def prepare_result_directory(self, include_date = True):
		"""
		Creates a time stamped directory within the result folder.
		Returns path as string.
		"""

		result_dir_ = self.my_dir + '/results' 
		if include_date:
			now = datetime.now().strftime('%d%m%Y %H%M')
			self.result_dir = os.path.join(result_dir_, '{} - {}'.format(now, self.name))
		else:
			self.result_dir = os.path.join(result_dir_, '{}'.format(self.name))
		if not os.path.exists(self.result_dir): 
			os.makedirs(self.result_dir)
			
		if self.debug  != 0:
			print('\nResults directory')
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
	def calculate_typical_days(self):
		"""
		Calculates typical days based on Tamb timeseries.
		Based on Nahmmacher et al. (2016), Carpe diem: A novel approach to select 
		representative days for long-term power system modeling.
		"""
		
		if self.debug  != 0:
			print('Typical days: {}'.format(self.number_of_typ_days))
		
		# 1. Divide data into days
		# ------------------------
		## initialize variables
		days_in_year = len(self.Tamb) // ((24 * 60) // self.resolution)
		data_in_days = np.zeros((days_in_year, (24 * 60) // self.resolution), dtype = float)
		row_start    = np.zeros((days_in_year), dtype = int)
		row_end      = np.zeros((days_in_year), dtype = int)

		## arrange data
		data = self.Tamb
		for day in range(days_in_year):
			if (day == 0):
				row_start[day]   = 0
				row_end[day]     = 24 * 60 / self.resolution
			else:
				row_start[day]   = row_end[day - 1]
				row_end[day]     = 24 * 60 / self.resolution + row_start[day]

			data_in_days[day, :] = np.reshape(data[row_start[day]:row_end[day]], (24,))

		
		# 2. Hierarchical clustering
		# --------------------------
		# Linkage matrix using Ward's method
		# Ward's algorithm groups similar days based on a dstance measure that leads to clusters with a minimum inner-cluster 
		# variance.
		Z = linkage(data_in_days, 'ward')

		
		# 3. Number of clusters
		# -----------------------
		# Number of clusters
		max_c = self.number_of_typ_days
		clusters = fcluster(Z, max_c, criterion = 'maxclust')

		# Distribution of historical days sorted by months and clusters
		## Number of clusters
		number_of_clusters = max(clusters)

		## Months
		month_names        = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','','Year'];
		days_in_month      = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
		month_cumsum       = np.cumsum(days_in_month);
							
		## Sorted clusters by temperature
		### sort clusters
		mean_clusters     = np.zeros((number_of_clusters, 24))
		max_mean_clusters = np.zeros(number_of_clusters)
		for c in range(1, number_of_clusters + 1):
			indices_in_cluster = np.where(clusters == c)[0]
			mean_clusters[c - 1, :] = np.mean(data_in_days[indices_in_cluster, :], axis = 0)
			max_mean_clusters[c - 1] = np.max(mean_clusters[c - 1, :])

		# sorted_clusters = np.argsort(max_mean_clusters)[::-1] # Hot days first
		sorted_clusters = np.argsort(max_mean_clusters)			# Cold days first

		### substitute ordered clusters
		clusters_copy = copy.deepcopy(clusters)	
		for iii, c in enumerate(sorted_clusters):
			indices_in_cluster = np.where(clusters_copy == c + 1)[0]
			clusters[indices_in_cluster] = iii + 1
			
		## Annual share per cluster
		clusters_per_month = np.zeros((number_of_clusters, 14))
		for c in range(number_of_clusters):
			clusters_per_month[c, 13] = np.sum(clusters == c + 1) / days_in_year
			
			for m in range(12): # month
				# group per month
				if m == 0:
					c_start = 0
				else:
					c_start = month_cumsum[m - 1]
				c_end = month_cumsum[m]
				days_per_month = np.zeros(days_in_month[m])
				days_per_month = clusters[c_start:c_end]  
				clusters_per_month[c, m] = np.sum(days_per_month == c + 1) / days_in_month[m]
				
		# 4. Deriving representative days
		# -------------------------------
		# All historical days grouped into the same cluster will be repersented by the 
		# same representative day. The representative day for each cluster is the 
		# historical day that with the min distance to the average day.

		min_distance_index_c = np.zeros(number_of_clusters, dtype = int)
		min_distance_index_y = np.zeros(number_of_clusters, dtype = int)
		min_distance_day     = np.zeros((number_of_clusters, 24))
		avg_day              = np.zeros((number_of_clusters, 24))

		for c in range(1, number_of_clusters + 1):
			indices_in_cluster = np.where(clusters == c)[0]
			
			if len(indices_in_cluster) > 1:
						
				# Representative day
				## append average and historical days
				avg_day[c - 1, :]  = np.mean(data_in_days[indices_in_cluster, :], axis = 0)
				historical_days    = data_in_days[indices_in_cluster, :]
				data_to_compare    = np.append(np.array([avg_day[c - 1, :]]), historical_days, axis = 0)
				
				## calculate distance matrix (eudclidean)
				distance_to_avg    = squareform(pdist(data_to_compare, 'euclidean'))[0, :]
				
				## identify day with minimum distance to average day
				min_distance_index_c[c - 1] = np.where(distance_to_avg == np.min(distance_to_avg[distance_to_avg > 0]))[0][0]
				min_distance_index_y[c - 1] = indices_in_cluster[min_distance_index_c[c - 1] - 1]
				min_distance_day[c - 1, :]  = historical_days[min_distance_index_c[c - 1] - 1, :]
			
			# if there is only one day in cluster, then this is the representative day
			else:
				avg_day[c - 1, :]           = data_in_days[indices_in_cluster, :]
				min_distance_day[c - 1, :]  = data_in_days[indices_in_cluster, :]
			
					
		# 5. Weighting representative days
		# --------------------------------
		# Representative days are weighten according to its cluster size to reflect 
		# the fact that large clusters, containing many days, represent common events, 
		# while smaller clusters represent occasional patterns.
		clusters_weight = np.zeros(number_of_clusters, dtype = int)
		for c in range(number_of_clusters):
			clusters_weight[c] = len(np.where(clusters == c + 1)[0])

			
		# 6. Deriving time series with representative days
		# ------------------------------------------------
		# A time series with the defined number of clusters is created and used for 
		# the simulation. 
		# Advantages: Less computational time, intraday behavior well represented 
		# (thermal storage during the same day)
		# Disadvantages: Interday day behavior ist represented as representative days 
		# are not consecutive and therefore do not represent the storage between days due to thermal capacity and/or solar gains.
		## Time series
		timeseries_min = min_distance_day.flatten()
		timeseries_avg = avg_day.flatten()
		
		## Time steps
		timesteps  = np.array([np.arange(row_start[d], row_end[d]) for d in min_distance_index_y]).flatten()
		
		if self.save >= 1:
			UrbanHeatPro.plot_typical_days(days_in_year, data_in_days, \
					  Z, self.number_of_typ_days, \
					  min_distance_day, avg_day, clusters, \
					  clusters_per_month, month_names, \
					  timeseries_min, timeseries_avg, self.result_dir)
			if self.debug  != 0:
				print('  {}\TypicalDays'.format(self.result_dir))
				
		return timesteps, clusters_weight	
#
	def calculate_dt_vector(self):
		"""
		Calculates a vector of datetime objects based on the raw dt_matrix of the
		form [Y, M, D, h, m] and the simulation time steps.
				
		returns:
			self.dt_vector  <list>	List of datetime objects	
		"""
		
		# read csv file
		filename   	   = '{}/input/Building Typology/{}'.format(self.my_dir, 'YMdhm.csv') 
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