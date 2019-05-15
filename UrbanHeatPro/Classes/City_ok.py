"""
City.py
AMC @ TUM ENS
"""

import numpy as np
import pandas as pd
from datetime import datetime
from scipy import interpolate
import multiprocessing

import UrbanHeatPro.Functions as UrbanHeatPro
from .Building import Building

class City():
# --------------------------------------------------------------------------------
	def __init__(self, SIMULATION, CITY, SPACE_HEATING, HOT_WATER, REPORTING):
		
		# SIMULATION
		self.region				  = SIMULATION[0][0]
		self.dt_vector    		  = SIMULATION[1][0] 		# Vector of time steps as datetime objects
		self.dt_vector_excel	  = SIMULATION[1][1] 		# Vector of time steps as excel date
		self.nts          		  = len(self.dt_vector)		# Number of time steps
		self.resolution			  = SIMULATION[2][0]		# Temporal resolution in min
		self.processes			  = SIMULATION[3][0]		# Number of parallel processes
		
		# CITY
		### Ambient conditions
		self.Tamb         		  = CITY[0][0]		 		# Ambient temperature vector in degC
		self.I	         	  	  = CITY[0][1]				# Solar radiation vector in W/m2 [I_Gh, I_Dh, I_ex, hs]
		### Building data
		self.buildings 			  = CITY[1][0]				# Building data
		self.building_stock_stats = CITY[1][1]				# Building data from statistics
		self.nb					  = len(self.buildings)		# Number of buildings
		self.building_typology	  = CITY[1][2]				# Building typology given? (as boolean)
		self.connection_factor	  = CITY[2][0]				# Share of buildings connected to the network
		### Flags
		self._space_heating		  = CITY[3][0]				# calculate space heating demand?
		self._hot_water		 	  = CITY[3][1]				# calculate hot water demand?
		### Base load
		self.base_load			  = np.ones([self.nts]) * CITY[4][0]	# Vector of base load in W
		
		# SPACE HEATING DEMAND
		### Flags
		self._internal_gains	  = SPACE_HEATING[0][0]		# consider internal gains?
		self._solar_gains		  = SPACE_HEATING[0][1]		# consider solar gains?
		self._active_population	  = SPACE_HEATING[0][2]		# consider active pop for occupancy vector?	
		### Transmission losses
		self.refurbishment_level  = SPACE_HEATING[1][0]		# Refurbishment level for all buildings
		self.Tb0_str			  = SPACE_HEATING[2][0]		# Initial building temperature as 'ambient' or 'Tset'?
		self.dTset			 	  = SPACE_HEATING[2][1]		# Delta temperature (for Tset_min, Tset_max)
		### Heating system
		self.eta				  = SPACE_HEATING[3][0]		# Heating process efficiency
		self.dT_per_hour		  = SPACE_HEATING[3][1]		# Maximum dT allowed in building per hour [degC]
		self.thermal_inertia	  = SPACE_HEATING[3][2]		# Thermal inertia of the heating system
		### DSM
		self._night_set_back  	  = SPACE_HEATING[4][0]		# Share of buildings with night set-back
		self.schedule_nsb		  = SPACE_HEATING[4][1]		# [start, end] of nsb in h
		self.T_nsb		  		  = SPACE_HEATING[4][2]		# Night set-back temperature in degC
		self.power_reduction 	  = SPACE_HEATING[4][3]		# Percentage of power reduced (as decimal)
		
		
		# HOT WATER DEMAND
		### Calculation
		self.Tw					  = HOT_WATER[0][0]			# Hot water temperature in degC
		self.hw_tank_limit	   	  = HOT_WATER[1][0]			# Hot water tank limit as perc (decimal)
		self.hw_flow		      = HOT_WATER[1][1]			# Flow to refill hot water tank in L/min
		self.dhw_prob			  = self.initialize_dhw_probabilities()	# Probabilities for calculation of DHW
		### Vectors
		self.day_vector			  = self.calculate_day_vector()	# Vector of days in time frame
		self.seasonal_vector	  = self.calculate_seasonal_variation_vector() # Sine vector to represent seasonal variations
		self.min_vector			  = self.calculate_min_vector() # Vector of minutes in time frame
		
		# RESULTS
		### Runs
		self.rid				  = SIMULATION[4][0]		# Run id
		self.result_dir			  = SIMULATION[4][1]		# Directory where results are stored
		### Space heating demand
		self.space_heating_power  = np.zeros([self.nts])    # City Space heating demand in W
		self.space_heating_energy = 0.0						# Aggregated city heating demand in Wh
		### Hot water demand
		self.hot_water_power      = np.zeros([self.nts])	# City hot water demand in W
		self.hot_water_energy	  = 0.0						# Aggregated city hot water demand in Wh
		### Total heating energy demand
		self.total_power   		  = np.zeros([self.nts])	# City total heat demand in W
		self.total_energy	  	  = 0.0						# Aggregated city total heat demand in Wh
		self.total_power_delayed  = np.zeros([self.nts])	# City total delayed heat demand in W
		self.total_energy_delayed = 0.0						# Aggregated city total delayed heat demand in Wh
		self.energy_per_building  = np.zeros([self.nb, 32]) # Matrix of energy demand per building
		
		# REPORTING
		self.plot				  = REPORTING[0]					# Plot level  [0, 1, 2]
		self.save 				  = REPORTING[1]					# Save level  [0, 1, 2]
		self.debug		   		  = REPORTING[2]					# Debug level [0, 1, 2]

		
# Synthetic city
# --------------------------------------------------------------------------------
	def create_synthetic_city(self):
		"""
		Create a synthetic city representing the building stock based on statistics. 			
		"""
		
		filename = '{}/SyntheticCity_{}.csv'.format(self.result_dir, self.rid)
		
		# extract building stock statistics
		self.FOOTPRINT			  = self.building_stock_stats[0][0]  # Footprint area of building typologies (TABULA)
		self.FLOORS				  = self.building_stock_stats[0][3]  # Number of floors
		self.REFURBISHED_RES	  = self.building_stock_stats[0][7]  # Percentage of residential refurbished buildings
		self.SINGLE_DWELLING	  = self.building_stock_stats[0][8]  # Percentage of single dwellings for SFH and TH
		self.AVG_DWELLING_SIZE	  = self.building_stock_stats[0][9]  # Average dwelling size in m2
		self.HOUSEHOLD_SIZE		  = self.building_stock_stats[0][10] # Household size for dwelling size categories
		self.STOCK_RES			  = self.building_stock_stats[0][16] # Building stock statistics for residential
		self.STOCK_NRES			  = self.building_stock_stats[0][17] # Building stock statistics for non-residential
					
		# Multiprocessing
		if (self.debug >= 1):
			print('      ***')
			print('      Starting multiprocessing with {}/{} threads...'.format(self.processes, multiprocessing.cpu_count()))
			
		### Writer queue
		if (self.debug >= 1):
			print('      Starting queue...')
		writerQueue = multiprocessing.Queue()
		writeProc   = multiprocessing.Process(target = self.write_to_synthetic_city, args = (writerQueue, filename))
		
		### Feeder queue
		buildings_list = range(len(self.buildings))
		feederQueue = multiprocessing.Queue()
		feedProc    = multiprocessing.Process(target = self.feed_building_to_process, args = (feederQueue, buildings_list))
		
		### Processes
		if (self.debug >= 1):
			print('      Starting processes...')
		calcProc = [multiprocessing.Process(target = self.call_create_synthetic_building,
					args = (feederQueue, writerQueue)) for iii in range(self.processes)]
		
		### start multiprocessing
		for p in calcProc:
			p.start()
		writeProc.start()
		
		for p in calcProc:
			p.join()
		writeProc.join()
#
	def feed_building_to_process(self, buildings_list, workerQueue):
		for iii in buildings_list:
			building = self.buildings.loc[iii, :]
			workerQueue.put(building)
#
	def call_create_synthetic_building(self, feederQueue, writerQueue):
		while True:
			try:
				building = feederQueue.get()
				result = self.create_synthetic_building(building)
				writerQueue.put(result)
				print('         {} | {}'.format(multiprocessing.current_process().name, result[0]))
			except:
				break
#
	def write_to_synthetic_city(self, writerQueue, filename):
		f = open(filename, 'a')
		while True:
			try:
				result = writerQueue.get(block = False)
				print(result)
				np.savetxt(f, result, delimiter = ';')
			except:
				break
		f.close()
#
	def create_synthetic_building(self, building):
		"""
		Calculate building properties based on statistics.
		It first checks if the property is already given in the input file.
		"""
		
		# RESIDENTIAL
		if (building.use == 3):
			
			# Construction year class and size class
			try:
				year_class = building.year_class
				size_class = building.size_class
			except:
				# categorize building according TABULA typology
				[year_class, size_class]   = self.categorize_building_residential(building.footprint_area, building.use)
			
			# Refurbishment level per element
			try: 
				ref_level_roof   = building.ref_level_roof
				ref_level_wall   = building.ref_level_wall 
				ref_level_floor  = building.ref_level_floor 
				ref_level_window = building.ref_level_window
			except:
				# compute refurbishment level
				[ref_level_roof, 
				 ref_level_wall, 
				 ref_level_floor, 
				 ref_level_window] = self.compute_refurbishment_level_residential(year_class, size_class)
			
			# Number of floors
			try:
				floors = building.floors
			except:
				floors = False
			# calculate floors and reference areas
			[area_correction_factor, floors, 
			storey_area, heated_area] = self.calculate_areas_residential(building.footprint_area, 
																		 year_class, size_class, floors)
			
			# Number of dwellings
			try:
				dwellings = building.dwellings
			except:
				dwellings = False
			# calculate number of dwellings per building
			dwellings, dwelling_size = self.calculate_number_of_dwellings(year_class, size_class, heated_area, dwellings)
			
			# Number of occupants
			try:
				occupants = building.occupants
			except:
				# calculate number of occupants per building
				occupants = self.calculate_number_of_occupants_residential(dwellings, dwelling_size)
			
		# NON-RESIDENTIAL
		else:
		
			# determine construction year class and size class
			try:
				year_class = building.year_class
				size_class = building.size_class
			except:
				# categorize building 
				[year_class, size_class] = self.categorize_building_non_residential()
			
			# Refurbishment level per element
			try: 
				ref_level_roof   = building.ref_level_roof
				ref_level_wall   = building.ref_level_wall 
				ref_level_floor  = building.ref_level_floor 
				ref_level_window = building.ref_level_window
			except:
				# compute refurbishment level
				[ref_level_roof, 
				 ref_level_wall, 
				 ref_level_floor, 
				 ref_level_window] = self.compute_refurbishment_level_non_residential()

			# Floors
			try:
				floors = building.floors
			except:
				floors = False
			# calculate floors and reference areas
			[area_correction_factor, floors, 
			storey_area, heated_area] = self.calculate_areas_non_residential(building.footprint_area, floors)
			
			# Number of occupants
			try:
				occupants = building.occupants
			except:
				# calculate number of occupants per building
				occupants     = self.calculate_number_of_occupants_non_residential(building.use, heated_area)
			dwellings     = -1
			dwelling_size = -1.
			
		# Distance to heat source
		try:
			dist_to_heat_source = building.dist_to_heat_source
		except:
			dist_to_heat_source = -1.
		
		# update synthetic city
		#self.append_building_to_synthetic_city(iii, building, dist_to_heat_source,
		#									year_class, size_class, 
		#									ref_level_roof, ref_level_wall, ref_level_floor, ref_level_window,
		#									floors, area_correction_factor, storey_area, heated_area,
		#									dwellings, dwelling_size, occupants)
		
		
		result = [building.bid, 
						year_class, size_class, 
						ref_level_roof, ref_level_wall, ref_level_floor, ref_level_window]
		result = np.resize(result, (1, len(result)))						
		return result
#
	def categorize_building_residential(self, footprint_area, use):
		"""
		Probabilistic categorization building according to TABULA typologies. 
		Construction year class and building type are calculated by comparing the residential
		building gross floor area (footprint_area) with the FOOTPRINT of typical buildings (from TABULA).
		Values are adapted to fit building stock statistics.
		
		returns:
			self.class_year	  <tuple>	(int, str)
			size_class		  <tuple>	(int, str)
		"""
		
		# get size of FOOTPRINT matrix
		rows 			= len(self.FOOTPRINT)	# construction year class
		cols 			= len(self.FOOTPRINT[0])	# building type
		
		# initialize distance vectors
		distance        = np.zeros(rows * cols)
		distance_inv    = np.zeros(rows * cols)				# 1/d
		row_col 	    = [[] for _ in range(rows * cols)]
		t_distance_inv  = 0									# total inv distance, sum(1/d)
		kkk 		    = 0									# typology counter
		
		
		for iii in range(0, rows):
			for jjj in range(0, cols):
				# calculate distance between building footprint_area and FOOTPRINT (TABULA) matrix
				distance[kkk]     = abs(self.FOOTPRINT[iii][jjj] - footprint_area) + 0.01 # "+ 0.001" to avoid having x/0
				# calculate inverse and weight by building stock statistics
				distance_inv[kkk] = (1.0 / distance[kkk]) * self.STOCK_RES[iii, jjj]
				row_col[kkk] 	  = [iii, jjj]
				t_distance_inv   += distance_inv[kkk] 
				kkk              += 1
		
		# normalize inverse distance vector
		norm_distance   = distance_inv / t_distance_inv

		# calculate cumulative density function
		cdf             = np.cumsum(norm_distance)

		# Probabilistic categorization
		rnd       	    = np.random.uniform(0, 1, 1)
		index 			= np.argmax(rnd < cdf)
		### Construction year class as int
		year_class 		= UrbanHeatPro.year_class_to_tuple(use, row_col[index][0])[0]
		### Building size class as int
		size_class 		= UrbanHeatPro.size_class_to_tuple(use, row_col[index][1])[0]
		
		return [year_class, size_class]
#	
	def categorize_building_non_residential(self):
		"""
		Probabilistic categorization of non-residential buildings according to the
		following construction year classes: 
		
		int	construction year class
		 0	 < 1918
		 1	 1919 - 1976
		 2	 1977 - 1983
		 3	 1984 - 1994
		 4	 > 1995

		>>> Statistics on the non-residential buildings stock are missing
		"""
		
		# assign random construction year class
		year_class = np.random.randint(0, 4)
		
		# set size class to Nan
		size_class	= -1

		return [year_class, size_class]
#
	def compute_refurbishment_level_residential(self, year_class, size_class):
		"""
		Determines the refurbishment level for the different building elements 
		[roof, wall, floor, window] depending on the refurbishment statistics or
		scenarios and the building typology.
		Refurbishment levels according to TABULA typology:
			1	National minimum requirement
			2	Improved standard
			3	Ambitious standard
		"""
		
		# check building size class
		if (size_class == 0 or size_class == 1): 		# SFH or TH
			# check building construction year class
			if (year_class <= 5): 						# <1978
				perc_refurbished = self.REFURBISHED_RES[0]
			elif (year_class > 5 and year_class <= 7): 	# 1978-1994
				perc_refurbished = self.REFURBISHED_RES[1]
			elif (year_class > 7 and year_class <= 9): 	# >1994
				perc_refurbished = self.REFURBISHED_RES[2]
		else: 											# MFH or AB
			# check building construction year class
			if year_class <= 5: 						# <1978
				perc_refurbished = self.REFURBISHED_RES[3]
			elif year_class > 5 and year_class <= 7: 	# 1978-1994
				perc_refurbished = self.REFURBISHED_RES[4]
			elif year_class > 7 and year_class <= 9: 	# >1994
				perc_refurbished = self.REFURBISHED_RES[5]
				
		# assign refurbishment level for every building element
		# if random number is smaller than percentage, then refurfishment level of
		# element is 2, otherwise is 1
		refurbishment_level = np.ones(4, dtype = int)
		for element in range(4): # [roof, wall, floor, window]	
			rand_num = np.random.uniform(0, 1, 1)
			if rand_num <= perc_refurbished[element]:
				refurbishment_level[element] = 2
				
		return refurbishment_level
#
	def compute_refurbishment_level_non_residential(self):
		"""
		Determines the refurbishment level for the different building components
		[roof, wall, floor, window] depending on the refurbishment statistics or
		scenarios and the building typology.
		Refurbishment levels according to TABULA typology:
			1	National minimum requirement
			2	Improved standard
			3	Ambitious standard
			
		Refurbishment level in non residential buildings is assumed to be 1
		>>> Statistics on refurbishment in non-residential buildings
		"""
		
		# all components are assumed to have refurbishment level 1.
		return [1, 1, 1, 1]
#
	def calculate_areas_residential(self, footprint_area, year_class, size_class, floors):
		"""
		Calculate storey area and heated/conditioned area based on definitions
		from VDI 3807.
		"""
			
		# Area correction factor (VDI 3807-2 Section 7.1.1 Table 1)
		area_correction_factor = 0.84
		
		# Floors
		if not floors:
			floors = self.calculate_number_of_floors_residential(year_class, size_class, area_correction_factor)
		
		# Storey area (VDI 3807-1 Section 5.4.2)
		storey_area = footprint_area * floors
		
		# Heated area (VDI 3807-1 Section 5.4.2)
		heated_area = footprint_area * area_correction_factor * floors
		
		return [area_correction_factor, floors, storey_area, heated_area]
#
	def calculate_areas_non_residential(self, footprint_area, floors):
		"""
		Calculate storey area and heated/conditioned area based on definitions
		from VDI 3807.
		"""
			
		# Area correction factor (VDI 3807-2 Section 7.1.1 Table 1)
		area_correction_factor = np.random.randint(75, 85) / 100.
		
		# Floors
		if not floors:
			floors = self.calculate_number_of_floors_non_residential()
		
		# Storey area (VDI 3807-1 Section 5.4.2)
		storey_area = footprint_area * floors
		
		# Heated area (VDI 3807-1 Section 5.4.2)
		heated_area = footprint_area * area_correction_factor * floors
		
		return [area_correction_factor, floors, storey_area, heated_area]
#	
	def calculate_number_of_floors_residential(self, year_class, size_class, area_correction_factor):
		"""
		Calculates number of floors based on the TABULA typology.
		The number of floors calculated from TABULA are referenced to the conditioned
		or heated area but the number of floors are calculated using the storey area.
		"""
		
		# get number of floors
		floors 	= self.FLOORS[year_class][size_class] / (1 + (1 - area_correction_factor))
		
		return floors
#
	def calculate_number_of_floors_non_residential(self, left = 1, mode = 2, right = 3):
		"""
		Calculates number of floors as random sample number from the triangular
		distribution with lower limit left, peak at mode and upper limit right.
		
		>>> Non-residential buildings are assumed to have two floors as mode and
			a maximum of three floors
			Source missing
		"""

		# calculate number of floors
		floors = np.random.triangular(left, mode, right)
		
		return floors
#
	def calculate_number_of_dwellings(self, year_class, size_class, heated_area, dwellings):
		"""
		Calculates number of dwellings based on the building living area and mean dwelling size.
		It is assumed that SFH and TH have only 1 or 2 dwellings which is determined using the
		single-dwelling buildings statistics. For MFH and AB, the number of dweelings is calculated
		based on the average dwelling size.
		"""
				
		if not dwellings:
			# For SFH and TH
			# compare random number to single-dwelling buildings statistics
			if (size_class == 0 or size_class == 1): # SFH or TH
				perc_single_dwelling = self.SINGLE_DWELLING[year_class][size_class]
				rand_num = np.random.uniform(0, 1, 1)
				if rand_num <= perc_single_dwelling:
					dwellings = 1
				else:
					dwellings = 2
				
			# For MFH and AB
			# divide total living area by dwelling size
			else:
				# calculate random dwelling size from normal dist with mean AVG_DWELLING_SIZE
				# and sigma AVG_DWELLING_SIZE/2
				dwelling_size = 0
				while (dwelling_size <= 25): # min dwelling size [m2]
					dwelling_size = np.random.normal(self.AVG_DWELLING_SIZE, self.AVG_DWELLING_SIZE / 10, 1)
				# calculate number of dwellings
				dwellings = int(np.ceil(heated_area / dwelling_size))
		
		# Dwelling size
		dwelling_size = heated_area / dwellings
			
		return [dwellings, dwelling_size]
#
	def determine_dwelling_size_category(self, dwelling_size):
		"""
		Determine dwelling size category based on statistics 
		https://ergebnisse.zensus2011.de/#StaticContent:091840148148,GWZ_4_3_2,m,table
		"""

		# determine size category
		if dwelling_size > 0 and dwelling_size <= 40:
			dwelling_size_cat = 0
		elif dwelling_size > 40 and dwelling_size <= 60:
			dwelling_size_cat = 1
		elif dwelling_size > 60 and dwelling_size <= 80:
			dwelling_size_cat = 2
		elif dwelling_size > 80 and dwelling_size <= 100:
			dwelling_size_cat = 3
		elif dwelling_size > 100 and dwelling_size <= 120:
			dwelling_size_cat = 4
		elif dwelling_size > 120 and dwelling_size <= 140:
			dwelling_size_cat = 5
		elif dwelling_size > 140 and dwelling_size <= 160:
			dwelling_size_cat = 6
		elif dwelling_size > 160 and dwelling_size <= 180:
			dwelling_size_cat = 7
		elif dwelling_size > 180 and dwelling_size <= 200:
			dwelling_size_cat = 8
		else:
			dwelling_size_cat = 9
		
		return dwelling_size_cat
#
	def calculate_number_of_occupants_residential(self, dwellings, dwelling_size):
		"""
		Calculates number of occupants based on household size and number of dwellings 
		statistics.
		"""
		
		# calculate dwelling size category
		dwelling_size_cat = self.determine_dwelling_size_category(dwelling_size)
		
		# CFD for household size
		### x-values
		household_size   = range(1, 7)
		### CFD
		cfd = np.cumsum(self.HOUSEHOLD_SIZE[dwelling_size_cat])
		
		# calculate number of occupants
		occupants_in_building = 0
		for dwelling in range(dwellings):		
			# get random_number and obtain x-value from cfd
			rand_num = np.random.uniform(0, 1, 1)
			index = np.argmax(rand_num <= cfd)
			# get household size (number of occupants)
			occupants_in_dwelling = household_size[index]
					
			# add occupants to building
			occupants_in_building += occupants_in_dwelling
		
		return occupants_in_building
#
	def calculate_number_of_occupants_non_residential(self, use, heated_area):
		"""
		Calculates random number of occupants in the building based on the 
		recommended area per person for different building types from
		https://www.engineeringtoolbox.com/number-persons-buildings-d_118.html.
		Buildings are assumed to be always occupied at their max capacity.
		"""
				
		if   (use == 0): # commercial
			# retail, supermarket, department stores
			area_per_person = np.random.randint(3, 10) 	# in m2
			
		elif (use == 1): # industrial
			# light manufacturing, heavy manufacturing
			area_per_person = np.random.randint(10, 30)	# in m2
		
		elif (use == 2): # public
			# municipal buildings, library, museum
			area_per_person = np.random.randint(3, 10)	# in m2
			
		# calculate number of occupants
		occupants_in_building = heated_area / area_per_person
		
		return occupants_in_building	

		
# City heat demand
# --------------------------------------------------------------------------------		
	def calculate_city_heat_demand(self):
		"""
		Paralellizes the calculation of heating energy demand per building using a 
		given number of processes. Every process modifies a shared dictionary where
		the heat demand is stored as power and energy.
		"""
		
		if (self.debug >= 1):
			print('\n***\nMultiprocessing using {}/{} available processors'.format(self.processes, multiprocessing.cpu_count()))
								
		# save dhw debug matrix
		save_hwd_debug = False
		if (self.save == 3):
			save_hwd_debug = True
					
		# define multiprocessing pool
		if (self.debug >= 1):
			print('   Starting multiprocessing pool...')
		pool    = multiprocessing.Pool(self.processes)
		
		# define city dictionary as multiprocessing manager
		# this dictionary is shared with all the processes
		if (self.debug >= 1):
			print('   Starting multiprocessing manager...')
		manager = multiprocessing.Manager()
		self.results = manager.dict(power_sh 	   = np.zeros([self.nts]),
									  energy_sh    = 0.,
									  power_hw 	   = np.zeros([self.nts]),
									  energy_hw    = 0.,
									  power_total  = np.zeros([self.nts]),
									  energy_total = 0.,
									  buildings	   = np.zeros([self.nb, 32]))

		# set tasks for multiprocessing
		# Task: calculate heat demand for every building in the city
		if (self.debug >= 1):
			print('   Setting tasks...')
		TASKS   = [(self.calculate_building_heat_demand, 
					([self.buildings[iii, :], iii, save_hwd_debug])) for iii in range(len(self.buildings))]
		
		if (self.debug >= 1):
			print('   Multiprocessing started correctly.\n***\n')
		
		# assign tasks to processes
		results = pool.imap_unordered(self.calculatestar, TASKS)
		for r in results:
			if (self.debug >= 1):
				print('   ', r)
				
		# save results per city
		if (self.save >= 1):
			self.save_csv_power()
			self.save_csv_energy()
#						
	def calculate(self, func, args):		# for Multiprocessing
		result = func(*args)
		return '%s | %s' % (
			multiprocessing.current_process().name, result)
#
	def calculatestar(self, args):	 		# for Multiprocessing
		return self.calculate(*args)

		
# Building heat demand
# --------------------------------------------------------------------------------	
	def calculate_building_heat_demand(self, building, iii, save_hwd_debug):
		"""
		
		Extracts building information needed to create a Building object. 
		If the building is connected to the district heating network, then a Building
		object is created and the heat demand is calculted. If it is not, then the
		heat demand is set to zero.
		
		args:
			building	dataframe with building information
			iii			building counter
		"""
		
		# extract building data
		building_data = self.extract_building_data(building)
		
		# check if building is connected to the network
		rand_num = np.random.uniform(0, 1, 1)[0]
		if (rand_num <= self.connection_factor):
	
			my_building = self.create_building_object(building_data)
	
			# calculate space heating demand
			if self._space_heating:
				my_building.calculate_space_heating_demand()
				# add building space heating demand to city space heating demand
				self.results['power_sh']  += my_building.space_heating_power
				self.results['energy_sh'] += my_building.space_heating_energy
			
			# calculate hot water demand
			if self._hot_water:
				my_building.calculate_hot_water_demand(save_hwd_debug)
				# add building hot water demand to city hot water demand
				self.results['power_hw']  += my_building.hot_water_power
				self.results['energy_hw'] += my_building.hot_water_energy
			
			# calculate total heat demand
			my_building.calculate_total_heat_demand()
			
			# add building demand to city demand
			self.results['power_total']   += my_building.total_power
			self.results['energy_total']  += my_building.total_energy
			
			# Save				
			### City results
			if (self.save >= 1):
				self.append_building_results_to_city_array(my_building, iii)

			### save results per building
			if (self.save >= 2):
				my_building.save_csv()
				#self.my_building.save_load_duration_curve()		
			### save results per timestep
			if (self.save == 3):
				if hot_water:
					my_building.save_dhw_debug_csv()
			
			# plot results per building
			if (self.plot == 2):
				my_building.plot_timeseries(space_heating = self._space_heating, Tb = True, \
										hot_water = self._hot_water, total = True)
				
		return building_data[0]
#																			
	def extract_building_data(self, building):
		"""
		Extracts and order building data from dataframe to create Building object
		"""
		
		bid 	   	   = int(building[0])  		# building id
		footprint_area = building[1]  	  		# in m2
		use 	   	   = int(building[2])  		# as int
		free_walls 	   = building[3]	  		# number of walls in contact with ambient temperature
		lat			   = building[4]	 		# latitude in degrees
		lon			   = building[5]	   		# longitude in degrees
		distance2hp	   = building[6]	   		# in m
		if self.building_typology:
			year_class_int = int(building[7])   # construction year class (0, 9)
			btype_int	   = int(building[8])   # building type (0, 3)
		else:
			year_class_int = None
			btype_int 	   = None
		
		return [bid, use, footprint_area, free_walls, lat, lon, distance2hp, year_class_int, btype_int]
#
	def create_building_object(self, building_data):
		"""
		Creates instance of class Object
		"""
		
		bid			   = building_data[0]
		use			   = building_data[1]
		footprint_area = building_data[2]
		free_walls	   = building_data[3]
		lat			   = building_data[4]
		lon			   = building_data[5]
		distance2hp	   = building_data[6]
		year_class_int = building_data[7]
		btype_int	   = building_data[8]
			
		# create Building object
		my_building = Building([self.dt_vector, self.dt_vector_excel], self.resolution, self.Tamb, self.I, \
						self._space_heating, self._hot_water, \
						self.Tb0_str, self.dTset, self.dT_per_hour, self.eta, self.thermal_inertia, self.building_stock_stats, \
						bid, use, footprint_area, free_walls, lat, lon, distance2hp, self.building_typology, year_class_int, btype_int, \
						self.refurbishment_level, self._active_population, self._solar_gains, self._internal_gains, \
						self._night_set_back, self.schedule_nsb, self.T_nsb, self.power_reduction, \
						self.Tw, self.dhw_prob, self.hw_tank_limit, self.hw_flow, \
						self.day_vector, self.seasonal_vector, self.min_vector, \
						self.result_dir_buildings, self.plot, self.save, self.debug)
							
		return my_building
#				
	def append_building_results_to_city_array(self, my_building, iii):
		"""
		Appends all building simulation data including thermal properties, 
		synthetic population and heat gains and losses to the city matrix.
		
		args:
			my_building		Building object with results
			iii				building counter
		"""
		
		# retrieve the shared dictionary and modify it
		b = self.results['buildings']
		b[iii, :] = [
			 my_building.bid, 							# 0
			 my_building.use[0], 						# 1
			 my_building.year_class[0], 				# 2
			 my_building.btype[0], 						# 3
			 my_building.env_areas[2], 					# 4
			 my_building.heated_area,					# 5
			 my_building.living_area,					# 6
			 my_building.floors, 						# 7
			 my_building.free_walls,					# 8
			 np.mean(my_building.refurbishment_level),  # 9
			 my_building.Tset, 							# 10
			 my_building.Tb[0], 						# 11
			 my_building.active_hours[0][0], 			# 12
			 my_building.active_hours[0][1],			# 13
			 my_building.active_hours[1][0], 			# 14
			 my_building.active_hours[1][1],			# 15
			 my_building.occupied,						# 16
			 my_building.U, 							# 17
			 my_building.V, 							# 18
			 my_building.C, 							# 19
			 my_building.Tau / 3600.,					# 20
			 my_building.space_heating_energy, 			# 21
			 my_building.space_heating_energy_per_area, # 22						 
			 my_building.solar_gains_per_area,			# 23
			 my_building.internal_gains_per_area,		# 24
			 my_building.daily_DHW,  	 				# 25
			 sum(my_building.hot_water_m3), 			# 26
			 my_building.hot_water_energy, 				# 27
			 my_building.total_energy,					# 28
			 my_building.distance2hp,					# 29
			 my_building.delay,							# 30
			 my_building.total_energy_delayed]			# 31
		
		# set the shared dictionary
		self.results['buildings'] = b
#	
	def calculate_total_heat_demand(self):
		"""
		Calculates the time series of the total heating energy demand for a given set of 
		buildings as GeoDataFrame. Space heating demand and domestic hot water (DHW)
		are included. 
		
		args:
			_space_heating	<boolean>
			_hot_water		<boolean>
			plot			<int>		0  none
										1  per run
										2  per building and run
			save			<int>		0  none
										1  per run
										2  per building and run
			debug			<int>		0  none
										1  basic
										2  complete
		"""
		
		# save dhw debug matrix
		save_hwd_debug = False
		if (self.save == 3):
			save_hwd_debug = True
			
		# Building loop
		for iii in range(len(self.buildings)):
		
			# Building data from csv
			bid 	   	   = int(self.buildings[iii][0])  # building id
			footprint_area = self.buildings[iii][1]  	  # in m2
			use 	   	   = int(self.buildings[iii][2])  # as int
			free_walls 	   = self.buildings[iii][3]		  # number of walls in contact with ambient temperature
			lat			   = self.buildings[iii][4]	  	  # latitude in degrees
			lon			   = self.buildings[iii][5]	  	  # longitude in degrees
			distance2hp	   = self.buildings[iii][6]		  # in m
			if self.building_typology:
				year_class_int = int(self.buildings[iii][7])  # construction year class (0, 9)
				btype_int	   = int(self.buildings[iii][8])  # building type (0, 3)
			else:
				year_class_int = None
				btype_int 	   = None
			
			if (self.debug >= 1):
				print('  ' + '({}) Building {}/{}'.format(self.rid, iii + 1, len(self.buildings)))
					
			# determine if building is connected to the network
			rand_num = np.random.uniform(0, 1, 1)[0]
			if (rand_num <= self.connection_factor):
			
				# create building object
				self.create_building_object(bid, use, footprint_area, free_walls, lat, lon, distance2hp, year_class_int, btype_int)
							
				# calculate space heating demand
				if self._space_heating:
					self.my_building.calculate_space_heating_demand()
					# add building space heating demand to city space heating demand
					self.space_heating_power  += self.my_building.space_heating_power
					self.space_heating_energy += self.my_building.space_heating_energy
				
				# calculate hot water demand
				if self._hot_water:
					self.my_building.calculate_hot_water_demand(save_hwd_debug)
					# add building hot water demand to city hot water demand
					self.hot_water_power  	  += self.my_building.hot_water_power
					self.hot_water_energy 	  += self.my_building.hot_water_energy
				
				# calculate total heat demand
				self.my_building.calculate_total_heat_demand()
			
				# add building energy results to city matrix
				# all in units W, s, J
				self.energy_per_building[iii] = [self.my_building.bid, 		# 0
							 self.my_building.use[0], 						# 1
							 self.my_building.year_class[0], 				# 2
							 self.my_building.btype[0], 					# 3
							 self.my_building.env_areas[2], 				# 4
							 self.my_building.heated_area,					# 5
							 self.my_building.living_area,					# 6
							 self.my_building.floors, 						# 7
							 self.my_building.free_walls,					# 8
							 np.mean(self.my_building.refurbishment_level), # 9
							 self.my_building.Tset, 						# 10
							 self.my_building.Tb[0], 						# 11
							 self.my_building.active_hours[0][0], 			# 12
							 self.my_building.active_hours[0][1],			# 13
							 self.my_building.active_hours[1][0], 			# 14
							 self.my_building.active_hours[1][1],			# 15
							 self.my_building.occupied,						# 16
							 self.my_building.U, 							# 17
							 self.my_building.V, 							# 18
							 self.my_building.C, 							# 19
							 self.my_building.Tau / 3600.,					# 20
							 self.my_building.space_heating_energy, 		# 21
							 self.my_building.space_heating_energy_per_area,# 22						 
							 self.my_building.solar_gains_per_area,			# 23
							 self.my_building.internal_gains_per_area,		# 24
							 self.my_building.daily_DHW,  	 				# 25
							 sum(self.my_building.hot_water_m3), 			# 26
							 self.my_building.hot_water_energy, 			# 27
							 self.my_building.total_energy,					# 28
							 self.my_building.distance2hp,					# 29
							 self.my_building.delay,						# 30
							 self.my_building.total_energy_delayed]			# 31
										
				if (self.debug == 2):
					print('      ' + 'RESULTS')
					print('      ' + '==========================================================================================================')
					print('      ' + 'Building id			{}'.format(self.my_building.bid))
					print('      ' + 'Building use			{}'.format(self.my_building.use))
					print('      ' + 'Year class			{}'.format(self.my_building.year_class))
					print('      ' + 'Building type			{}'.format(self.my_building.btype))
					print('      ' + 'FOOTPRINT [m2]				{:.3f}'.format(self.my_building.env_areas[2]))
					print('      ' + 'Heated area [m2]			{:.3f}'.format(self.my_building.heated_area))
					print('      ' + 'Living area [m2]			{:.3f}'.format(self.my_building.living_area))
					print('      ' + 'No. floors			{:.3f}'.format(self.my_building.floors))
					print('      ' + 'No. free walls		 	{}'.format(self.my_building.free_walls))
					print('      ' + 'Refurbishment level	 	{}'.format(self.my_building.refurbishment_level))
					print('      ' + 'Tset [*C]				{}'.format(self.my_building.Tset))
					print('      ' + 'Tb0  [*C]				{}'.format(self.my_building.Tb[0]))
					print('      ' + 'Active hours [h]		 	{} - {}, {} - {}'.format(self.my_building.active_hours[0][0], self.my_building.active_hours[0][1], self.my_building.active_hours[1][0], self.my_building.active_hours[1][1]))
					print('      ' + 'Household vector			{}'.format(self.my_building.household_vector))
					print('      ' + 'Occupied?				{}'.format(self.my_building.occupied))
					print('      ' + 'U [kW/K]				{:.3f}'.format(self.my_building.U / 1e3))
					print('      ' + 'V [kW/K]				{:.3f}'.format(self.my_building.V / 1e3))
					print('      ' + 'C [MJ/K]				{:.3f}'.format(self.my_building.C / 1e6))
					print('      ' + 'Tau [h]				{:.3f}'.format(self.my_building.Tau / 3600.0))
					print('      ' + 'Space heating demand [kWh]	{:.3f}'.format(self.my_building.space_heating_energy / 1e3))
					print('      ' + 'Daily DHW demand [m3/day]		{:.3f}'.format(self.my_building.daily_DHW))
					print('      ' + 'DHW demand per liv_area [m3/m2]   {:.3f}'.format(sum(self.my_building.hot_water_m3) / self.my_building.living_area))
					print('      ' + 'Hot water demand [m3]		{:.3f}'.format(sum(self.my_building.hot_water_m3)))
					print('      ' + 'Hot water demand [kWh]		{:.3f}'.format(self.my_building.hot_water_energy / 1e3))
					print('      ' + 'Agg. heat demand [kWh]		{:.3f}'.format(self.my_building.total_energy / 1e3))
					print('      ' + 'Distance to heat plant [m]	{:.3f}'.format(self.my_building.distance2hp))
					print('      ' + 'Delay due to distance [min]	{:.3f}'.format(self.my_building.delay))
					print('      ' + 'Agg. delayed heat demand [kWh]	{:.3f}'.format(self.my_building.total_energy_delayed / 1e3))
					print('      ' + '==========================================================================================================\n')
				
				# plot results per building
				if (self.plot == 2):
					self.my_building.plot_timeseries(space_heating = self._space_heating, Tb = True, \
											hot_water = self._hot_water, total = True)
				
				# save results per building
				if (self.save >= 2):
					self.my_building.save_csv()
					#self.my_building.save_load_duration_curve()
				
				if (self.save == 3):
					if hot_water:
						self.my_building.save_dhw_debug_csv()
				
				# add building demand to city demand
				self.total_power  		 += self.my_building.total_power
				self.total_power_delayed += self.my_building.total_power_delayed
							
		# add base load
		self.total_power 		 += self.base_load
		self.total_energy 		  = self.total_power.sum() * self.resolution * 1 / 60
		self.total_power_delayed += self.base_load
		self.total_energy_delayed = self.total_power_delayed.sum() * self.resolution * 1 / 60
		
		# plot results per city
		if (self.plot >= 1):
			self.plot_timeseries(space_heating = self._space_heating, hot_water = self._hot_water, total = True)
			
		# save results per city
		if (self.save >= 1):
			self.save_csv_power()
			self.save_csv_energy()

			
# Domestic hot water demand
# --------------------------------------------------------------------------------
	def initialize_dhw_probabilities(self):
		"""
		Calculates dhw probabilities (daily consumption, event loads, flow rate and 
		duration as interpolate objects.
		"""
		
		# Daily specific dhw consumption [m3/m2 of living area]	as cdf	
		x 		= [i[0] for i in self.building_stock_stats[1][0]]
		p 		= [i[1] for i in self.building_stock_stats[1][0]]
		dhw_cdf = UrbanHeatPro.create_interpolated_cdf(x, p)
		
		# Probability distribution of the DHW-load happening at specific time of the day
		#
		### Probabilities in 0.2h-steps
		p_shower 	 = [i[1] for i in self.building_stock_stats[1][1]]
		p_bath 		 = [i[2] for i in self.building_stock_stats[1][1]]
		p_med_small  = [i[3] for i in self.building_stock_stats[1][1]]
		#
		### Interpolators
		x_min    	 = 0.
		x_max    	 = 23.8
		x_steps  	 = 120
		x	     	 = np.linspace(x_min, x_max, num = x_steps, endpoint = False)
		shower_ptime = interpolate.interp1d(x, p_shower)
		bath_ptime 	 = interpolate.interp1d(x, p_bath)
		medium_ptime = interpolate.interp1d(x, p_med_small)
		small_ptime  = medium_ptime
		#
		prob_time 	 = [shower_ptime, bath_ptime, medium_ptime, small_ptime]
		
		# Factor for probability distribution of the DHW-load happening at specific weekday
		shower_pwday = [i[1] for i in self.building_stock_stats[1][2]]
		bath_pwday   = [i[2] for i in self.building_stock_stats[1][2]]
		medium_pwday = [i[3] for i in self.building_stock_stats[1][2]]
		small_pwday  = medium_pwday
		prob_wday	 = [shower_pwday, bath_pwday, medium_pwday, small_pwday]
		
		# Load flow rate CDF
		x_min    	 = 3.
		x_max    	 = 20.
		x_steps  	 = 1000
		x	     	 = np.linspace(x_min, x_max, num = x_steps, endpoint = False)
		### Shower
		mean  		 = self.building_stock_stats[1][3][0][0]
		sigma 		 = self.building_stock_stats[1][3][1][0]
		norm_dist	 = UrbanHeatPro.create_normal_distribution(x, mean, sigma)
		shower_f_cdf = UrbanHeatPro.create_interpolated_cdf(x, norm_dist)
		#
		### Bath
		mean 		 = self.building_stock_stats[1][3][0][1]
		sigma		 = self.building_stock_stats[1][3][1][1]
		norm_dist	 = UrbanHeatPro.create_normal_distribution(x, mean, sigma)
		bath_f_cdf   = UrbanHeatPro.create_interpolated_cdf(x, norm_dist)
		#
		x_min    	 = 0.
		x_max    	 = 10.
		x_steps  	 = 1000
		x	     	 = np.linspace(x_min, x_max, num = x_steps, endpoint = False)
		### Medium
		mean 		 = self.building_stock_stats[1][3][0][2]
		sigma		 = self.building_stock_stats[1][3][1][2]
		norm_dist	 = UrbanHeatPro.create_normal_distribution(x, mean, sigma)
		medium_f_cdf = UrbanHeatPro.create_interpolated_cdf(x, norm_dist)
		#
		### Small
		mean 		 = self.building_stock_stats[1][3][0][3]
		sigma		 = self.building_stock_stats[1][3][1][3]
		norm_dist	 = UrbanHeatPro.create_normal_distribution(x, mean, sigma)
		small_f_cdf  = UrbanHeatPro.create_interpolated_cdf(x, norm_dist)
		#
		flowRate_cdf = [shower_f_cdf, bath_f_cdf, medium_f_cdf, small_f_cdf]
		
		# Load duration CDF
		x_min    	 = 3.
		x_max    	 = 15.
		x_steps  	 = 1000
		x	     	 = np.linspace(x_min, x_max, num = x_steps, endpoint = False)
		### Shower
		mean  		 = self.building_stock_stats[1][4][0][0]
		sigma 		 = self.building_stock_stats[1][4][1][0]
		norm_dist	 = UrbanHeatPro.create_normal_distribution(x, mean, sigma)
		shower_d_cdf = UrbanHeatPro.create_interpolated_cdf(x, norm_dist)
		#
		### Bath
		mean 		 = self.building_stock_stats[1][4][0][1]
		sigma		 = self.building_stock_stats[1][4][1][1]
		norm_dist	 = UrbanHeatPro.create_normal_distribution(x, mean, sigma)
		bath_d_cdf   = UrbanHeatPro.create_interpolated_cdf(x, norm_dist)
		#
		x_min    	 = 0.
		x_max    	 = 5.
		x_steps  	 = 1000
		x	     	 = np.linspace(x_min, x_max, num = x_steps, endpoint = False)
		### Medium
		mean 		 = self.building_stock_stats[1][4][0][2]
		sigma		 = self.building_stock_stats[1][4][1][2]
		norm_dist	 = UrbanHeatPro.create_normal_distribution(x, mean, sigma)
		medium_d_cdf = UrbanHeatPro.create_interpolated_cdf(x, norm_dist)
		#
		### Small
		mean 		 = self.building_stock_stats[1][4][0][3]
		sigma		 = self.building_stock_stats[1][4][1][3]
		norm_dist	 = UrbanHeatPro.create_normal_distribution(x, mean, sigma)
		small_d_cdf  = UrbanHeatPro.create_interpolated_cdf(x, norm_dist)
		#
		duration_cdf = [shower_d_cdf, bath_d_cdf, medium_d_cdf, small_d_cdf]
		
		dhw_prob = [dhw_cdf, prob_time, prob_wday, flowRate_cdf, duration_cdf]

		return dhw_prob
#
	def calculate_seasonal_variation_vector(self, amplitude = 0.1, max_point = 45):
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
		
		# calculate number of days in simulation
		start 	   = self.dt_vector[0].timetuple().tm_yday
		end   	   = self.dt_vector[-1].timetuple().tm_yday
		num_days   = (end - start) + 1
		
		# initialize day_vector
		day_vector = np.zeros([num_days, 3], dtype = int)
		
		# day_vector: [day in year, start time step]
		prev_day   = 0
		day_count  = 0
		for iii, date in enumerate(self.dt_vector):
			if not(date.timetuple().tm_yday == prev_day):
				day_vector[day_count][0] = date.timetuple().tm_yday
				day_vector[day_count][1] = iii	# index of dt_vector where a new day starts
				day_count += 1			
			prev_day   = date.timetuple().tm_yday
		
		# add end time step
		for iii, day in enumerate(day_vector):
			try:
				day_vector[iii][2] = day_vector[iii + 1][1] - 1
			except: # last day
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
		min_vector = np.zeros([self.nts], dtype = int)
		
		# calculate minutes of year
		for iii, date in enumerate(self.dt_vector):
			day_in_year = date.timetuple().tm_yday
			hour_in_day = date.timetuple().tm_hour
			min_in_hour = date.timetuple().tm_min
			
			min_vector[iii] = (day_in_year - 1) * 24 * 60 + hour_in_day * 60 + min_in_hour
			
		return min_vector

		
# Results		
# --------------------------------------------------------------------------------
	def plot_timeseries(self, space_heating = True, hot_water = True, total = True):
		"""
		"""
		
		if space_heating:
			fig_name = '{}/SpaceHeatingDemand_{}.png'.format(self.result_dir, self.rid)
			
			UrbanHeatPro.plot_timeseries(self.dt_vector, \
					[self.space_heating_power], ['Space heating demand'], \
					fig_name, xticks = ('month', 3), \
					ynumticks = 'auto', ylabel = 'Power [kW]', ylim0 = True, yfactor = 1e3)
													
		if hot_water:
			fig_name = '{}/HotWaterDemand_{}.png'.format(self.result_dir, self.rid)
			
			UrbanHeatPro.plot_timeseries(self.dt_vector, \
					[self.hot_water_power], ['Hot water demand'], \
					fig_name, xticks = ('month', 3), \
					ynumticks = 'auto', ylabel = 'Power [kW]', ylim0 = True, yfactor = 1e3)

		if total:
			fig_name = '{}/TotalHeatDemand_{}.png'.format(self.result_dir, self.rid)
			
			UrbanHeatPro.plot_stacked_timeseries(self.dt_vector, \
					[self.hot_water_power, self.space_heating_power], \
					['Hot water', 'Space heating'], \
					fig_name, xticks = ('month', 3), \
					ynumticks = 'auto', ylabel = 'Power [kW]', ylim0 = True, yfactor = 1e3)	
#
	def save_csv_syn_city(self, chunk_id, chunk_data):
		"""
		Saves key building parameters of every chunk.
		"""
			
		filename = '{}/SyntheticCity_{}.csv'.format(self.result_dir, self.rid)
		
		# Header
		if chunk_id == 0:

			with open(filename, 'w') as text_file:
				# Region
				text_file.write(self.region)
				# Building data
				text_file.write('bid;footprint_area;use;free_walls;' + \
								'lat;lon;dist_to_heat_source;' + \
								'NUTS_0;NUTS_1;NUTS_2;NUTS_3;' + \
								'year_class;size_class;' + \
								'ref_level_roof;ref_level_wall;ref_level_floor;ref_level_window;' + \
								'floors;area_corr_factor;storey_area;heated_area;' + \
								'dwellings;dwelling_size;occupants')
		
		# Building data
		with open(filename, 'a') as f:
			np.savetxt(f, chunk_data, delimiter = ';', 
						fmt = ['%d', '%.2f' ,'%d', '%d', \
							   '%.3f', '%.3f', '%.2f', \
							   '%s', '%s', '%s', '%s', \
							   '%d', '%d', \
							   '%d', '%d', '%d', '%d', \
							   '%.2f', '%.2f', '%.2f', '%.2f', \
							   '%d', '%.2f', '%d'])
#
	def save_csv_energy(self):
		"""
		Saves key building parameters and heat energy demand (space heating, hot water and 
		total).
		"""
			
		filename = '{}/EnergyPerBuilding_{}.csv'.format(self.result_dir, self.rid)
		
		# Header
		with open(filename, 'w') as text_file:
			text_file.write('bid;use;year_class;btype;footprint_area_m2;heatedArea_m2;livingArea_m2;' + \
							'floors;freeWalls;refurbishment;Tset_degC;' + 
							'Tb0_degC;actStart1_h;actEnd1_h;actStart2_h;actEnd2_h;occupied;' + 
							'U_W/K;V_W/K;C_J/K;Tau_h;SpaceHeatingDemand_Wh;' +
							'SpaceHeatingDemandPerHArea_kWh/m2;' +
							'SolarGainsPerHArea_kWh/m2;InternalGainsPerHArea_kWh/m2' +
							'DailyDHWLimit_m3/day;HotWaterDemand_m3;' +
							'HotWaterDemand_Wh;TotalHeatDemand_Wh;'+
							'distance2hp_m;delay_min;TotalHeatDemand_Delayed_Wh\n')
				
		# Building data
		with open(filename, 'a') as f:
			np.savetxt(f, self.results['buildings'], delimiter = ';', fmt='%.4f')
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
				 self.total_power_delayed]).transpose()
				
		filename = '{}/HeatDemandProfile_{}.csv'.format(self.result_dir, self.rid)
		
		# Header
		with open(filename, 'w') as text_file:
			text_file.write('Run id;{};\n'.format(self.rid))
			text_file.write('Total buildings;{};\n'.format(len(self.buildings)))
			text_file.write('Commercial;{};\n'.format(np.sum(self.energy_per_building[:,1] == 0)))
			text_file.write('Industrial;{};\n'.format(np.sum(self.energy_per_building[:,1] == 1)))
			text_file.write('Public;{};\n'.format(np.sum(self.energy_per_building[:,1] == 2)))
			text_file.write('Residential;{};\n'.format(np.sum(self.energy_per_building[:,1] == 3)))
			text_file.write('SpaceHeatingDemand_GWh;{};\n'.format(self.space_heating_energy / 1e9))
			text_file.write('HotWaterDemand_GWh;{};\n'.format(self.hot_water_energy / 1e9))
			text_file.write('TotalHeatDemand_GWh;{};\n'.format(self.total_energy / 1e9))
			text_file.write('datenum;Tamb_degC;SpaceHeatingDemand_MW;HotWaterDemand_MW;' +
							'TotalHeatDemand_MW;TotalHeatDemand_Delayed_MW\n')
		
		# Building data
		with open(filename, 'a') as f:
			np.savetxt(f, array_to_save, delimiter = ';', fmt='%.3f')