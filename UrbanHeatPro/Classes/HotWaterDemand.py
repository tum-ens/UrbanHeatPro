'''
HotWaterDemand.py
A. Molar-Cruz @ TUM ENS
'''

import os
import sys
import numpy as np
from datetime import datetime

class HotWaterDemand():
# --------------------------------------------------------------------------------
	def __init__(self, dt_vector, resolution, \
					day_vector, seasonal_vector, activity_vector, \
					Tw, daily_DHW, dhw_prob, \
					hw_tank_capacity, hw_tank_limit, hw_tank_volume_t0, hw_flow, \
					result_dir, use, year_class, btype, bid,
					debug, save_debug):
		
		# Mean daily consumption
		self.daily_DHW		   = daily_DHW				  # Mean daily hot water consumption [m3]
		
		# Hot water tank
		self.hw_tank_capacity  = hw_tank_capacity		  # Hot water tank capacity in m3
		self.hw_tank_limit	   = hw_tank_limit * self.hw_tank_capacity # Hot water tank limit in m3
		self.hw_tank_volume_t0 = hw_tank_volume_t0		  # Initial state of hot water tank in m3
		self.hw_tank_volume	   = self.hw_tank_volume_t0	  # State of hot water tank in m3
		self.hw_flow		   = hw_flow				  # Flow to refill hot water tank in L/min
		
		# Hot Water values
		self.Tw		      	   = Tw						  # Supply temperature of water in degC
		self.density	  	   = 1000					  # Water density in kg/m3
		self.c_p		  	   = 4186					  # Heat capacity of water in J/(kg degC)
		
		# DHW-loads probabilities
		self.flow_rate_cdf	   = dhw_prob[3]			  # CDF of dhw-load flow rate
		self.duration_cdf	   = dhw_prob[4]			  # CDF of dhw-load duration
		self.prob_daytime	   = dhw_prob[1]			  # Probability of dhw-load happening at specific time of day
		self.prob_weekday	   = dhw_prob[2]			  # Probability factor of dhw-load for dif days of the weekday
		
		# Simulation time frame
		self.dt_vector   	   = dt_vector	 			  # Vector of time steps as datetime objects
		self.nts         	   = len(self.dt_vector) 	  # Number of time steps
		self.resolution  	   = resolution				  # resolution in min
		self.day_vector		   = day_vector				  # Vector with day of year in simulaiton time frame
		self.seasonal_vector   = seasonal_vector		  # Sinusoidal function for seasonal variation of dhw consumption
		self.activity_vector   = activity_vector		  # Building activity vector (boolean)
			
		# Results
		self.dhw_m3	 	 	   = np.zeros([self.nts])  	  # Instantaneous DHW demand in m3
		self.dhw_tank_m3	   = np.zeros([self.nts])  	  # Hot water tank demand in m3
		self.dhw_power	 	   = np.zeros([self.nts])     # Domestic hot water demand in W
		self.dhw_debug	 	   = np.zeros([self.nts, 25]) # Domestic hot water debug data per time step
		
		#
		self.debug			   = debug
		self.save_debug		   = save_debug
		
		###
		self.agg_consumption = np.zeros([5])			  # Aggregated consumption

# Hot water demand
# --------------------------------------------------------------------------------		
	def calculate(self):
		'''
		Calculates time series of domestic hot water demand in m3/min and W for the 
		four dhw-load categories (0. Shower, 1. Bath, 2. Medium load and 3. Small load) 
		for the whole simulation time.
		'''
		
		if self.debug == 3:
			print('      ' + '= HOT WATER DEMAND')
			print('      ' + '==========================================================================================================')
				
		# calculate flow rate in m3
		for iii, day_info in enumerate(self.day_vector):
			day_in_year					     = day_info[0]
			start						     = day_info[1]
			end							     = day_info[2]
			
			# check if building is occupied
			occupied = sum(self.activity_vector[start:end + 1])
			
			# if building is occupied, create dhw demand events
			if occupied > 0:
			
				# initialize day consumption, consumption limit and flag
				self.daily_limit = self.daily_DHW * self.seasonal_vector[day_in_year - 1]
				self.daily_consumption_reached   = False
				self.daily_consumption 		     = 0.
				
				# add day data to debug matrix
				if self.save_debug:
					self.dhw_debug[start:end + 1, 0] = range(start, end + 1)
					self.dhw_debug[start:end + 1, 2] = iii
					self.dhw_debug[start:end + 1, 3] = day_in_year
					self.dhw_debug[start:end + 1, 4] = self.daily_DHW
					self.dhw_debug[start:end + 1, 5] = self.daily_limit
				
				#
				if self.debug == 3:
					print('			' + 'day: {}\tDay_in_year: {}\tMean_load: {:.3f}\tSeasonal_factor: {:.3f}\tDaily_limit [m3]: {:.3f}\HW Tank [m3]: {:.3f}'.format(iii, \
											day_in_year, self.daily_DHW, self.seasonal_vector[day_in_year - 1], self.daily_limit, self.hw_tank_volume))
				
				# Loop until daily consumption is reached
				loop_counter = 1
				while not(self.daily_consumption_reached):
					
					#
					if self.debug == 3:				
						print('          ' + '----------------------------------------')
						print('          ' + 'Loop\tLimit\tConsumption\tTankState\tReached?')
						print('          ' + '{}\t{:.3f}\t{:.3f}\t\t{:.3f}\t\t{}'.format(loop_counter, self.daily_limit, \
												self.daily_consumption, self.hw_tank_volume/self.hw_tank_capacity, self.daily_consumption_reached))
						print('            ' + '____________________________________________________________________________________________________')
						print('            ' + 'time_step\t\tActive?\tprob\trndnum\tload\tevent?\tflow\tdur\tdur_ts\tcons\t\tagg.cons')
					
					# Loop through timesteps in day
					for jjj, time_step in enumerate(self.dt_vector[start:end + 1]):
					
						# Time step info
						if self.save_debug:
							self.dhw_debug[start + jjj][1]  = time_step.hour + time_step.minute/60.0
						
						# calculate dhw-events for every time step
						self.calculate_dhw_consumption_in_time_step(time_step, start + jjj, iii, day_in_year)		
						
						if self.daily_consumption_reached:
							break
					loop_counter += 1
									
				# calculate convergence of daily flow [m3/day]
				if self.save_debug:
					self.dhw_debug[start:end + 1, 24] 	   = self.agg_consumption[4] / day_in_year
						
			# calculate power consumption in W
			# power [W] = flow_rate [m3/min] * duration [min] * density [kg/m3] * c_p [J/(kg degC)] * delta_temperature_water [degC]
			Tw_cold = 6 # degC
			Tw_hot  = self.Tw
			dTw = Tw_hot - Tw_cold
			for iii, date in enumerate(self.dt_vector):
				#self.dhw_power[iii] = self.dhw_m3[iii] * self.density * self.c_p * dTw / (60 * self.resolution)
				self.dhw_power[iii] = self.dhw_tank_m3[iii] * self.density * self.c_p * dTw / (60 * self.resolution)
			#
			if self.debug == 3:
				print('      ' + '==========================================================================================================\n')
		
#	
	def calculate_dhw_consumption_in_time_step(self, time_step, dt_vector_index, day, day_in_year):
		'''
		Calculates hot water consumption events in every time step for the load categories:
			0   shower
			1   bath
			2   medium load
			3   small load
			
		Args:
			time_step 			time step as datetime object
			dt_vector_index		index of time step in simulation time frame dt_vector
			day					number of day in day_vector
			day_in_year			number of day in year (1-366)
		'''
		
		# check if building is active
		if self.activity_vector[dt_vector_index]:
		
			# add occupancy to debug matrix
			if self.save_debug:
				self.dhw_debug[dt_vector_index][6] = 1
		
			# Loop for dhw-load categories
			load_order = [0, 3, 2, 1] # shower, small, medium, bath
			for load in load_order:
				# get time step information
				hour    = time_step.hour + time_step.minute/60.0  # Float value [0,24)
				weekday = time_step.weekday()   # Integer value [0: Monday, 6: Sunday]
				
				# calculate probability of dhw-load event happening
				# Value is scaled down with the weekday factor because the value stored in prob_daytime
				# corresponds to Sunday, the day with the highest probability of water usage.
				prob_event  = self.prob_daytime[load](hour) * self.prob_weekday[load][weekday] / self.prob_weekday[load][6]

				# calculate if dhw consumption event is taking place
				rand_num_e = np.random.uniform(0, 1, 1)[0]
				if rand_num_e <= prob_event: 

					# calculate event flow rate in m3/min
					rand_num    = np.random.uniform(0, 1, 1)
					flow_rate   = self.flow_rate_cdf[load](rand_num)[0] / 1000
					
					# calculate event duration in min
					rand_num    = np.random.uniform(0, 1, 1)
					duration    = self.duration_cdf[load](rand_num)[0]
					
					if duration > 0:

						# calculate event hot water consumption in m3
						consumption 		        = flow_rate * duration

						### scale-up consumption to reduce number of loops
						min_share 					= 0.1	# min share of daily limit
						consumption 				= max(consumption, self.daily_limit * min_share)
						
						# calculate event duration in time_steps
						duration_in_time_steps      = int(np.ceil(duration / self.resolution))
										
						# calculate volume per time_step
						volume_time_step 	        = consumption / duration_in_time_steps
						
						# add event consumption to time series
						self.dhw_m3[dt_vector_index:dt_vector_index + duration_in_time_steps] += volume_time_step					
												
						# add event consumption to daily hot water consumption
						self.daily_consumption     += consumption
						self.agg_consumption[load] += consumption
						self.agg_consumption[4]    += consumption
						
						# add event to debug matrix
						if self.save_debug:
							self.dhw_debug[dt_vector_index, 7  + load] = 1
							self.dhw_debug[dt_vector_index, 11 + load] = flow_rate
							self.dhw_debug[dt_vector_index, 15 + load] = duration
							self.dhw_debug[dt_vector_index:dt_vector_index + duration_in_time_steps, 19 + load]   = volume_time_step
						
						# subtract volume from hot water tank
						self.hw_tank_volume        -= consumption
						
						# if tank has reached limit, refill tank
						if self.hw_tank_volume     <= self.hw_tank_limit:
							
							# calculate volume to refill per time step
							volume_to_refill          = self.hw_tank_capacity - self.hw_tank_volume
							dur_refill_in_time_steps  = int(np.ceil(volume_to_refill / self.hw_flow * (self.resolution / 1000.)))
							vol_refill_time_step 	  = volume_to_refill / dur_refill_in_time_steps
							
							# add refill event to time series
							self.dhw_tank_m3[dt_vector_index:dt_vector_index + dur_refill_in_time_steps] += vol_refill_time_step
							
							# set hw_tank_volume to 100% capacity
							self.hw_tank_volume = self.hw_tank_capacity
							
						#
						if self.debug == 3:
							print('            ' + '{}\t{}\t{:.3f}\t{:.3f}\t{}\t{}\t{:.3f}\t{:.3f}\t{}\t{:.3f}\t\t{:.3f}'.format(time_step, 1, prob_event, rand_num_e, load, 1, flow_rate, duration, duration_in_time_steps, consumption, self.daily_consumption))
						
						# check daily limit (80% of nominal load)
						if self.daily_consumption >= self.daily_limit * 0.8:
							self.daily_consumption_reached = True
							
							#
							if self.debug == 3:
								print('          ' + '______________________________________________________________________________________________________')
								print('            ' + '*** Daily consumption reached: {:.3f} >= {:.3f}***\n'.format(self.daily_consumption, self.daily_limit))
							
				# No dhw-load event because rand_num_e > prob_event
				else:
					if self.debug == 3:
						print('            ' + '{}\t{}\t{:.3f}\t{:.3f}\t{}\t{}\t{}\t{}\t{}\t{}\t\t{}'.format(time_step, 1, prob_event, rand_num_e, load, 0, 0, 0, 0, 0, 0))
				
				if self.daily_consumption_reached:
					if self.save_debug:
						self.dhw_debug[dt_vector_index, 23] 	   = self.daily_consumption
					break
		
		else:
			#
			if self.debug == 3:
				print('            ' + '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t\t{}'.format(time_step, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
