"""
plot.py
AMC @ TUM ENS
"""

import numpy as np
import os

from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.dates import MonthLocator, DayLocator, HourLocator, DateFormatter


# --------------------------------------------------------------------------------
def plot_timeseries(dt_vector, timeseries, legend, fig_name, \
					xticks = ('month', 3), \
					ynumticks = 'auto', ylabel = 'Power [kW]', ylim0 = True, yfactor = 1e3):
	"""
	Plots timeseries as steps
	
	args:
		dt_vector	<list>		list of datetime objects
		timeseries	<list>		timeseries to plot [ts1, ts2, ts3]
		legend		<list>		list of legends [leg1, leg2, leg3]
		fig_name	<string>	figure name
		xticks		<tuple>		('month', X), ('day', X), ('hour', X)
								Every X months/days/hours
		ynumticks
		ylabel
		ylim0
		yfactor
		
	"""
	
	# define plot style
	style_dir = '{}/input/Styles/'.format(os.getcwd())
	TUM = style_dir + 'TUM.mplstyle'
	presentation = style_dir + 'presentation.mplstyle'
	plt.style.use([TUM, presentation])
	
	# define figure and axes
	fig, ax = plt.subplots()
	
	# plot
	for iii, ts in enumerate(timeseries):
		ax.plot(dt_vector, ts / yfactor, label = legend[iii])
	ax.legend(bbox_to_anchor = (0., 1.02, 1., .102), loc = 1, borderaxespad = 0., \
								ncol = len(legend))
	
	# x-axis
	if xticks[0]   == 'month':
		major_ticks = MonthLocator(range(1, 13), bymonthday = 1, interval = xticks[1])
		minor_ticks = MonthLocator(range(1, 13), bymonthday = 1, interval = 1)
		ticksFmt    = DateFormatter("%b '%y")
	elif xticks[0] == 'day':
		major_ticks	= DayLocator(range(1, 32), interval = xticks[1])
		minor_ticks	= DayLocator(range(1, 32), interval = 1)
		ticksFmt    = DateFormatter("%d.%m")
	elif xticks[0] == 'hour':
		major_ticks = HourLocator(range(24), interval = xticks[1])
		minor_ticks = HourLocator(range(24), interval = 1)
		ticksFmt    = DateFormatter("%H")
	
	ax.xaxis.set_major_locator(major_ticks)
	ax.xaxis.set_minor_locator(minor_ticks)
	ax.xaxis.set_major_formatter(ticksFmt)
	ax.autoscale_view()
	
	# y-axis
	if ylim0:
		start, end = ax.get_ylim()
		ax.set_ylim([0, end])
	if ynumticks != 'auto':
		start, end = ax.get_ylim()
		step = int(round((end - start) / ynumticks, -1))
		ax.yaxis.set_ticks(np.arange(start, step * ynumticks, step))
	ax.set_ylabel(ylabel)
	
	# save figure
	fig.savefig(fig_name, dpi = 300)
	
	# close figure
	plt.close(fig)
	
	# back to classic plot style
	plt.style.use('classic')
#	
def plot_stacked_timeseries(dt_vector, timeseries, legend, fig_name, \
					xticks = ('month', 3), \
					ynumticks = 'auto', ylabel = 'Power [kW]', ylim0 = True, yfactor = 1e3):
	"""
	Plots timeseries as steps
	
	args:
		dt_vector	<list>		list of datetime objects
		timeseries	<list>		timeseries to stack [ts1, ts2]
		legend		<list>		list of legends [leg1, leg2, leg3]
		fig_name	<string>	figure name
		xticks		<tuple>		('month', X), ('day', X), ('hour', X)
								Every X months/days/hours
		ynumticks
		ylabel
		ylim0
		yfactor
		
	"""
	# define plot style
	style_dir = '{}/input/Styles/'.format(os.getcwd())
	TUM = style_dir + 'TUM.mplstyle'
	presentation = style_dir + 'presentation.mplstyle'
	plt.style.use([TUM, presentation])
	
	# define figure and axes
	fig, ax = plt.subplots()
	
	# plot
	ax.stackplot(dt_vector, timeseries)
	ax.legend(legend, bbox_to_anchor = (0., 1.02, 1., .102), loc = 1, borderaxespad = 0., \
								ncol = len(legend))
	
	# x-axis
	if xticks[0]   == 'month':
		major_ticks = MonthLocator(range(1, 13), bymonthday = 1, interval = xticks[1])
		minor_ticks = MonthLocator(range(1, 13), bymonthday = 1, interval = 1)
		ticksFmt    = DateFormatter("%b '%y")
	elif xticks[0] == 'day':
		major_ticks	= DayLocator(range(1, 32), interval = xticks[1])
		minor_ticks	= DayLocator(range(1, 32), interval = 1)
		ticksFmt    = DateFormatter("%d.%m")
	elif xticks[0] == 'hour':
		major_ticks = HourLocator(range(24), interval = xticks[1])
		minor_ticks = HourLocator(range(24), interval = 1)
		ticksFmt    = DateFormatter("%H")
	
	ax.xaxis.set_major_locator(major_ticks)
	ax.xaxis.set_minor_locator(minor_ticks)
	ax.xaxis.set_major_formatter(ticksFmt)
	ax.autoscale_view()
	
	# y-axis
	if ylim0:
		start, end = ax.get_ylim()
		ax.set_ylim([0, end])
	if ynumticks != 'auto':
		start, end = ax.get_ylim()
		step = int(round((end - start) / ynumticks, -1))
		ax.yaxis.set_ticks(np.arange(start, step * ynumticks, step))
	ax.set_ylabel(ylabel)
	
	# save figure
	fig.savefig(fig_name, dpi = 300)
	
	# close figure
	plt.close(fig)
	
	# back to classic plot style
	plt.style.use('classic')
#
def plot_histogram(values, ylabel, fig_name, factor = 1e3, statistics = []):
	"""
	Plots simple histogram
	"""
	
	# define plot style
	style_dir = '{}/input/Styles/'.format(os.getcwd())
	TUM = style_dir + 'TUM.mplstyle'
	presentation = style_dir + 'presentation.mplstyle'
	plt.style.use([TUM, presentation])
	
	# define figure and axes
	fig, ax = plt.subplots(figsize = (10, 5))
	
	# plot hist
	ax.hist(values / factor, color = (0.8, 0.8, 0.8))
	ax.set_ylabel(ylabel)
				
	# mean
	mean = (values / factor).mean()
	ax.axvline(mean, linestyle = 'dashed', linewidth = 2)
					
	# save figure
	fig.savefig(fig_name, dpi = 300)
	
	# close figure
	plt.close(fig)
	
	# back to classic plot style
	plt.style.use('classic')
#
def plot_histogram_table(use, thermal_property, title, fig_name,
								factor = 1e3, statistics = [],
								figsize = (30, 25)):
	"""
	Plots a histogram showing the values of the thermal properties for all buildings in 
	the city (only residential). A histogram per year construction class and building type 
	is shown.
	"""
	
	thermal_property = np.array(thermal_property)
	
	if use == 3: # RESIDENTIAL

		# define figures and axes
		fig, ax 	= plt.subplots(10, 4, figsize = figsize, tight_layout = True)
		
		# Title
		fig.suptitle(title)
		
		# Subplots
		for iii in range(10):
			for jjj in range(4):
										
				# histogram
				if (thermal_property[iii][jjj].size != 0):
					ax[iii, jjj].hist(thermal_property[iii][jjj] / factor, color = (0.8, 0.8, 0.8))
				
					# mean
					mean = (thermal_property[iii][jjj] / factor).mean()
					ax[iii, jjj].axvline(mean, color = 'b', linestyle = 'dashed', linewidth = 2)
					
					# statistics
					if statistics != []:
						stat = (statistics[iii][jjj]).mean()
						ax[iii, jjj].axvline(stat, color = 'r', linestyle = 'dashed', linewidth = 2)
					
	else: # NON-RESIDENTIAL

		# define figures and axes
		fig, ax 	= plt.subplots(5, 1, tight_layout = True)
		
		# Title
		fig.suptitle(title)
		
		# Subplots
		for iii in range(5):
													
			# plot
			if (thermal_property[iii].size != 0):
				ax[iii].hist(thermal_property[iii] / factor, color = (0.8, 0.8, 0.8))
				
				# mean
				mean = (thermal_property[iii] / factor).mean()
				ax[iii].axvline(mean, color = 'b', linestyle = 'dashed', linewidth = 2)
				
				# statistics
				if statistics != []:
					stat = (statistics[iii]).mean()
					ax[iii].axvline(stat, color = 'r', linestyle = 'dashed', linewidth = 2)

	# save figure
	fig.savefig(fig_name, dpi = 300)
	
	# close figure
	plt.close(fig)
#
def plot_imshow_comparison(use, sim_stock, stat_stock, fig_name, cmap = 'RdBu'):
	"""
	Shows a figure with two tables:
		Left	The distribution of residential buildings in the diff categories (year_class, btype)
				according to the statistics used to generate the synthetic building stock
		Right 	The distribution of residential buildings in the diff categories (year_class, btype)
				in the synthetic building stock
	"""
	
	# Ticks
	
	if use == 3: # RESIDENTIAL
		xticks_value = [0, 1, 2, 3]
		xticks_label = ['SFH', 'TH', 'MFH', 'AB']
		yticks_value = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
		yticks_label = ['<1859', 
						'1860-1918',
						'1919-1948',
						'1949-1957',
						'1958-1968',
						'1969-1978',
						'1979-1983',
						'1984-1994',
						'1995-2001',
						'2002-2009']
		
		max = 0.1
	
	
	else:
		xticks_value = [0.5]
		xticks_label = ['All']
		yticks_value = [0, 1, 2, 3, 4]
		yticks_label = ['<1918',
						'1919-1976',
						'1977-1983',
						'1984-1994',
						'>1995']
		
		max = 0.25
		
		# create [5,2] array to plot
		# 1D arrays cannot be used with plt.imshow()
		stat_stock = np.array([stat_stock, stat_stock]).transpose()
		sim_stock  = np.array([sim_stock, sim_stock]).transpose()

	# define figure and axes
	fig, ax = plt.subplots(1, 2, subplot_kw = {'xticks': xticks_value, 'yticks': yticks_value}, sharey = True)
	
	# Subplots
	# Ax1 - Building stock from statistics
	
	im = ax[0].imshow(stat_stock, interpolation = 'none', cmap = cmap, vmin = 0, vmax = max, aspect = 'equal')
	ax[0].set_title('Statistics')
	ax[0].set_xticklabels(xticks_label)
	ax[0].set_yticklabels(yticks_label)
	
	# Ax2 - Building stock from simulation
	ax[1].imshow(sim_stock, interpolation = 'none', cmap = cmap, vmin = 0, vmax = max, aspect = 'equal')
	ax[1].set_title('Simulation')
	ax[1].set_xticklabels(xticks_label)
	ax[1].set_yticklabels(yticks_label)

	# Colorbar
	fig.colorbar(im, ax = ax.ravel().tolist())
		
	# save figure
	fig.savefig(fig_name, dpi = 300)
	
	# close figure
	plt.close(fig)
#

	