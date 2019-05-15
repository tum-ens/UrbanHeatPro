"""
probabilistic.py
AMC @ TUM ENS
"""

import numpy  as np
from datetime import datetime
from scipy 	  import interpolate

# --------------------------------------------------------------------------------
def create_normal_distribution(x, mean, sigma):
	"""
	Returns a vector of len(x) values distributed normally around the specific
	mean and standard deviation sigma.
	"""
	
	# initialize vector
	normdist = np.zeros(len(x))
	
	# calculate normally distributed values
	for iii, value in enumerate(x):
		normdist[iii] = (((2. * np.pi * sigma ** 2)) ** -1) * np.exp(-((value - mean) ** 2) / (2. * sigma ** 2))
	
	return normdist
#	
def create_interpolated_cdf(x, p):
	"""
	Returns the cummulative density function (cdf) as a vector of len(x) values
	given a probability function, p. The cdf interpolated in case the probability 
	function has less values than the given x vector.
	"""
	
	# calculate cdf
	cdf 	= np.cumsum(p)
	
	# normalize
	cdf_max = max(cdf)
	for iii, val in enumerate(cdf):
		cdf[iii] = val / cdf_max

	# create interpolators so that there is a cdf value for every x-value
	cdf	 	= interpolate.interp1d(cdf, x, 
				fill_value = (x[0], x[-1]), bounds_error = False) 
	
	return cdf