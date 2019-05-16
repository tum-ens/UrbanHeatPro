"""
to_tuple.py
AMC @ TUM ENS
"""

# --------------------------------------------------------------------------------
def building_use_to_tuple(use):
	"""
	"""
	if (type(use) == str):
		USE = {'commercial': 0, 'industrial': 1, 'public': 2, 'residential': 3}
		return (USE[use], str(use))
	elif type(use) == int:
		USE = {0: 'commercial', 1: 'industrial', 2: 'public', 3: 'residential'}
		return (use, USE[use])
# --------------------------------------------------------------------------------
def year_class_to_tuple(use_int, year_class_int):
	"""
	"""
	if (use_int == 3): # residential
		YEAR_CLASS = ['<1859',
					  '1860-1918',
					  '1919-1948',
					  '1949-1957',
					  '1958-1968',
					  '1969-1978',
					  '1979-1983',
					  '1984-1994',
					  '1995-2001',
					  '2002-2009',
					  '>2009']
	else: # non-residential
		YEAR_CLASS = ['<1918', 
			 		  '1919-1976',
					  '1977-1983',
					  '1984-1994',
					  '>1995']

	return (year_class_int, YEAR_CLASS[year_class_int])
# --------------------------------------------------------------------------------
def size_class_to_tuple(use_int, size_class_int):
	"""
	"""
	if (use_int == 3): # residential
		SIZE_CLASS 		= ['SFH', 'TH', 'MFH', 'AB']
		size_class_str  = SIZE_CLASS[size_class_int]
	
	else: # non-residential
		size_class_str = 'None'
		
	return (size_class_int, size_class_str)


	