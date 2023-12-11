"""
to_tuple.py
A. Molar-Cruz @ TUM ENS
"""


# --------------------------------------------------------------------------------
def building_use_to_tuple(use):
    """
    Map building use to tuple (use_int, use_str)
    
    :param use: building use as string or integer
    :returns: tuple (use_int, use_str)
    """
    if type(use) == str:
        USE = {'commercial': 0, 'industrial': 1, 'public': 2, 'residential': 3}
        return USE[use], str(use)
    elif type(use) == int:
        USE = {0: 'commercial', 1: 'industrial', 2: 'public', 3: 'residential'}
        return use, USE[use]


# --------------------------------------------------------------------------------
def year_class_to_tuple(use_int, year_class_int):
    """
    Map year class to tuple (year_class_int, year_class_str)
    
    :param use_int: building use as integer
    :param year_class_int: year class as integer
    :returns: tuple (year_class_int, year_class_str)
    """
    if use_int == 3:  # residential
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
    else:  # non-residential
        YEAR_CLASS = ['<1918',
                      '1919-1976',
                      '1977-1983',
                      '1984-1994',
                      '>1995']

    return year_class_int, YEAR_CLASS[year_class_int]


# --------------------------------------------------------------------------------
def size_class_to_tuple(use_int, size_class_int):
    """
    Map size class to tuple (size_class_int, size_class_str)
    
    :param use_int: building use as integer
    :param size_class_int: size class as integer
    :returns: tuple (size_class_int, size_class_str)
    """
    if use_int == 3:  # residential
        SIZE_CLASS = ['SFH', 'TH', 'MFH', 'AB']
        size_class_str = SIZE_CLASS[size_class_int]

    else:  # non-residential
        size_class_str = 'None'

    return size_class_int, size_class_str
