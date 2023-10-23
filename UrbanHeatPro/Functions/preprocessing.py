"""
preprocessing.py
A. Molar-Cruz @ TUM ENS
"""

import geopandas as gpd
import numpy as np


# --------------------------------------------------------------------------------
def shp_to_csv(my_dir, filename_buildings, debug):
    """
    Reads shapefile and saves a csv file with dbf data.
    """

    filename, buildings = read_building_data_from_shp(my_dir, filename_buildings, debug)

    with open(filename[:-4] + '.csv', 'a') as f:
        np.savetxt(f, buildings, fmt='%.3f, %s', delimiter='...', header='area; use', comments='Unterhaching')


def read_building_data_from_shp(my_dir, filename_buildings, debug):
    """
    Reads shapefile and calculates missing fields (bid, area, free_walls)
    """

    filename = '{}\\input\\Shapefiles\\{}'.format(my_dir, filename_buildings)

    if debug != 0:
        print('\nBuilding data (shp)')
        print('  ' + filename)

    # read shapefile
    buildings = gpd.read_file(filename)

    # project shapefile to epsg: 32632 (Germany) so that area is in m2
    buildings = buildings.to_crs({'init': 'epsg:32632'})

    # calculates missing fields
    flag_modified = False

    # calculate building id
    try:
        buildings['bid']
    except:
        calculate_bid(buildings)
        flag_modified = True
        if debug != 0:
            print('      ' + "'bid' field added")

    # calculate building area
    try:
        buildings['area']
    except:
        calculate_areas(buildings)
        flag_modified = True
        if debug != 0:
            print('      ' + "'area' field added")

    # calculate free walls
    try:
        buildings['free_walls']
    except:
        calculate_free_walls(buildings)
        flag_modified = True
        if debug != 0:
            print('      ' + "'free_walls' field added")

    # calculate latitude and longitude from building centroid
    try:
        buildings['lat']
    except:
        calculate_lat_lon(buildings)
        flag_modified = True
        if debug != 0:
            print('      ' + "'lat' and 'lon' fields added")

    # calculate distance to heat plant
    try:
        buildings['dist2hp']
    except:
        calculate_distance2hp(buildings, my_dir, filename_heat_plant)
        flag_modified = True
        if debug != 0:
            print('      ' + "'dist2hp' field added")

    # save shapefile with new fields
    if flag_modified:
        filename_new = filename[:-4] + '_mod.shp'
        buildings.to_file(filename_new)
        if debug != 0:
            print('  ' + "Shapefile with modifications saved as")
            print('    ' + filename_new)

    return filename, buildings.drop(columns='geometry').values


#
def calculate_bid(buildings):
    """
    Adds field "bid" with a consecutive building number as id
    """

    # add field "bid"
    buildings['bid'] = 0

    # add bid
    for iii in range(len(buildings)):
        buildings.ix[iii, 'bid'] = iii


#
def calculate_areas(buildings):
    """
    Adds field "area" with the polygon area in m2
    """

    buildings['area'] = buildings.geometry.area


#
def calculate_free_walls(buildings):
    """
    Adds field "free_walls" with the number of walls in direct contact with
    ambient temperature. It is assumed that all buildings have only four walls.
    """

    # add field "free walls"
    buildings['free_walls'] = 4

    # calculate free walls
    for iii in range(len(buildings)):
        intersected = buildings[buildings.geometry.intersects(buildings.ix[iii].geometry.buffer(0.1))]
        if len(intersected) > 1:  # the analyzed building is taken as intersected
            free_walls = 4 - (len(intersected) - 1)
            buildings.loc[iii, 'free_walls'] = free_walls


#
def calculate_lat_lon(buildings):
    """
    Adds field "lat" and "lot" with the latitude and longitude values of the
    building centroid in degrees.
    """

    # add fields "lat" and "lon"
    buildings['lat'] = 1.
    buildings['lon'] = 1.

    # project buildings to have lat/lon in degrees
    proj_buildings = buildings.to_crs(epsg=4326)

    # calculate centroids
    centroids = proj_buildings.centroid

    # calculate lat/lon
    for iii in range(len(centroids)):
        lat = centroids[iii].coords.xy[1][0]
        lon = centroids[iii].coords.xy[0][0]

        # add to building layer
        buildings.loc[iii, 'lat'] = lat
        buildings.loc[iii, 'lon'] = lon


#
def calculate_distance2hp(buildings, my_dir, filename_heat_plant):
    """
    Adds field "dist2hp" with the distance between building and heat plant.
    The heat plant location is given in a shapefile.
    """

    # add field "dist2hp"
    buildings['dist2hp'] = 0

    # import shapefile with heat plant
    filename = '{}/input/Shapefiles/{}'.format(my_dir, filename_heat_plant)
    hp = gpd.read_file(filename)
    hp = hp.to_crs({'init': 'epsg:32632'})
    hp = hp.ix[0].geometry

    # calculate distance to heat plant
    for iii in range(len(buildings)):
        centroid = buildings.ix[iii].geometry.centroid
        distance2hp = centroid.distance(hp)
        buildings.loc[iii, 'dist2hp'] = distance2hp
