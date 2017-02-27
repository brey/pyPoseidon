import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.mlab import griddata
from scipy.interpolate import griddata

from geopandas import GeoDataFrame
from shapely.geometry import Polygon, MultiPolygon

def collec_to_gdf(collec_poly):
    """Transform a `matplotlib.contour.QuadContourSet` to a GeoDataFrame"""
    polygons = []
    for i, polygon in enumerate(collec_poly.collections):
        mpoly = []
        for path in polygon.get_paths():
            try:
                path.should_simplify = False
                poly = path.to_polygons()
                # Each polygon should contain an exterior ring + maybe hole(s):
                exterior, holes = [], []
                if len(poly) > 0 and len(poly[0]) > 3:
                    # The first of the list is the exterior ring :
                    exterior = poly[0]
                    # Other(s) are hole(s):
                    if len(poly) > 1:
                        holes = [h for h in poly[1:] if len(h) > 3]
                mpoly.append(Polygon(exterior, holes))
            except:
                print('Warning: Geometry error when making polygon #{}'
                      .format(i))
        if len(mpoly) > 1:
            mpoly = MultiPolygon(mpoly)
            polygons.append(mpoly)
        elif len(mpoly) == 1:
            polygons.append(mpoly[0])
    return GeoDataFrame(geometry=polygons, crs={'init': 'epsg:4326'})

def convert_to_hex(rgba_color) :
    red = str(hex(int(rgba_color[0]*255)))[2:].capitalize()
    green = str(hex(int(rgba_color[1]*255)))[2:].capitalize()
    blue = str(hex(int(rgba_color[2]*255)))[2:].capitalize()

    if blue=='0':
        blue = '00'
    if red=='0':
        red = '00'
    if green=='0':
        green='00'

    return '#'+ red + green + blue
