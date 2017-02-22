import numpy as np
import sys
from netCDF4 import Dataset

def write2nc(filename,lat,lon,zy,zx,tstamp,t0,val):

 ni=lon.shape[0]
 nj=lat.shape[0]

 rootgrp = Dataset(filename, 'w', format='NETCDF3_64BIT')
 lats = rootgrp.createDimension('LAT', nj)
 lons = rootgrp.createDimension('LON', ni)
 time = rootgrp.createDimension('TIME', None)


 longitudes = rootgrp.createVariable('LON','f8',('LON',))
 latitudes = rootgrp.createVariable('LAT','f8',('LAT',))
 lon_masked = rootgrp.createVariable('ZX','f8',('LON','LAT'))
 lat_masked = rootgrp.createVariable('ZY','f8',('LON','LAT'))
 times = rootgrp.createVariable('TIME','f8',('TIME',))
 levels = rootgrp.createVariable('HA','f8',('TIME','LON','LAT'))

 rootgrp.description = ''
 rootgrp.history = 'DELFT3D - JRC Ispra European Commission'
 rootgrp.source = 'netCDF4 python module tutorial'
 latitudes.units = 'degrees_north'
 latitudes.point_spacing = 'even'
 longitudes.units = 'degrees_east'
 longitudes.point_spacing = 'even'
 lat_masked.units = 'degrees_north'
 lat_masked.point_spacing = 'even'
 lon_masked.units = 'degrees_east'
 lon_masked.point_spacing = 'even'
 levels.units = 'm'
 times.units = 'seconds since {}'.format(t0)


 levels[:]=val
 times[:]=tstamp
 latitudes[:]=lat
 longitudes[:]=lon
 lat_masked[:]=zy
 lon_masked[:]=zx

 rootgrp.close()


if __name__ == "__main__":
    t=sys.argv[1]
    t0=sys.argv[2]
    ha=sys.argv[3]
    lat=sys.argv[4]
    lon=sys.argv[5]
    zx=sys.argv[6]
    zy=sys.argv[7]
    fname=sys.argv[8]
    
    write2nc(fname,lat,lon,zy,zx,t,t0,ha)
