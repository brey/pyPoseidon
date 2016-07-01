import numpy as np
from netCDF4 import Dataset
from grid import *


rootgrp = Dataset('test.nc', 'w', format='NETCDF3_64BIT')
lat = rootgrp.createDimension('lat', 260)
lon = rootgrp.createDimension('lao', 700)
time = rootgrp.createDimension('time', None)
level = rootgrp.createDimension('level', None)


times = rootgrp.createVariable('time','f8',('time',))
levels = rootgrp.createVariable('level','i4',('level',))
latitudes = rootgrp.createVariable('latitude','f4',('lat',))
longitudes = rootgrp.createVariable('longitude','f4',('lon',))

rootgrp.description = ''
rootgrp.history = 'DELFT3D - JRC Ispra European Commission'
rootgrp.source = 'netCDF4 python module tutorial'
latitudes.units = 'degrees north'
latitudes.point_spacing = 'even'
longitudes.units = 'degrees east'
longitudes.point_spacing = 'even'
levels.units = 'm'
times.units = 'seconds'


latitudes[:]=lat
longitudes[:]=lon
levels[:]=ha
times[:]=tsn


rootgrp.close()
