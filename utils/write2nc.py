import numpy as np
import sys
from netCDF4 import Dataset

def writenc(filename,lat,lon,t,t0,pval,uval,vval):#,sval,spval):

 nj,ni =lon.shape

 rootgrp = Dataset(filename, 'w', format='NETCDF3_64BIT')
 lats = rootgrp.createDimension('ny_grid', nj)
 lons = rootgrp.createDimension('nx_grid', ni)
 time = rootgrp.createDimension('time', None)


 times = rootgrp.createVariable('time','f8',('time',))
 longitudes = rootgrp.createVariable('lon','f8',('nx_grid','ny_grid'))
 latitudes = rootgrp.createVariable('lat','f8',('nx_grid','ny_grid',))
 p = rootgrp.createVariable('prmsl','f8',('time','nx_grid','ny_grid'))
 u = rootgrp.createVariable('uwind','f8',('time','nx_grid','ny_grid'))
 v = rootgrp.createVariable('vwind','f8',('time','nx_grid','ny_grid'))
 stmp = rootgrp.createVariable('stmp','f8',('time','nx_grid','ny_grid'))
 spfh = rootgrp.createVariable('spfh','f8',('time','nx_grid','ny_grid'))

 rootgrp.description = ''
 rootgrp.history = 'JRC Ispra European Commission'
 rootgrp.source = 'netCDF4 python module'

 latitudes.units = 'degrees_north'
 latitudes.long_name = 'Latitude'
 latitudes.standard_name = 'latitude'

 longitudes.units = 'degrees_east'
 longitudes.long_name = 'Longitude'
 longitudes.standard_name = 'longitude'

 p.units = 'Pa'
 p.long_name = 'Pressure reduced to MSL'
 p.standard_name = 'air_pressure_at_sea_level'

 u.units = 'm/s'
 u.long_name = 'Surface Eastward Air Velocity'
 u.standard_name = 'eastward_wind'

 v.units = 'm/s'
 v.long_name = 'Surface Northward Air Velocity'
 v.standard_name = 'northward_wind'

 times.units = 'seconds since {}'.format(t0)
 times.long_name = 'Time'
 times.standard_name = 'time'

 spfh.units = '1'
 spfh.long_name = 'Surface Specific Humidity (2m AGL)'
 spfh.standard_name = 'specific_humidity'

 stmp.units = 'degrees'
 stmp.long_name = 'Surface Temperature'
 stmp.standard_name = 'surface temperature'

 # transpose for complying with lon,lat format

 pf=np.zeros((t.shape[0],ni,nj))
 uf=np.zeros((t.shape[0],ni,nj))
 vf=np.zeros((t.shape[0],ni,nj))
 for k in range(t.shape[0]):
   pf[k,:,:]=pval[k,:,:].T
   uf[k,:,:]=uval[k,:,:].T
   vf[k,:,:]=vval[k,:,:].T

 times[:]=t
 longitudes[:]=lon.T
 latitudes[:]=lat.T
 p[:]=pf
 u[:]=uf
 v[:]=vf
#stmp[:]=sval
#spfh[:]=spval

 rootgrp.close()

