import numpy as np
from netCDF4 import Dataset
import glob


d = Dataset('../gebco08.nc')
ntilesx=d.variables['ntilesx'][:]
ntilesy=d.variables['ntilesy'][:]

#global array
ni=43200
nj=21600

lats=np.empty((ntilesx[0],ntilesy[0]),dtype=object)
lons=np.empty((ntilesx[0],ntilesy[0]),dtype=object)
depths=np.empty((ntilesx[0],ntilesy[0]),dtype=object)

for i in range(1,ntilesx[0]+1):

   files=glob.glob('*.{:05}.*.nc'.format(i))
   files.sort()

   j=0
   for f in files:
     d = Dataset(f)
     lon=d.variables['lon'][:]
     lat=d.variables['lat'][:]
     depth=d.variables['depth'][:]

     lats[i-1,j]=lat 
     lons[i-1,j]=lon 
     depths[i-1,j]=depth
     j+=1



y=np.concatenate(lats[0,:])

x=np.concatenate(lons[:,0])

nrows = nj
ncols = ni

#grp = Dataset('out.nc', 'w', format='NETCDF3_64BIT')
grp = Dataset('out.nc', 'w', format='NETCDF4')
 
grp.createDimension('lat', nrows)
grp.createDimension('lon', ncols)
 
longitude = grp.createVariable('lon', 'f4', ('lon',))
latitude = grp.createVariable('lat', 'f4', ('lat',))
 
chunk_size=(nj,300)
elevation  = grp.createVariable('elevation', 'int16', ('lat', 'lon'), chunksizes=chunk_size)
 
longitude[:] = x
latitude[:] = y
grp.sync()
for i in range(144):
   elevation[:,300*i:300*(i+1)] = np.concatenate(depths[i,:])
   grp.sync()
 
grp.close()

