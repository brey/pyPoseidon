import numpy as np

from dep import *
from grid import *

from netCDF4 import Dataset

PATH='../TIDES/'


grd=Grid.fromfile('../../../tmp2/med.grd')
dlat=grd.y[1,0]-grd.y[0,0]
dlon=grd.x[0,1]-grd.x[0,0]


ba=Dep.read('../../../tmp2/med.dep',grd.shape)

west=[]
for j in range(grd.y.shape[0]):
  if ba.val[j,0] > 0. : west.append(j)

with open('../../../tmp2/med.bnd', 'w') as f:
   f.write('West1                Z A     1    {}     1    {}   0.0000000e+00 West1A West1B'.format(west[0]+1,west[1]+1)) # fortran index ??


ll=['A','B']


#read from med.nc
dmed=Dataset(PATH+'med.nc')

lat=dmed['lat'][:]
lon=dmed['lon'][:]

tidal_c=dmed['tidal_constituents'][:]

tidal_c=[''.join(k).upper().strip() for k in tidal_c]

with open('../../../tmp2/med.bca', 'w') as f:

 for k in range(2):

  plon=grd.x[west[k],0]
  plat=grd.y[west[k],0]
  i=np.abs(lon-np.float(plon)).argmin()
  j=np.abs(lat-np.float(plat)).argmin()

  amp=dmed['tidal_amplitude_h'][i,j,:]
  ph=dmed['tidal_phase_h'][i,j,:]

  f.write('West1{}\n'.format(ll[k]))
  for a,b,c in zip(tidal_c,amp,ph):
     f.write('{}         {:.7e}   {:.7e}\n'.format(a,b,c))


#

ini=np.zeros(grd.x.shape)     
np.savetxt('../../../tmp2/med.ini',ini)
with open('../../../tmp2/med.ini', 'a') as f:
    np.savetxt(f,ini)
    np.savetxt(f,ini)

