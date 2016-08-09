import numpy as np

from dep import *
from grid import *

from netCDF4 import Dataset

from bilinear import bilinear_interpolation

PATH='../TIDES/'


grd=Grid.fromfile('../../../tmp2/med.grd')
dlat=grd.y[1,0]-grd.y[0,0]
dlon=grd.x[0,1]-grd.x[0,0]


ba=Dep.read('../../../tmp2/med.dep',grd.shape)

west=[]
for j in range(grd.y.shape[0]):
  if ba.val[j,0] > 0. : west.append(j)

# take concecutive points
swest=[]
swest.append(west[0])
for k in xrange(1,np.size(west)):
  if west[k]-1 in swest : swest.append(west[k])


with open('../../../tmp2/med.bnd', 'w') as f:
   f.write('West1                Z A     1    {}     1    {}   0.0000000e+00 West1A West1B'.format(swest[0]+1,swest[-1]+1)) # fortran index ??


le=['A','B']


#read from med.nc
dmed=Dataset(PATH+'med.nc')

lat=dmed['lat'][:]
lon=dmed['lon'][:]

tidal_c=dmed['tidal_constituents'][:]

tidal_c=[''.join(k).upper().strip() for k in tidal_c]

amp=dmed['tidal_amplitude_h']
ph=dmed['tidal_phase_h']

with open('../../../tmp2/med.bca', 'w') as f:

 for k,l in zip([swest[0],swest[-1]],le):

  plon=grd.x[k,0]
  plat=grd.y[k,0]
  i=np.abs(lon-np.float(plon)).argmin()
  j=np.abs(lat-np.float(plat)).argmin()

  x0,y0 = lon[i], lat[j]

  # retrieve the 4 nearest diagonal points of the i,j
  lon4=lon[i-1:i+2:2]
  lon4=np.vstack([lon4,lon4])
  lat4=lat[j-1:j+2:2]
  lat4=np.vstack([lat4,lat4]).T
# print '==================='
  # define the quandrant
  A=plon-x0
  B=plat-y0
  for xx,yy in zip(lon4.flatten(),lat4.flatten()):
#        print xx,yy
         C=np.sign([xx-plon,A])
         D=np.sign([yy-plat,B])
         if C[0]==C[-1] and D[0] == D[-1] : 
           corner=[xx,yy]
           indx=[np.abs(lon-xx).argmin(),np.abs(lat-yy).argmin()]

  phv=np.zeros(ph.shape[-1])
  amv=np.zeros(amp.shape[-1])
  for k in range(amp.shape[-1]):

      p1=[x0,y0,amp[i,j,k]]
      p2=[x0,corner[1],amp[i,indx[1],k]]
      p3=[corner[0],corner[1],amp[indx[0],indx[1],k]]
      p4=[corner[0],y0,amp[indx[0],j,k]]
  
      points=[p1,p2,p3,p4]

      amv[k]= bilinear_interpolation(plon,plat,points)

      p1=[x0,y0,ph[i,j,k]]
      p2=[x0,corner[1],ph[i,indx[1],k]]
      p3=[corner[0],corner[1],ph[indx[0],indx[1],k]]
      p4=[corner[0],y0,ph[indx[0],j,k]]
  
      points=[p1,p2,p3,p4]

      phv[k]= bilinear_interpolation(plon,plat,points)

  f.write('West1{}\n'.format(l))
  for a,b,c in zip(tidal_c,amv,phv):
     f.write('{}         {:.7e}   {:.7e}\n'.format(a,b,c))


