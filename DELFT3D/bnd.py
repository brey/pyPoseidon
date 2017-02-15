import numpy as np

from dep import *
from grid import *

from netCDF4 import Dataset

from scipy import interpolate

from bilinear import bilinear_interpolation

PATH='../TIDES/'
RPATH='/home/critechuser/COAST/EUR/20131101.00/'
basename='eur'


grd=Grid.fromfile(RPATH+basename+'.grd')
dlat=grd.y[1,0]-grd.y[0,0]
dlon=grd.x[0,1]-grd.x[0,0]


ba=Dep.read(RPATH+basename+'.dep',grd.shape)

#read from tide file
dmed=Dataset(PATH+'tpxo72.nc')

lat=dmed['lat'][:]
lon=dmed['lon'][:]

tidal_c=dmed['tidal_constituents'][:]

tidal_c=[''.join(k).upper().strip() for k in tidal_c]

amp=dmed['tidal_amplitude_h']
ph=dmed['tidal_phase_h']

#Identify WEST boundary 

west=[]
for j in range(grd.y.shape[0]):
  if ba.val[j,0] > 0. : west.append(j)

# split into branches    
s0=[x for x in west if x-1 not in west]
s1=[x for x in west if x+1 not in west]

# strings to be used 
le=['A','B']


#iterate over all branches
for g,h in zip(s0,s1):

  swest = west[g-1:h+1] if g not in [0,1] else west[1:h+1]

# make chunks of n points if range is large
  n=10
  chunks = [ swest[i:i+n] for i in xrange(0, len(swest),n if n < len(swest) else len(swest)) ]
  if len(chunks[-1]) < 3 : 
    chunks[-2]=chunks[-2]+chunks[-1] # merge the last chunk if too small
    chunks = chunks[:-1] # eliminate the last chunk

#write bnd file
  with open(RPATH+basename+'.bnd', 'wa') as f:

   for k in chunks:
     f.write('West{}                Z A     1    {}     1    {}   0.0000000e+00 West{}A West{}B\n'.format(k[0],k[0]+1,k[-1]+1,k[0],k[0])) # fortran index ??


#write bca file

  with open(RPATH+basename+'.bca', 'wa') as f:

    for ch in chunks:

      for k,l in zip([ch[0],ch[-1]],le):

        if l == 'A' : label = k

        plon=grd.x[k,0]
        if plon < 0 : plon = plon + 360.
        plat=grd.y[k,0]
        i=np.abs(lon-np.float(plon)).argmin()
        j=np.abs(lat-np.float(plat)).argmin()

        xx = lon[i-1:i+2]
        yy = lat[j-1:j+2]


        phv=np.zeros(ph.shape[-1])
        amv=np.zeros(amp.shape[-1])
        for m in range(amp.shape[-1]):
     
           zz = amp[i-1:i+2,j-1:j+2,m]
           fa=interpolate.RectBivariateSpline(xx,yy,zz,kx=2,ky=2)
           amv[m]= fa(plon,plat)

           zz = ph[i-1:i+2,j-1:j+2,m]
           fa=interpolate.RectBivariateSpline(xx,yy,zz,kx=2,ky=2)
           phv[m]= fa(plon,plat)

        f.write('West{}{}\n'.format(label,l))
        for a,b,c in zip(tidal_c,amv,phv):
          f.write('{}         {:.7e}   {:.7e}\n'.format(a,b,c))


