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
lons=grd.x[0,:]
lats=grd.y[:,0]

ba=Dep.read(RPATH+basename+'.dep',grd.shape)

#read from tide file
dmed=Dataset(PATH+'tpxo72.nc')

lat=dmed['lat'][:]
lon=dmed['lon'][:]

tidal_c=dmed['tidal_constituents'][:]

tidal_c=[''.join(k).upper().strip() for k in tidal_c]

amp=dmed['tidal_amplitude_h']
ph=dmed['tidal_phase_h']

#adjust lon according to the grid window
if lons[0]<0 :
    ii=np.abs(lon-180.).argmin()
    topo = lon[ii:]-360.
    lon = np.hstack([topo,lon[:ii]])
    topo = amp[ii:,:,:]
    amp = np.vstack([topo,amp[:ii,:,:]])
    topo = ph[ii:,:,:]
    ph = np.vstack([topo,ph[:ii,:,:]])


# strings to be used 
le=['A','B']

n=10  # points in chunks of the boundary

dic={'West':ba.val[:-1,0],'East':ba.val[:-1,-2],'South':ba.val[0,:-1],'North':ba.val[-2,:-1]}
idic={'West':[0,-99],'East':[lons.shape[0],-99],'South':[-99,0],'North':[-99,lats.shape[0]]}

#initiate output files
f = open(RPATH+basename+'.bnd', 'w')
f.close()
f = open(RPATH+basename+'.bca', 'w')
f.close()

def getboundary(bound,llon,llat):

#Identify WEST boundary 
   bv = dic[bound] 

   v=np.isfinite(bv).nonzero()
   branches = [list(x) for x in v]

#iterate over all branches
   for b in branches:

    s = b if b[0] not in [0,1] else b[1:]

# make chunks of n points if range is large
    chunks = [ s[i:i+n] for i in xrange(0, len(s),n if n < len(s) else len(s)) ]
    if len(chunks[-1]) < 3 : 
      chunks[-2]=chunks[-2]+chunks[-1] # merge the last chunk if too small
      chunks = chunks[:-1] # eliminate the last chunk

#write bnd file
    with open(RPATH+basename+'.bnd', 'a') as f:
      for ch in chunks:
        k1 = [ch[0] if x==-99 else x for x in idic[bound]] 
        k2 = [ch[-1] if x==-99 else x for x in idic[bound]] 
        f.write('{}{}                Z A     {}    {}     {}    {}   0.0000000e+00 {}{}A {}{}B\n'.format(bound,ch[0],k1[0]+1,k1[1]+1,k2[0]+1,k2[1]+1,bound,ch[0],bound,ch[0])) # fortran index ??


#write bca file

    with open(RPATH+basename+'.bca', 'a') as f:

     for ch in chunks:

      for k,l in zip([ch[0],ch[-1]],le):

        if l == 'A' : label = k

        plon=llon[k]
        plat=llat[k]
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

        f.write('{}{}{}\n'.format(bound,label,l))
        for a,b,c in zip(tidal_c,amv,phv):
          f.write('{}         {:.7e}   {:.7e}\n'.format(a,b,c))


getboundary('West',np.ones(lats.shape)*lons[0],lats)
getboundary('East',np.ones(lats.shape)*lons[-1],lats)
getboundary('North',lons,np.ones(lons.shape)*lats[-1])
getboundary('South',lons,np.ones(lons.shape)*lats[0])
