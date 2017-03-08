import numpy as np

from dep import *
from grid import *

from netCDF4 import Dataset

from scipy import interpolate

from bilinear import bilinear_interpolation

import sys

# strings to be used 
le=['A','B']

nm = ['Z', 'A']

def grouper(iterable):
    prev = None
    group = []
    for item in iterable:
        if not prev or item - prev <= 1:
            group.append(item)
        else:
            yield group
            group = [item]
        prev = item
    if group:
        yield group


def getboundary(bound,llon,llat,n,RPATH,basename):

#Identify WEST boundary 
   bv = dic[bound] 

   v=np.isfinite(bv).nonzero()
   branches = dict(enumerate(grouper(list(v[0])), 1))

   idx=0
#iterate over all branches
   for b in branches:

    if len(branches[b]) > 0:
     s = branches[b][:] if branches[b][0] not in [0,1] else branches[b][1:]

# make chunks of n points if range is large
     chunks = [ s[i:i+n] for i in xrange(0, len(s),n if n < len(s) else len(s)) ]
#    if len(chunks[-1]) < 3 : 
#     chunks[-2]=chunks[-2]+chunks[-1] # merge the last chunk if too small
#     chunks = chunks[:-1] # eliminate the last chunk

#write bnd file
     f1 = open(RPATH+basename+'.bnd', 'a') 
     f2 = open(RPATH+basename+'.bca', 'a')
     for ch in chunks:
        idx=idx+1
        k1 = [ch[0] if x==-99 else x for x in idic[bound]] 
        k2 = [ch[-1] if x==-99 else x for x in idic[bound]] 
        bname=bound+str(idx)
        f1.write('{0:<10s}{1:>12s}{2:>2s}{3:>6d}{4:>6d}{5:>6d}{6:>6d}   0.0000000e+00 {7:<s}{8:<g}A {9:<s}{10:<g}B\n'.format(bname,nm[0],nm[1],k1[0]+1,k1[1]+1,k2[0]+1,k2[1]+1,bound,idx,bound,idx)) # fortran index ??

#write bca file

        for k,l in zip([ch[0]-1,ch[-1]],le):

         if l == 'A' : label = idx

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

         f2.write('{}{}{}\n'.format(bound,label,l))
         for a,b,c in zip(tidal_c,amv,phv):
           f2.write('{0:<3s}        {1:<.7e}   {2:<.7e}\n'.format(a,b,c))


def tidebound(RPATH,basename,grd,ba,n):

  global dic, idic, lon, lat, ph, amp, tidal_c

  dlat=grd.y[1,0]-grd.y[0,0]
  dlon=grd.x[0,1]-grd.x[0,0]
  lons=grd.x[0,:]
  lats=grd.y[:,0]

  #read from tide file
  dmed=Dataset('../TIDES/tides.nc')

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



  dic={'West':ba.val[:-1,0],'East':ba.val[:-1,-2],'South':ba.val[0,:-1],'North':ba.val[-2,:-1]}
  idic={'West':[0,-99],'East':[lons.shape[0],-99],'South':[-99,0],'North':[-99,lats.shape[0]]}

  #initiate output files
  f = open(RPATH+basename+'.bnd', 'w')
  f.close()
  f = open(RPATH+basename+'.bca', 'w')
  f.close()

  getboundary('North',lons,np.ones(lons.shape)*lats[-1],n,RPATH,basename)
  getboundary('South',lons,np.ones(lons.shape)*lats[0],n,RPATH,basename)
  getboundary('West',np.ones(lats.shape)*lons[0],lats,n,RPATH,basename)
  getboundary('East',np.ones(lats.shape)*lons[-1],lats,n,RPATH,basename)



if __name__ == "__main__":
    path='/home/critechuser/test/20170227.00/'
    basename='tide'
    n=10  # points in chunks of the boundary
#   path=sys.argv[1]
#   basename=sys.argv[2]
#   n=np.int(sys.argv[3])
# read grd file
    grd=Grid.fromfile(path+basename+'.grd')
# read dep file
    ba=Dep.read(path+basename+'.dep',grd.shape)

    tidebound(path,basename,grd,ba,n)
