import numpy as np
import datetime
import sys
from netCDF4 import Dataset
from grid import *
from dep import *
import itertools


def get(t0,t1,path0,basename,**kwargs):


  if 'lon0' in kwargs.keys():
     lonmin = kwargs['lon0']
  if 'lon1' in kwargs.keys():
     lonmax = kwargs['lon1']
  if 'lat0' in kwargs.keys():
     latmin = kwargs['lat0']
  if 'lat1' in kwargs.keys():
     latmax = kwargs['lat1']

  tstart=datetime.datetime.strptime(t0,'%Y%m%d.%H')
  tend=datetime.datetime.strptime(t1,'%Y%m%d.%H')
  dt=(tend-tstart).total_seconds()
  ndt=dt/(3600*12)
  ndt=np.int(ndt)+1

  path=path0+'{}/'.format(t0)

  # read computational grid
  grid = Grid.fromfile(path+basename+'.grd')
  lon=grid.x[0,:].data
  lat=grid.y[:,0].data
  bath=Dep.read(path+basename+'.dep',grid.shape)
  bath.val=bath.val.T[:-1,:-1]
  dx=lon[1]-lon[0]
  dy=lat[1]-lat[0]
  x=lon-dx/2.
  y=lat-dy/2.
  
  
  i0=0
  i1=lon.size
  j0=0
  j1=lat.size

  if 'lon0' in kwargs.keys(): i0=np.abs(lon-lonmin).argmin()
  if 'lon1' in kwargs.keys(): i1=np.abs(lon-lonmax).argmin()
  if 'lat0' in kwargs.keys(): j0=np.abs(lat-latmin).argmin()
  if 'lat1' in kwargs.keys(): j1=np.abs(lat-latmax).argmin()


  combined=[]
  tw=[]

  for it in range(ndt-1): 
    idate=tstart+datetime.timedelta(hours=12*it)
    idt=datetime.datetime.strftime(idate,'%Y%m%d.%H')
    path=path0+'{}/'.format(idt)

 #  print path
    d = Dataset(path+'trim-'+basename+'.nc')
    h=d.variables['S1'][:12,i0:i1,j0:j1]
    

    time=d.variables['time'][:12]

    tm=(idate-tstart).total_seconds()+(time-time[0])

    tstamp=[]
    for l in tm : tstamp.append(tstart+datetime.timedelta(0,int(l)))

    tw.append(tstamp)
    combined.append(h)#.mean(axis=(1,2)))

  tcw=np.array(tw).flatten()
  cw=np.array(combined).flatten()
# Add the last one with forecasting

  idate=tend
  idt=datetime.datetime.strftime(idate,'%Y%m%d.%H')
  path=path0+'{}/'.format(idt)
# print path

  d = Dataset(path+'trim-'+basename+'.nc')
  h=d.variables['S1'][:,i0:i1,j0:j1]
  xz = d.variables['XZ'][i0:i1,j0:j1]
  yz = d.variables['YZ'][i0:i1,j0:j1]
    

  time=d.variables['time'][:]

  tm=(idate-tstart).total_seconds()+(time-time[0])

  tstamp=[]
  for l in tm : tstamp.append(tstart+datetime.timedelta(0,int(l)))

  tcw=np.append(tcw,np.array(tstamp))
  cw=np.append(cw,np.array(h))#.mean(axis=(1,2))))
  
  return tcw,cw.reshape(tcw.shape[0],xz.shape[0],xz.shape[1]),x[i0:i1],y[j0:j1],xz,yz

if __name__ == "__main__":
    tstart=sys.argv[1]
    tend=sys.argv[2]
    path=sys.argv[3]
    basename=sys.argv[4]
    t,ha,x,y,zx,zy=get(tstart,tend,path,basename)

