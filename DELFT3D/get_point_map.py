import numpy as np
import datetime
import sys
from netCDF4 import Dataset
from grid import *
from dep import *


#path0='/home/critechuser/2016/'

def get(t0,t1,path0,basename,plat,plon):

  plat=float(plat)
  plon=float(plon)

  print plon,plat

  tstart=datetime.datetime.strptime(t0,'%Y%m%d.%H')
  tend=datetime.datetime.strptime(t1,'%Y%m%d.%H')
  dt=(tend-tstart).total_seconds()
  ndt=dt/(3600*12)
  ndt=np.int(ndt)+1


  path=path0+'{}/{}/{:02d}/'.format(tstart.month,tstart.day,tstart.hour)
  grid = Grid.fromfile(path+basename+'.grd')
  lon=grid.x[0,:].data
  lat=grid.y[:,0].data
  bath=Dep.read(path+basename+'.dep',grid.shape)
  bath.val=bath.val.T[:-1,:-1]


  i=np.abs(lon-plon).argmin()
  j=np.abs(lat-plat).argmin()
  print i,j,grid.x.T[i,j],grid.y.T[i,j],bath.val[i,j]
  if bath.val[i,j]<0. :
    lon4=grid.x.T[i-5:i+6,j-5:j+6]
    lat4=grid.y.T[i-5:i+6,j-5:j+6]
    print zip(lon4.flatten(),lat4.flatten())
    print '==================='
    for k,l in zip(lon4.flatten(),lat4.flatten()):
         print k,l
         i4=np.abs(lon-k).argmin()
         j4=np.abs(lat-l).argmin()
         print i4,j4,bath.val[i4,j4]
         if bath.val[i4,j4] > 0.:
             i=i4
             j=j4
             print i,j,grid.x.T[i,j],grid.y.T[i,j],bath.val[i,j]
             break



  combined=[]
  tw=[]

  for it in range(ndt-1): 
    idate=tstart+datetime.timedelta(hours=12*it)
    path=path0+'{}/{}/{:02d}/'.format(idate.month,idate.day,idate.hour)

    print path
    d = Dataset(path+'trim-'+basename+'.nc')
    h=d.variables['S1'][:12,i,j]
    

    time=d.variables['time'][:12]

    tm=(idate-tstart).total_seconds()+(time-time[0])

    tstamp=[]
    for l in tm : tstamp.append(tstart+datetime.timedelta(0,int(l)))

    tw.append(tstamp)
    combined.append(h)

  tcw=np.array(tw).flatten()
  cw=np.array(combined).flatten()
# Add the last one with forecasting

  idate=tend
  path=path0+'{}/{}/{:02d}/'.format(idate.month,idate.day,idate.hour)

  d = Dataset(path+'trim-'+basename+'.nc')
  h=d.variables['S1'][:,i,j]
    

  time=d.variables['time'][:]

  tm=(idate-tstart).total_seconds()+(time-time[0])

  tstamp=[]
  for l in tm : tstamp.append(tstart+datetime.timedelta(0,int(l)))

  tcw=np.append(tcw,np.array(tstamp))
  cw=np.append(cw,np.array(h))

# w1=tcw<datetime.datetime.strptime('20160228.00','%Y%m%d.%H')
# cwm = cw[w1].mean()
  
  return tcw,cw

if __name__ == "__main__":
    tstart=sys.argv[1]
    tend=sys.argv[2]
    basename=sys.argv[3]
    plat =sys.argv[4]
    plon =sys.argv[5]
    t,ha=get(tstart,tend,basename,plat,plon)

