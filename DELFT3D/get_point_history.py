import numpy as np
import datetime
import sys
from netCDF4 import Dataset
import pandas as pd


def pget(t0,t1,path0,basename,point):

  tstart=datetime.datetime.strptime(t0,'%Y%m%d.%H')
  tend=datetime.datetime.strptime(t1,'%Y%m%d.%H')
  dt=(tend-tstart).total_seconds()
  ndt=dt/(3600*12)
  ndt=np.int(ndt)+1

  path=path0+'{}/'.format(t0)

  obs=pd.read_csv(path+basename+'.obs',delim_whitespace=True,header=None, names=['ID','i','j'])
#  obs=obs.set_index(['ID'])
  ind,i,j=obs.xs(point)
  d = Dataset(path+'trim-'+basename+'.nc')
  xz=d['XZ'][i-1,j-1] # fortran/python conversion
  yz=d['YZ'][i-1,j-1] # fortran/python conversion

  combined=[]
  tw=[]

  for it in range(ndt-1): 
    idate=tstart+datetime.timedelta(hours=12*it)
    idt=datetime.datetime.strftime(idate,'%Y%m%d.%H')
    path=path0+'{}/'.format(idt)

   #print path
    d = Dataset(path+'trih-'+basename+'.nc')

    # 12 hours / dt 
    k=12*60/1

    h=d.variables['ZWL'][:k,point]
    time=d.variables['time'][:k]

    tm=(idate-tstart).total_seconds()+(time-time[0])

    tstamp=[]
    for l in tm : tstamp.append(tstart+datetime.timedelta(0,int(l)))

    tw.append(tstamp)
    combined.append(h)

# Add the last one with forecasting

  tcw=np.array(tw).flatten()
  cw=np.array(combined).flatten()

  idate=tend
  idt=datetime.datetime.strftime(idate,'%Y%m%d.%H')
  path=path0+'{}/'.format(idt)

   #print path
  d = Dataset(path+'trih-'+basename+'.nc')

  h=d.variables['ZWL'][:,point]
# hlon=d.variables['XSTAT'][:,point]
# hlat=d.variables['YSTAT'][:,point]
  time=d.variables['time'][:]

  tm=(idate-tstart).total_seconds()+(time-time[0])

  tstamp=[]
  for l in tm : tstamp.append(tstart+datetime.timedelta(0,int(l)))

  tcw=np.append(tcw,np.array(tstamp))
  cw=np.append(cw,np.array(h))

  return tcw,cw,yz,xz  

if __name__ == "__main__":
    tstart=sys.argv[1]
    tend=sys.argv[2]
    path0=sys.argv[3]
    basename=sys.argv[4]
    point =sys.argv[5]
    t,ha=pget(tstart,tend,path0,basename,point)

