import numpy as np
import datetime
import sys
from netCDF4 import Dataset


def pget(tstart,tend,path0,basename,point):

  tstart=datetime.datetime.strptime(tstart,'%Y%m%d.%H')
  tend=datetime.datetime.strptime(tend,'%Y%m%d.%H')
  dt=(tend-tstart).total_seconds()
  ndt=dt/(3600*12)
  ndt=np.int(ndt)+1

  combined=[]
  tw=[]

  for it in range(ndt-1): 
    idate=tstart+datetime.timedelta(hours=12*it)
    path=path0+'{}/{}/{:02d}/'.format(idate.month,idate.day,idate.hour)

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
  path=path0+'{}/{}/{:02d}/'.format(idate.month,idate.day,idate.hour)

   #print path
  d = Dataset(path+'trih-'+basename+'.nc')

  h=d.variables['ZWL'][:,point]
  time=d.variables['time'][:]

  tm=(idate-tstart).total_seconds()+(time-time[0])

  tstamp=[]
  for l in tm : tstamp.append(tstart+datetime.timedelta(0,int(l)))

  tcw=np.append(tcw,np.array(tstamp))
  cw=np.append(cw,np.array(h))

  return tcw,cw  

if __name__ == "__main__":
    tstart=sys.argv[1]
    tend=sys.argv[2]
    basename=sys.argv[3]
    point =sys.argv[4]
    t,ha=pget(tstart,tend,basename,point)

