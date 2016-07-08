import numpy as np
import datetime
import sys
from netCDF4 import Dataset
from pmap import getmap
import lxml.etree as et

##########################################################################################
  #get hyflux data
##########################################################################################

def hget(t0,t1,basename,plat,plon):

  tinit=datetime.datetime(2016,1,1,0,0)


  plat=float(plat)
  plon=float(plon)

  hfpath='/mnt/ECMWF/processed/2016/FIX_MED_SEA/'
  tstart=datetime.datetime.strptime(t0,'%Y%m%d.%H')
  tend=datetime.datetime.strptime(t1,'%Y%m%d.%H')
  dt=(tend-tstart).total_seconds()
  ndt=dt/(3600*12)
  ndt=np.int(ndt)+1

  filename=hfpath+'calc_{}/bathymetry.tif'.format(t0)
  grid = getmap(filename)

  gt=grid.GeoTr

  width=grid.NCOLS
  height=grid.NROWS

  minx = gt[0]
  miny = gt[3] + width*gt[4] + height*gt[5] 
  maxx = gt[0] + width*gt[1] + height*gt[2]
  maxy = gt[3] 

  lon=np.linspace(minx,maxx,width,endpoint=True)
  lat=np.linspace(miny,maxy,height,endpoint=True)

  i=np.abs(lat-plat).argmin()
  j=np.abs(lon-plon).argmin()


  combined=[]
  tw=[]

  for it in range(ndt-1): 
    idate=tstart+datetime.timedelta(hours=12*it)
    dstamp=datetime.datetime.strftime(idate,'%Y%m%d.%H')
    print dstamp
    try:
       tree0=et.ElementTree(file=hfpath+'calc_{}/locations.xml'.format(dstamp))
       for elem in list(tree0.getiterator()):
         if elem.tag == 'pubDate': 
                 tn=datetime.datetime.strptime(elem.text,'%d %b %Y %H:%M')
                 break
    except:
       print 'problem with locations.xml file'
       pass  

    if tn != tinit :
         print 'restart date :',tn
         tinit=tn

    try:

      filename=hfpath+'calc_{}/NETCDF_H.nc'.format(dstamp)

      d =  Dataset(filename)
      ha=d.variables['HA'][:19,i,j]
      t=d.variables['TIME'][:19]

      tp=[]
#     print tinit,t[0]
      for ts in t:
        if t[0] == 0 :
          tp.append(idate+datetime.timedelta(seconds=ts))
          tinit=idate
        else:
          tp.append(tinit+datetime.timedelta(seconds=ts))

      iw1=np.argwhere(np.array(tp)==idate).flatten()[0]
      iw2=np.argwhere(np.array(tp)==idate+datetime.timedelta(hours=12)).flatten()[0]
      tw.append(tp[iw1:iw2+1])
      combined.append(ha[iw1:iw2+1])

    except:
      print '{} not available, skiping'.format(dstamp)

  htcw=np.array(tw).flatten()
  hcw=np.array(combined).flatten()

# Add the last one with forecasting
  idate=tend
  dstamp=datetime.datetime.strftime(idate,'%Y%m%d.%H')
  print dstamp
  try:
       tree0=et.ElementTree(file=hfpath+'calc_{}/locations.xml'.format(dstamp))
       for elem in list(tree0.getiterator()):
         if elem.tag == 'pubDate': 
                 tn=datetime.datetime.strptime(elem.text,'%d %b %Y %H:%M')
                 break
  except:
       print 'problem with locations.xml file'
       pass  

  if tn != tinit :
         print 'restart date :',tn
         tinit=tn

  try:

      filename=hfpath+'calc_{}/NETCDF_H.nc'.format(dstamp)

      d =  Dataset(filename)
      ha=d.variables['HA'][:,i,j]
      t=d.variables['TIME'][:]

      tp=[]
#     print tinit,t[0]
      for ts in t:
        if t[0] == 0 :
          tp.append(idate+datetime.timedelta(seconds=ts))
          tinit=idate
        else:
          tp.append(tinit+datetime.timedelta(seconds=ts))

      iw1=np.argwhere(np.array(tp)==idate).flatten()[0]
      iw2=np.argwhere(np.array(tp)==idate+datetime.timedelta(hours=72)).flatten()[0]
      tw.append(tp[iw1:iw2+1])
      combined.append(ha[iw1:iw2+1])

  except:
      print '{} not available, skiping'.format(dstamp)

  htcw=np.append(htcw,np.array(tp[iw1:iw2+1]))
  hcw=np.append(hcw,np.array(ha[iw1:iw2+1]))

  return  htcw,hcw

if __name__ == "__main__":
    tstart=sys.argv[1]
    tend=sys.argv[2]
    basename=sys.argv[3]
    plat =sys.argv[4]
    plon =sys.argv[5]
    t,ha=hget(tstart,tend,basename,plat,plon)

