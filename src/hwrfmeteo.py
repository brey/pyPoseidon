import numpy as np
import matplotlib.pyplot as plt
import datetime
from mpl_toolkits.basemap import Basemap, shiftgrid
import pygrib
import string
import sys
import urllib2, urllib
import os
import glob

PATH='/mnt/poseidon/HWRF/TC/'


def createf(hurName):

   dfiles=glob.glob(PATH+'*')

   for ipath in sorted(dfiles)[:1]:
     print ipath

     gstorm=glob.glob(ipath+'/*hwrfprs.core*')

     pt=[]
     ut=[]
     vt=[]

     for ifile in sorted(gstorm)[:3]: #
       print ifile
       try:
        data=pygrib.open(ifile)
       except:
        print 'no file {}'.format(ifile)
        return

       slp=data[1] # for storm
       u10=data[706]
       v10=data[707]
   #slp=data[41] # for d123
   #u10=data[45]
   #v10=data[46]
    # read lats,lons
       latitudes=u10.latlons()[0][:,0]
       longitudes=u10.latlons()[1][0,:]
       print latitudes.min(),latitudes.max()
       print longitudes.min(),longitudes.max()
    # reverse latitudes so they go from south to north.
       lons, lats = np.meshgrid(longitudes,latitudes)
    # get sea level pressure and 10-m wind data.
    # mult slp by 0.01 to put in units of hPa
       p = 0.01*slp.values
       u = u10.values
       v = v10.values

# APPEND TO ARRAY
       pt.append(p)
       ut.append(u)
       vt.append(v)

# END OF FOR

   #np.savez(PATH+cyclons[ic]+'.meteo.d123.{}'.format(dstamp), p=np.array(pt), u=np.array(ut), v=np.array(vt), lon=lons, lat=lats)

     return pt,ut,vt

if __name__ == "__main__":
  hurName=sys.argv[1]
  p,u,v=createf(hurName)
