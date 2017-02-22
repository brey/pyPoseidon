import numpy as np
import datetime
import glob
import sys
from time import sleep
from copy import deepcopy
import os

import pyximport
pyximport.install()
from gribapi import *
from redtoreg import _redtoreg
from pygrib import gaulats
from tqdm import tqdm
import time


def gridd(lon1,lat1,lon2,lat2,nlats):

        #   lon1, lat1 = self.longitude_first_gridpoint, self.latitude_first_gridpoint
        #   lon2, lat2 = self.longitude_last_gridpoint, self.latitude_last_gridpoint
        #   nlats = self.points_in_y_direction
            # ECMWF 'reduced' gaussian grid.
            nlons = 2*nlats
            delon = 360./nlons
        #   lons = np.arange(lon1,lon2,delon)
            lons = np.linspace(lon1,lon2,nlons)
            # compute gaussian lats (north to south)
            lats = gaulats(nlats)
            if lat1 > lat2 :
               lats = lats[::-1]
          # lons = lons[::-1]
            lons,lats = np.meshgrid(lons,lats) # make 2-d arrays

            return lons,lats


def getd(f,t):
        gid = grib_new_from_file(f)#,headers_only = True)
        if gid is None: 
            print 'time = {}, gid = None'.format(t)
            sys.exit(1)

        name=grib_get(gid, 'shortName')
        mv=grib_get(gid,'missingValue')

        lonfgp=grib_get(gid,'longitudeOfFirstGridPointInDegrees')
        latfgp=grib_get(gid,'latitudeOfFirstGridPointInDegrees')
        lonlgp=grib_get(gid,'longitudeOfLastGridPointInDegrees')
        latlgp=grib_get(gid,'latitudeOfLastGridPointInDegrees')

        if grib_get(gid,'gridType') == 'regular_gg':

          Ni=grib_get(gid,'Ni')
          Nj=grib_get(gid,'Nj')
          lat=grib_get_array(gid,'latitudes')
          lat=lat.reshape(Nj,Ni)
          lat=np.flipud(lat)
          lon=grib_get_array(gid,'longitudes')
          lon=lon.reshape(Nj,Ni)

          values=grib_get_values(gid)
          dat=values.reshape(Nj,Ni)
          dat=np.flipud(dat)
          
        elif grib_get(gid,'gridType') == 'reduced_gg' :

          ss=grib_get_array(gid,'pl')  # lons per lat for the reduced_gg grid
          lon,lat = gridd(lonfgp,latfgp,lonlgp,latlgp,ss.size)

          values=grib_get_values(gid)
          ny=2*np.size(ss)

          dat=_redtoreg(ny,ss,values,mv)
          dat=np.flipud(dat)

        grib_release(gid)

        return name,dat,lon,lat



def wmap(date,nt1,nt2,minlon,maxlon,minlat,maxlat):

  yyyy=date.year # this is for facilitating the folder format below
  mm=date.month
  dd=date.day
  hh=date.hour
  # set PATH of the database.
  PATHbase=os.environ['ECMWF']  # Local location of the grib files
  PATH=PATHbase+'{:04d}/{:02d}/{:02d}/'.format(yyyy,mm,dd)

  dpath=glob.glob(PATH+'*{:04d}{:02d}{:02d}.{:02d}.tropical_cyclone.grib'.format(yyyy,mm,dd,hh))

  try: 
   f = open(dpath[0])
  except:
     print 'no file in {}'.format(PATH)
     return

  for it in xrange(nt1):
        gid = grib_new_from_file(f)#,headers_only = True)
        grib_release(gid)

  pt=[]
  ut=[]
  vt=[]

  mxv=nt2-nt1-1
  try:
    for it in tqdm(range(nt1,nt2)): # nt + the 0 hour

        name,varin,ilon,ilat=getd(f,it)        

        lon=ilon[0,:]
        lat=ilat[:,0]

    # get sea level pressure and 10-m wind data.
    # mult slp by 0.01 to put in units of hPa
        if name == 'msl' : varin=varin*.01

        if minlon < 0. :
           lon=lon-180.

           i1=np.abs(lon-minlon).argmin()-2
           i2=np.abs(lon-maxlon).argmin()+2
           j1=np.abs(lat-minlat).argmin()-2
           j2=np.abs(lat-maxlat).argmin()+2

           lons, lats = np.meshgrid(lon[i1:i2],lat[j1:j2])

           zlon=lon.shape[0]

           data1 = deepcopy(varin[j1:j2,zlon/2+i1:])
           data2 = deepcopy(varin[j1:j2,:i2-zlon/2])
           data = np.hstack([data1,data2])
          #data = np.hstack([data,varin[j1:j2,:i2-zlon/2]])

        else:

           i1=np.abs(lon-minlon).argmin()-2
           i2=np.abs(lon-maxlon).argmin()+2
           j1=np.abs(lat-minlat).argmin()-2
           j2=np.abs(lat-maxlat).argmin()+2

           lons, lats = np.meshgrid(lon[i1:i2],lat[j1:j2])
           data = deepcopy(varin[j1:j2,i1:i2])



    # mask the window

        if name == 'msl' : 
                 pt.append(data)
        elif name == '10u':
                 ut.append(data)
        elif name == '10v':
                 vt.append(data)


# END OF FOR
    #--------------------------------------------------------------------- 
    sys.stdout.flush()
    sys.stdout.write('\n')
    sys.stdout.write('meteo done\n')
    #--------------------------------------------------------------------- 

  except Exception as e:
    print e
    print 'ERROR in meteo input'

  f.close()

  return np.array(pt), np.array(ut), np.array(vt), lats, lons



if __name__ == "__main__":
  try:

    lon0=sys.argv[1]
    lon1=sys.argv[2]
    lat0=sys.argv[3]
    lat1=sys.argv[4]
    runtime=sys.argv[5]
    nt=sys.argv[6]
    
    to=datetime.datetime.strptime(runtime,'%Y%m%d.%H' )

  except:
    print "usage: meteo.py lon0 lon1 lat0 lat1 '20160101.00' 72"

  nt=3*(int(nt)+1)

  p,u,v,lat,lon = wmap(to,0,nt,float(lon0),float(lon1),float(lat0),float(lat1))


