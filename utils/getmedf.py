import numpy as np
import matplotlib.pyplot as plt
import datetime
from mpl_toolkits.basemap import Basemap, shiftgrid
import glob
import string
#from shapely.geometry import Point, mapping
#from fiona import collection
import sys

import pyximport
pyximport.install()
from gribapi import *
from redtoreg import _redtoreg
from pygrib import gaulats


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
          
        elif grib_get(gid,'gridType') == 'reduced_gg' :

          ss=grib_get_array(gid,'pl')  # lons per lat for the reduced_gg grid
          lon,lat = gridd(lonfgp,latfgp,lonlgp,latlgp,ss.size)

          values=grib_get_values(gid)
          ny=2*np.size(ss)

          dat=_redtoreg(ny,ss,values,mv)

        grib_release(gid)

        return name,dat,lon,lat



def wmap(yyyy,mm,dd,hh,time,nt,lon0,lon1,lat0,lat1,ni,nj):
  hh=np.int(hh) # The variables are passed as string and needs to become integer to be used later
  yyyy=np.int(yyyy)
  mm=np.int(mm)
  dd=np.int(dd)
  # specify date to plot.
  date = datetime.datetime(yyyy,mm,dd,hh)

  # set PATH of the database.
  PATHbase="/mnt/ECMWF/grib/"  # Local mapping location for the above network drive
  PATH=PATHbase+"%04i/%02i/%02i/" % (yyyy,mm,dd)

  dpath=glob.glob(PATH+"*%04i%02i%02i.%02i.tropical_cyclone.grib" % (yyyy,mm,dd,hh))

  try: 
   f = open(dpath[0])
  except:
     print 'no file in {}'.format(dpath[0])
     return

  # make orthographic basemap.
  m = Basemap(projection='cyl',llcrnrlat=lat0,urcrnrlat=lat1,\
             llcrnrlon=lon0,urcrnrlon=lon1,resolution='l')

  pt=[]
  ut=[]
  vt=[]

  try:
    for it in range(nt): # nt + the 0 hour
#       gid = grib_new_from_file(f)
#       if gid is None: 
#           print 'time = {}, gid = None'.format(t)
#           break
#       t=t+1

#       name=grib_get(gid, 'shortName')

#       Ni=grib_get(gid,'Ni')
#       Nj=grib_get(gid,'Nj')
#       lat=grib_get_array(gid,'latitudes')
#       lat=lat.reshape(Nj,Ni)
#       lat=np.flipud(lat)
#       lon=grib_get_array(gid,'longitudes')
#       lon=lon.reshape(Nj,Ni)

#       values=grib_get_values(gid)
#       var=np.flipud(values.reshape(Nj,Ni))

#       grib_release(gid)

        name,varin,lon,lat=getd(f,it)        
        if name == 'msl' : varin=varin*.01
        print it

    # transform lon to -180:180
        
        longitudes=lon[0,:]
        latitudes=lat[:,0]
    # add cyclic points manually (could use addcyclic function)
        var= np.zeros((varin.shape[0],varin.shape[1]+1),np.float64)
        var[:,0:-1] = varin[::-1]; var[:,-1] = varin[::-1,0]
        longitudes=np.append(longitudes,360.)
    # first, shift grid so it goes from -180 to 180 (instead of 0 to 360)
        g,newlons = shiftgrid(180.,var,longitudes,start=False)
        x,y=np.meshgrid(newlons,latitudes)


    # mask the window
#       lons=xx[0,:]
#       lats=yy[:,0]
#       w1=(lats > lat0) & (lats < lat1) # mask latitudes
#       w2=(lons > lon0) & (lons < lon1) # mask longitudes
#       i=np.sum(w1)
#       j=np.sum(w2)

#       aw1=(yy > lat0) & (yy < lat1) # mask latitudes
#       aw2=(xx > lon0) & (xx < lon1) # mask longitudes
#       w=(aw1==True ) & (aw2==True)

#       v=g[w].reshape(i,j)


    # get sea level pressure and 10-m wind data.
    # mult slp by 0.01 to put in units of hPa
        if name == 'msl' : 
                      v,xx,yy = \
    m.transform_scalar(g,newlons,latitudes,ni,nj,returnxy=True,masked=True)
                      pt.append(v)
        elif name == '10u':
                      ugrid=g
        elif name == '10v':
                      vgrid=g
    # transform vectors to projection grid.
                      uproj,vproj,xx,yy = \
    m.transform_vector(ugrid,vgrid,newlons,latitudes,ni,nj,returnxy=True,masked=True)
                      ut.append(uproj)
                      vt.append(vproj)

# END OF FOR

  except:
    print 'ERROR in meteo input'

  f.close()

# ln=xx[w].reshape(i,j)
# la=yy[w].reshape(i,j)

# return np.array(pt), np.array(ut), np.array(vt), ln, la
  return np.array(pt), np.array(ut), np.array(vt), xx, yy



if __name__ == "__main__":
  try:
    yyyy=sys.argv[1]
    mm=sys.argv[2]
    dd=sys.argv[3]
    hh=sys.argv[4]

    yyyy=int(yyyy)
    mm=int(mm)
    dd=int(dd)
    hh=int(hh)
  except:
    to=datetime.datetime.today() 
    yyyy=to.year 
    mm=to.month 
    dd=to.day-1
    hh=0

  lon0=-10.
  lat0=28.
  lon1=48.
  lat1=48.

  runtime=datetime.datetime(yyyy,mm,dd,hh)

  nt=3*73
  ni=700
  nj=260

  p,u,v,lon,lat = wmap(yyyy,mm,dd,hh,0,nt,lon0,lon1,lat0,lat1,ni,nj)


