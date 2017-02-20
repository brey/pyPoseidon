import numpy as np
from readbufr import read_bufr
import datetime
from setup import setrun
from init2  import hinit
from grid import *
import os


SPATH='/home/critechuser/REPOS/StormS/DELFT3D/'
nt=72
resolution=0.1

def setw(hurName,path,rpath,buf=5.):


 #tc=read_bufr(path+hurName)
  tc=read_bufr(hurName)

  for idate in sorted(tc)[2:3]:

      print idate
      runtime=datetime.datetime.strptime(idate,'%Y%m%d%H%M') 
      print runtime

      idata=tc[idate]

      tlons=idata.plons
      tlats=idata.plats


      lon0=tlons.min()-buf
      lon1=tlons.max()+buf

      lat0=tlats.min()-buf
      lat1=tlats.max()+buf

      previous_run=runtime-datetime.timedelta(hours=12)
      PATH1=rpath+datetime.datetime.strftime(previous_run,'%Y%m%d.%H')
      PATH2=rpath+datetime.datetime.strftime(runtime,'%Y%m%d.%H')

      print PATH1,PATH2

#read previous grid
      grd=Grid.fromfile(PATH1+'/'+hurName+'.grd')
      lons=grd.x.T
      lats=grd.y.T
      dx=lons[1,0]-lons[0,0]
      dy=lats[0,1]-lats[0,0]

# adapt boundaries
      dlon0=lons.min()-lon0
      nlon0=int(dlon0/dx)
      lon0=lons.min()-nlon0*dx

      dlon1=lons.max()-lon1
      nlon1=int(dlon1/dx)
      lon1=lons.max()-nlon1*dx

      dlat0=lats.min()-lat0
      nlat0=int(dlat0/dy)
      lat0=lats.min()-nlat0*dy

      dlat1=lats.max()-lat1
      nlat1=int(dlat1/dy)
      lat1=lats.max()-nlat1*dy

      xg=np.arange(lon0,lon1,dx) 
      yg=np.arange(lat0,lat1,dy) 
      
      lon,lat=np.meshgrid(xg,yg)

      setrun(lon0,lon1,lat0,lat1,hurName,runtime,nt,resolution,rpath,'False',lon=lon,lat=lat)

      hinit(PATH1,PATH2,hurName,nt,lon,lat)

      os.chdir(PATH2)
      os.system('./run_flow2d3d.sh')


if __name__ == "__main__":
    path='/mnt/ECMWF/bufr/2016/'
    rpath='/home/critechuser/TC/MATTHEW/'
   #rpath='/mnt/rmdisk/MATTHEW/'
    hurName='MATTHEW'
    setw(hurName,path,rpath)
