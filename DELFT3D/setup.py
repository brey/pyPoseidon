import numpy as np
import datetime

from mpl_toolkits.basemap import Basemap, shiftgrid

import sys

from meteo import *
from grid import *
from dep import *
from dem import *
from idelft3d import meteo2delft3d
import mdf
import collections
from setobs import *


def  setrun(lon0,lon1,lat0,lat1,fname,to,path,ni,nj):

  yyyy=to.year
  mm=to.month
  dd=to.day
  hh=to.hour # 0 or 12
  time=0
  
  nt=72
  
  runtime=datetime.datetime(yyyy,mm,dd,hh)
  
  Tstop=60.*nt #minutes

  try: 
    p,u,v,lat,lon = wmap(yyyy,mm,dd,hh,0,3*(nt+1),lon0,lon1,lat0,lat1)
  except:
    dd=dd-1
    hh=12
    runtime=datetime.datetime(yyyy,mm,dd,hh)
    print 'running with input from the day before {} '.format( datetime.datetime.strftime(runtime,"%Y-%m-%d %H:00") )
    p,u,v,lat,lon = wmap(yyyy,mm,dd,hh,0,3*(nt+1),lon0,lon1,lat0,lat1)

  # Write meteo data
  dlat=lat[1,0]-lat[0,0]
  dlon=lon[0,1]-lon[0,0]
  mlat0=lat[0,0] 
  mlon0=lon[0,0] 
  
  meteo2delft3d(p,u,v,mlat0,mlon0,dlat,dlon,runtime,nt,path=path,curvi=False)
  
  # set the grid 
  x=np.linspace(lon0,lon1,ni)
  y=np.linspace(lat0,lat1,nj)
  lon,lat=np.meshgrid(x,y)
  
  #  GET bathymetry interpolated onto lon,lat
  bat = read_mod_gebco(lat0,lat1,lon0,lon1,lon,lat,True)
  
  # Create the GRID file
  grd = Grid()
  
  grd.shape = bat.shape
  grd.x = lon
  grd.y = lat
  grd.properties = {'Coordinate System': 'Spherical', 'alfori': 0.0, 'xori': 0.0, 'yori': 0.0}
  
  grd.write(path+fname+'.grd')
  
  # Write bathymetry file
  ba = Dep()
  # append the line/column of nodata 
  nodata=np.empty(ni)
  nodata.fill(np.nan)
  bat1=np.vstack((bat,nodata))
  nodata=np.empty((nj+1,1))
  nodata.fill(np.nan)
  bat2=np.hstack((bat1,nodata))
  ba.val = -bat2
  ba.shape = bat2.shape
  
  Dep.write(ba,path+fname+'.dep')
  
  
  # Write .enc file
  
  with open(path+fname+'.enc','w') as f:
      f.write('{:>5}{:>5}\n'.format(ni,1))
      f.write('{:>5}{:>5}\n'.format(ni,nj))
      f.write('{:>5}{:>5}\n'.format(1,nj))
      f.write('{:>5}{:>5}\n'.format(1,1))
      f.write('{:>5}{:>5}\n'.format(ni,1))
  
  f.close()
  
  # Write .obs file
  createf(path,fname,lat0,lat1,lon0,lon1,grd,ba)
  
  
  # Define the mdf input file
  # first read the default
  inp, ord = mdf.read('default.mdf')
  
  #define grid file
  inp['Filcco']=fname+'.grd'
  
  #define enc file
  inp['Filgrd']=fname+'.enc'
  
  #define dep file
  inp['Fildep']=fname+'.dep'
  
  #define obs file
  inp['Filsta']=fname+'.obs'
  
  
  # adjust ni,nj
  inp['MNKmax']=[ni,nj,1]
  
  # adjust iteration date
  inp['Itdate']=datetime.datetime.strftime(runtime.date(),'%Y-%m-%d')
  
  #set time unit
  inp['Tunit']='M'
  
  #adjust iteration stop
  inp['Tstop']=[Tstop]
  
  #adjust time step
  inp['Dt']=[1.]
  
  #adjust time for output
  step=60
  inp['Flmap']=[0,step,Tstop]
  inp['Flhis']=[0,1,Tstop]
  inp['Flpp']=[0,0,0]
  inp['Flrst']=[720]
  
  #time interval to smooth the hydrodynamic boundary conditions
  inp['Tlfsmo']=[0.]
  
  #SAVE mdf file
  mdf.write(inp, path+fname+'.mdf',selection=ord)
  
  
if __name__ == "__main__":
#ni=727
#nj=285
#lon0=-5.5
#lat0=28.5
#lon1=43.
#lat1=47.5
  try:

    lon0=np.float(sys.argv[1])
    lon1=np.float(sys.argv[2])
    lat0=np.float(sys.argv[3])
    lat1=np.float(sys.argv[4])
    basename=sys.argv[5]
    tstamp=sys.argv[6]
    date=datetime.datetime.strptime(tstamp,'%Y%m%d.%H') 
    path=sys.argv[7]
    ni=np.int(sys.argv[8])
    nj=np.int(sys.argv[9])
    setrun(lon0,lon1,lat0,lat1,basename,date,path,ni,nj)

  except:
    print 'usage: python setup minlon, maxlon, minlat, maxlat, basename, date (YYYYMMDD.HH), path, ni, nj'
