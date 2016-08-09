import numpy as np
import datetime
import sys

from meteo import wmap
from grid import *
from dep import *
#from dem import readgebco14
from dem import readem
from idelft3d import meteo2delft3d
import mdf
from setobs import createf


def  setrun(lon0,lon1,lat0,lat1,fname,runtime,path,ni,nj):

  nt=72
  
  Tstop=60.*nt #minutes

  try: 
    yyyy=runtime.year
    mm=runtime.month
    dd=runtime.day
    hh=runtime.hour
    p,u,v,elat,elon = wmap(yyyy,mm,dd,hh,0,3*(nt+1),lon0,lon1,lat0,lat1)
  except:
    print 'error'
    runtime=runtime-datetime.timedelta(hours=12)
    yyyy=runtime.year
    mm=runtime.month
    dd=runtime.day-1
    hh=runtime.hour
    print 'running with input from the day before {} '.format( datetime.datetime.strftime(runtime,"%Y-%m-%d %H:00") )
    p,u,v,elat,elon = wmap(yyyy,mm,dd,hh,0,3*(nt+1),lon0,lon1,lat0,lat1)

# # Write meteo data
  dlat=elat[1,0]-elat[0,0]
  dlon=elon[0,1]-elon[0,0]
  mlat0=elat[0,0] 
  mlon0=elon[0,0] 
  
  meteo2delft3d(p,u,v,mlat0,mlon0,dlat,dlon,runtime,nt,path=path,curvi=False)
  
  # set the grid 
  x=np.linspace(lon0,lon1,ni)
  y=np.linspace(lat0,lat1,nj)
  lon,lat=np.meshgrid(x,y)
  
  #  GET bathymetry interpolated onto lon,lat
  pathb='../BATHYMETRY/GLOBAL/gebco30_DELTARES.nc'
  bat = readem(lat0,lat1,lon0,lon1,pathb,lon,lat,True)
 #bat = readgebco(lat0,lat1,lon0,lon1,lon,lat,True)
 #blons,blats,bat = readgebco(lat0,lat1,lon0,lon1,lon,lat)
  
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

  ba.val[ba.val<0]=-999.  # mask all dry points

  Dep.write(ba,path+fname+'.dep')
  
  # Write .ini file

# ini=np.zeros(grd.x.shape)     
# np.savetxt(path+fname+'.ini',ini)
# with open(path+fname+'.ini', 'w') as f:
#   np.savetxt(f,ini)
#   np.savetxt(f,ini)

  
  # Write .enc file
  
  with open(path+fname+'.enc','w') as f:
      f.write('{:>5}{:>5}\n'.format(ni+1,1))  # add one like ddb
      f.write('{:>5}{:>5}\n'.format(ni+1,nj+1))
      f.write('{:>5}{:>5}\n'.format(1,nj+1))
      f.write('{:>5}{:>5}\n'.format(1,1))
      f.write('{:>5}{:>5}\n'.format(ni+1,1))
  
  f.close()
  
  # Write .obs file
  createf(path,fname,lat0,lat1,lon0,lon1,grd,ba)
  
  
  # Define the mdf input file
  # first read the default
  inp, order = mdf.read('default.mdf')
  
  #define grid file
  inp['Filcco']=fname+'.grd'
  
  #define enc file
  inp['Filgrd']=fname+'.enc'
  
  #define dep file
  inp['Fildep']=fname+'.dep'
  
  #define obs file
  inp['Filsta']=fname+'.obs'
  
  
  # adjust ni,nj
  inp['MNKmax']=[ni+1,nj+1,1]  # add one like ddb
  
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

  # specify ini file
# if 'Filic' not in order: order.append('Filic')
# inp['Filic']=fname+'.ini'

  # netCDF output
  if 'FlNcdf' not in order: order.append('FlNcdf')
  inp['FlNcdf'] = 'map his'

  
  #SAVE mdf file
  mdf.write(inp, path+fname+'.mdf',selection=order)
  
  
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
    fname=sys.argv[5]
    tstamp=sys.argv[6]
    runtime=datetime.datetime.strptime(tstamp,'%Y%m%d.%H') 
    path=sys.argv[7]
    ni=np.int(sys.argv[8])
    nj=np.int(sys.argv[9])
    setrun(lon0,lon1,lat0,lat1,fname,runtime,path,ni,nj)

  except Exception as e:
    print e
    print 'usage: python setup minlon, maxlon, minlat, maxlat, basename, date (YYYYMMDD.HH), path, ni, nj'
    print "ex: python setup -5.5 47.5 28.5 48. 'med' '20160620.00' '../../../tmp2/' 727 285"
