import numpy as np
import datetime
import sys
from shutil import copy2
import pickle

from meteo import wmap
from grid import *
from dep import *
#from dem import readgebco14
from dem import readem
from idelft3d import meteo2delft3d
import mdf
from setobs import createf
import glob
import os
import xml.dom.minidom as md


def setrun(lon0,lon1,lat0,lat1,basename,runtime,nt,resolution,path,force=False,**kwargs):

  resmin=resolution*60

  # computei ni,nj / correct lat/lon

  if 'lon' not in kwargs.keys() :

    ni=int((lon1-lon0)/resolution)+1
    nj=int((lat1-lat0)/resolution)+1
  
    lon1=lon0+ni*resolution
    lat1=lat0+nj*resolution

  else:

    lon=kwargs.lon
    lat=kwargs.lat

    nj,ni=lon.shape
    lon0=lon.min()
    lon1=lon.max()
    lat0=lat.min()
    lat1=lat.max()
   
  sys.stdout.write('\n')
  sys.stdout.write('run attributes lons=[{}, {}], lats=[{}, {}], ni={}, nj={}, resolution={} decimal degrees ({} minutes)\n'.format(lon0,lon1,lat0,lat1,ni,nj,resolution,resmin))
  sys.stdout.flush()
  sys.stdout.write('\n')
  
  # save attributes dic
  dic={}
  dic['bname']=basename
  dic['lon0']=lon0
  dic['lon1']=lon1
  dic['lat0']=lat0
  dic['lat1']=lat1
  dic['ni']=ni
  dic['nj']=nj
  dic['nt']=nt

  case=path.split('/')[-2]


  yyyy=runtime.year
  mm=runtime.month
  dd=runtime.day
  hh=runtime.hour

  Tstart=60.*hh #minutes

  Tstop=Tstart+60.*nt #minutes

# check if path exists
  if not os.path.exists(path):
    os.makedirs(path)

  with open(path+case+'.pkl', 'w') as f:
    pickle.dump(dic,f)


  rpath=datetime.datetime.strftime(runtime,'%Y%m%d.%H' )

  calc_dir=path+'{}/'.format(rpath)  #run subfolder  
  if not os.path.exists(calc_dir):
    os.makedirs(calc_dir)


# ffiles=glob.glob(calc_dir+'*')
# mfiles=[calc_dir+'u.amu',calc_dir+'v.amv',calc_dir+'p.amp']
# check=np.all([mf in ffiles for mf in mfiles])

  if force == 'True':

#   try: 
    p,u,v,elat,elon = wmap(yyyy,mm,dd,hh,0,3*(nt+1),lon0,lon1,lat0,lat1)
#   except Exception as e:
#     print e

# # Write meteo data
    dlat=elat[1,0]-elat[0,0]
    dlon=elon[0,1]-elon[0,0]
    mlat0=elat[0,0] 
    mlon0=elon[0,0] 

  # convert to DELFT3D ascii format
    meteo2delft3d(p,u,v,mlat0,mlon0,dlat,dlon,runtime,nt,path=calc_dir,curvi=False)
  

  sys.stdout.write('\n')
  sys.stdout.write('Set bathymetry')
  sys.stdout.flush()
  sys.stdout.write('\n')

  if 'lon' not in kwargs.keys() :
    sys.stdout.write('create grid')
    sys.stdout.flush()
    sys.stdout.write('\n')
  # set the grid 
    x=np.linspace(lon0,lon1,ni)
    y=np.linspace(lat0,lat1,nj)
    lon,lat=np.meshgrid(x,y)

  #  GET bathymetry interpolated onto lon,lat
  bathf='../BATHYMETRY/dem.nc'
  bat = readem(lat0,lat1,lon0,lon1,bathf,lon,lat,plot=False,interpolate=True)
  
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

  Dep.write(ba,calc_dir+basename+'.dep')
  
  sys.stdout.write('\n')
  sys.stdout.write('Bathymetry done')
  sys.stdout.flush()
  sys.stdout.write('\n')
  sys.stdout.write('Set grid')
  
  # Create the GRID file
  grd = Grid()
  
  grd.shape = bat.shape
  grd.x = lon
  grd.y = lat
  grd.properties = {'Coordinate System': 'Spherical', 'alfori': 0.0, 'xori': 0.0, 'yori': 0.0}
  
  grd.write(calc_dir+basename+'.grd')

  sys.stdout.write('\n')
  sys.stdout.write('Grid done')
  sys.stdout.flush()
  sys.stdout.write('\n')
  sys.stdout.write('Set enc')

  # Write .ini file

# ini=np.zeros(grd.x.shape)     
# np.savetxt(calc_dir+basename+'.ini',ini)
# with open(calc_dir+basename+'.ini', 'w') as f:
#   np.savetxt(f,ini)
#   np.savetxt(f,ini)

  
  # Write .enc file
  
  with open(calc_dir+basename+'.enc','w') as f:
      f.write('{:>5}{:>5}\n'.format(ni+1,1))  # add one like ddb
      f.write('{:>5}{:>5}\n'.format(ni+1,nj+1))
      f.write('{:>5}{:>5}\n'.format(1,nj+1))
      f.write('{:>5}{:>5}\n'.format(1,1))
      f.write('{:>5}{:>5}\n'.format(ni+1,1))
  
  f.close()

  sys.stdout.write('\n')
  sys.stdout.write('Enc done')
  sys.stdout.flush()
  sys.stdout.write('\n')
  sys.stdout.write('Set observation points')
  
  # Write .obs file
  createf(calc_dir,basename,lat0,lat1,lon0,lon1,grd,ba)
  
  sys.stdout.write('\n')
  sys.stdout.write('Observation points done')
  sys.stdout.flush()
  sys.stdout.write('\n')
  sys.stdout.write('Set mdf')
  
  # Define the mdf input file
  # first read the default
  inp, order = mdf.read('default.mdf')
  
  #define grid file
  inp['Filcco']=basename+'.grd'
  
  #define enc file
  inp['Filgrd']=basename+'.enc'
  
  #define dep file
  inp['Fildep']=basename+'.dep'
  
  #define obs file
  inp['Filsta']=basename+'.obs'
  
  
  # adjust ni,nj
  inp['MNKmax']=[ni+1,nj+1,1]  # add one like ddb
  
  # adjust iteration date
  inp['Itdate']=datetime.datetime.strftime(runtime.date(),'%Y-%m-%d')
  
  #set time unit
  inp['Tunit']='M'

  #adjust iteration start
  inp['Tstart']=[Tstart]
  
  #adjust iteration stop
  inp['Tstop']=[Tstop]
  
  #adjust time step
  inp['Dt']=[1.]
  
  #adjust time for output
  step=60
  inp['Flmap']=[Tstart,step,Tstop]
  inp['Flhis']=[Tstart,1,Tstop]
  inp['Flpp']=[0,0,0]
  inp['Flrst']=[720]
  
  #time interval to smooth the hydrodynamic boundary conditions
  inp['Tlfsmo']=[0.]

  # specify ini file
# if 'Filic' not in order: order.append('Filic')
# inp['Filic']=basename+'.ini'

  # netCDF output
  if 'FlNcdf' not in order: order.append('FlNcdf')
  inp['FlNcdf'] = 'map his'

  
  #SAVE mdf file
  mdf.write(inp, calc_dir+basename+'.mdf',selection=order)

  sys.stdout.write('\n')
  sys.stdout.write('mdf done')
  sys.stdout.flush()
  sys.stdout.write('\n')

  # edit and save config file
  copy2('config_d_hydro.xml',calc_dir+'config_d_hydro.xml')

  xml=md.parse(calc_dir+'config_d_hydro.xml')

  xml.getElementsByTagName('mdfFile')[0].firstChild.replaceWholeText(basename+'.mdf')

  with open(calc_dir+'config_d_hydro.xml','w') as f:
      xml.writexml(f)

  copy2('run_flow2d3d.sh',calc_dir+'run_flow2d3d.sh')
  
  sys.stdout.write('\n')
  sys.stdout.write('All done')
  sys.stdout.flush()
  sys.stdout.write('\n')
  
if __name__ == "__main__":
  try:

    lon0=np.float(sys.argv[1])
    lon1=np.float(sys.argv[2])
    lat0=np.float(sys.argv[3])
    lat1=np.float(sys.argv[4])
    basename=sys.argv[5]
    tstamp=sys.argv[6]
    runtime=datetime.datetime.strptime(tstamp,'%Y%m%d.%H') 
    nt=np.int(sys.argv[7])
    resolution=np.float(sys.argv[8])
    path=sys.argv[9]
    force=sys.argv[10]

    setrun(lon0,lon1,lat0,lat1,basename,runtime,nt,resolution,path,force)


  except Exception as e:
    print e
    print 'usage: python setup.py minlon, maxlon, minlat, maxlat, basename, date (YYYYMMDD.HH), number of forecasts, resolution(decimal degrees), path, compute uvp(T|F)'
    print "ex: python setup.py -5.5 47.5 28.5 48. 'med' '20160620.00'  72  .05  '../../../tmp/' True"
