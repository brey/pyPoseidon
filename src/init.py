import numpy as np
import pyresample
from netCDF4 import Dataset
from grid import *
import mdf
import os


def hinit(PATH1,PATH2,basename,nt):
  # previous run data
  d = Dataset(PATH1+'/trim-'+basename+'.nc')


  hz = d.variables['S1'][13,:,:]
  uw = d.variables['U1'][13,:,:]
  vw = d.variables['V1'][13,:,:]
  xz = d.variables['XZ'][:]
  yz = d.variables['YZ'][:]

#read grid
  grd=Grid.fromfile(PATH1+'/'+basename+'.grd')
  lons=grd.x.T
  lats=grd.y.T
  dx=lons[1,0]-lons[0,0]
  dy=lons[0,1]-lats[0,0]

# mask the values
  xzm=np.ma.masked_where(xz==0,xz)
  yzm=np.ma.masked_where(yz==0,yz)
  hzm=np.ma.masked_where(hz==0,hz)
  uwm=np.ma.masked_array(uw,hz==0)
  vwm=np.ma.masked_array(vw,hz==0)

#set up resampling 
  pz=pyresample.geometry.SwathDefinition(lons=xzm,lats=yzm)
  pu=pyresample.geometry.SwathDefinition(lons=xzm+dx/2.,lats=yzm)
  pv=pyresample.geometry.SwathDefinition(lons=xzm,lats=yzm+dy/2.)


#run current run for one minute in order to create output
  inp, order = mdf.read(PATH2+'/'+basename+'.mdf')

#read iteration start
  Tstart=inp['Tstart'][0]

  Tstop=Tstart+1.

#adjust iteration stop
  inp['Tstop']=[Tstop]
  
# current  run data
#adjust time for output
  step=1
  inp['Flmap']=[Tstart,step,Tstop]
  inp['Flhis']=[0,0,0]
  inp['Flpp']=[0,0,0]
  inp['Flrst']=[720]

#SAVE mdf file
  mdf.write(inp, PATH2+'/'+basename+'.mdf',selection=order)


#run interim 
  os.chdir(PATH2)
#subprocess.call(rpath+'run_flow2d3d.sh',shell=True)
  os.system('./run_flow2d3d.sh')


#current data
  d_ = Dataset(PATH2+'/trim-'+basename+'.nc')


  xz_ = d_.variables['XZ'][:]
  yz_ = d_.variables['YZ'][:]

#read grid
  grd_=Grid.fromfile(PATH2+'/'+basename+'.grd')
  lons_=grd_.x.T
  lats_=grd_.y.T
  dx_=lons_[1,0]-lons_[0,0]
  dy_=lons_[0,1]-lats_[0,0]


# mask the values
  xzm_=np.ma.masked_where(xz_==0,xz_)
  yzm_=np.ma.masked_where(yz_==0,yz_)

#set up resampling 
  cz=pyresample.geometry.SwathDefinition(lons=xzm_,lats=yzm_)
  cu=pyresample.geometry.SwathDefinition(lons=xzm_+dx_/2.,lats=yzm_)
  cv=pyresample.geometry.SwathDefinition(lons=xzm_,lats=yzm_+dy_/2.)


#resample
  hzm_=pyresample.kd_tree.resample_nearest(pz,hzm,cz,radius_of_influence=500000,fill_value=0)
  uwm_=pyresample.kd_tree.resample_nearest(pu,uwm,cu,radius_of_influence=500000,fill_value=0)
  vwm_=pyresample.kd_tree.resample_nearest(pv,vwm,cv,radius_of_influence=500000,fill_value=0)

# Write .ini file

  ini=np.zeros(hzm_.shape)     
  with open(PATH2+'/'+basename+'.ini', 'w') as f:
   np.savetxt(f,hzm_)
   np.savetxt(f,uwm_)
   np.savetxt(f,vwm_)

#update mdf for complete run
  inp, order = mdf.read(PATH2+'/'+basename+'.mdf')

#read iteration start
  Tstart=inp['Tstart'][0]

  Tstop=Tstart+60.*nt #minutes

#adjust iteration stop
  inp['Tstop']=[Tstop]
  
# current  run data
#adjust time for output
  step=60
  inp['Flmap']=[Tstart,step,Tstop]
  inp['Flhis']=[Tstart,1,Tstop]
  inp['Flpp']=[0,0,0]
  inp['Flrst']=[720]

# specify ini file
  if 'Filic' not in order: order.append('Filic')
  inp['Filic']=basename+'.ini'

  
#SAVE mdf file
  mdf.write(inp, PATH2+'/'+basename+'.mdf',selection=order)


if __name__ == "__main__":
    PATH1='/home/critechuser/TC/MATTHEW/20160929.00/'
    PATH2='/home/critechuser/TC/MATTHEW/20160929.12/'
    basename='MATTHEW'
    nt=72
    hinit(PATH1,PATH2,basename,nt)
