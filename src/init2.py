import numpy as np
from netCDF4 import Dataset
from grid import *
import mdf
import os

def extendv(ar,ni,nj):
  nodata=np.zeros(nj)
  ar_1=np.vstack((ar,nodata))
  nodata=np.zeros((ni+1,1))
  ar_2=np.hstack((ar_1,nodata))
  return ar_2

def hinit(PATH1,PATH2,basename,nt,lons,lats):
  # previous run data
  d = Dataset(PATH1+'/trim-'+basename+'.nc')


  hz = d.variables['S1'][13,:,:].T
  uw = d.variables['U1'][13,0,:,:].T
  vw = d.variables['V1'][13,0,:,:].T
  xc = d.variables['XCOR'][:]
  yc = d.variables['YCOR'][:]

#read grid
  grd=Grid.fromfile(PATH1+'/'+basename+'.grd')
  plons=grd.x
  plats=grd.y
  dx=plons[0,1]-plons[0,0]
  dy=plats[1,0]-plats[0,0]

# roundup the grid for creating the map

  lons=np.around(lons,5)
  lats=np.around(lats,5)
  plons=np.around(plons,5)
  plats=np.around(plats,5)

# create the .ini file

  ni,nj=lons.shape

  hz_=np.zeros(lons.shape) 
  uw_=np.zeros(lons.shape) 
  vw_=np.zeros(lons.shape) 

  mask1=lons>=plons.min()
  mask2=lons<=plons.max()
  masklon=mask1 & mask2

  pmask1=plons>=lons.min()
  pmask2=plons<=lons.max()
  pmasklon=pmask1 & pmask2
  
  mask1=lats>=plats.min()
  mask2=lats<=plats.max()
  masklat=mask1 & mask2

  pmask1=plats>=lats.min()
  pmask2=plats<=lats.max()
  pmasklat=pmask1 & pmask2
  
# mask the common window
  hzm=np.ma.masked_array(hz[:-1,:-1],mask=np.logical_not(pmasklon & pmasklat))
  hzm_=np.ma.masked_array(hz_,mask=np.logical_not(masklon & masklat))
  uwm=np.ma.masked_array(uw[:-1,:-1],mask=np.logical_not(pmasklon & pmasklat))
  uwm_=np.ma.masked_array(uw_,mask=np.logical_not(masklon & masklat))
  vwm=np.ma.masked_array(vw[:-1,:-1],mask=np.logical_not(pmasklon & pmasklat))
  vwm_=np.ma.masked_array(vw_,mask=np.logical_not(masklon & masklat))
# copy values
  hzm_[np.logical_not(hzm_.mask)]=hzm[np.logical_not(hzm.mask)]
  uwm_[np.logical_not(uwm_.mask)]=uwm[np.logical_not(uwm.mask)]
  vwm_[np.logical_not(vwm_.mask)]=vwm[np.logical_not(vwm.mask)]

# add a line/column for compatibility
  fz=extendv(hz_,ni,nj)
  fu=extendv(uw_,ni,nj)
  fv=extendv(vw_,ni,nj)
# Write .ini file

  ini=np.zeros(hz_.shape)     
  with open(PATH2+'/'+basename+'.ini', 'w') as f:
   np.savetxt(f,fz)
   np.savetxt(f,fu)
   np.savetxt(f,fv)

#update mdf for complete run
  inp, order = mdf.read(PATH2+'/'+basename+'.mdf')

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
