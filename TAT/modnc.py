import numpy as np
from netCDF4 import Dataset
import glob
import datetime
import os
from grid import *
from dep import *
import matplotlib.pyplot as plt
import pandas
import copy
from osgeo import gdal,gdal_array
from locations import locfile
from writenc import writenc
import sys
import logging


from mpl_toolkits.basemap import Basemap, shiftgrid
from pmap import *

path0='/mnt/web/brey/2016H/'
basename='med'
SAVEPATH='D3/'
SAVEPATH='/mnt/ECMWF/processed/2016/FIX_MED_SEA_D3/'
RUNPATH='/home/critechuser/DELFT3D/python/'

#logging.basicConfig(filename='med.log',level=logging.INFO)

# start of the runs (usually 1st of the year)
t_ref=datetime.datetime.strptime('20160101.00','%Y%m%d.%H')


def fTAT(date=None):


  if not date:
   #--------------------------------------------------------------------
   # Evaluate the last run with 2 ways
   #--------------------------------------------------------------------
   # based on day
   flist=glob.glob(path0+'*/*/*')
   fflist=[x[-7:] for x in flist]
   fflist.sort()
   newest=fflist[-1]
   rtag=newest.split('/')
   rstamp=str(t_ref.year)+''.join(rtag)
  
  # based on time stamp
   #newest=max(glob.iglob(path0+'*/*/*'), key=os.path.getctime)
   #rtag=newest.split('/')[-4:]
   #rtag[0]=rtag[0][:4] # get only the year
   #rstamp=''.join(rtag)
  
   #rstamp='2016061712'
  #--------------------------------------------------------------------

  else:
   rstamp=date
  
  tend=datetime.datetime.strptime(rstamp,'%Y%m%d.%H')
  #tend=tend-datetime.timedelta(days=1)
  logging.info(tend)
  
  # 10 days back for the computation of the baseline
  tstart=tend-datetime.timedelta(days=10)
  
  
  # compute the number of rus within those days
  dt=(tend-tstart).total_seconds()
  ndt=dt/(3600*12)
  ndt=np.int(ndt)+1
  
  
  # read grid for lat/lon
  grid = Grid.fromfile(path0+'{}/{}/{:02d}/med.grd'.format(tend.month,tend.day,tend.hour))
  lon=grid.x[0,:].data
  lat=grid.y[:,0].data
  
  # read bathymetry
  deb=Dep.read(path0+'{}/{}/{:02d}/med.dep'.format(tend.month,tend.day,tend.hour),grid.shape)
  b=deb.val[:-1,:-1]
  w=np.isnan(b)
  w=b<0.
  
  
  
  # Variables to store h,u,v
  hcombined=[]
  ucombined=[]
  vcombined=[]
  tw=[]
  
  # read and store values from the first 12 hours forecasting
  for it in range(ndt): 
      idate=tstart+datetime.timedelta(hours=12*it)
      path=path0+'{}/{}/{:02d}/'.format(idate.month,idate.day,idate.hour)
  
      logging.info(path)
      try:
        d = Dataset(path+'trim-'+basename+'.nc')
        h=d.variables['S1'][:12,:,:]
        u=d.variables['WINDU'][:12,:,:]
        v=d.variables['WINDV'][:12,:,:]
        time=d.variables['time'][:12]
  
        tm=(idate-tstart).total_seconds()+(time-time[0])
        tstamp=[]
        for l in tm : tstamp.append(tstart+datetime.timedelta(0,int(l)))
      except:
        logging.info('no available data')
        continue
  
      tw.append(tstamp)
      hcombined.append(h)
      ucombined.append(u)
      vcombined.append(v)
  
  #create and flatten numpy arrays
  tcw=np.array(tw).flatten()
  cw=np.array(hcombined)
  uw=np.array(ucombined)
  vw=np.array(vcombined)
  
  #Read the forecasting beyond the 12 hours for the last run
  hf=d.variables['S1'][12:,:,:]
  uf=d.variables['WINDU'][12:,:,:]
  vf=d.variables['WINDV'][12:,:,:]
  timef=d.variables['time'][12:]
  
  
  # append to the corresponding variables
  tm=(idate-tstart).total_seconds()+(timef-time[0])
  tstamp=[]
  for l in tm : tstamp.append(tstart+datetime.timedelta(0,int(l)))
  
  tcw=np.append(tcw,tstamp)
  cw=np.append(cw,hf)
  uw=np.append(uw,uf)
  vw=np.append(vw,vf)
  
  # reshape properly
  cw=cw.reshape(tcw.size,lon.size,lat.size)
  uw=uw.reshape(tcw.size,lon.size,lat.size)
  vw=vw.reshape(tcw.size,lon.size,lat.size)
  
  # transpose prorerly
  ha=np.transpose(cw,axes=(0,2,1))
  ua=np.transpose(uw,axes=(0,2,1))
  va=np.transpose(vw,axes=(0,2,1))
  
  # adapt the shape (transpose)
  ts=[]
  
  for i in range(tcw.size):
     secs=(tcw[i]-t_ref).total_seconds()
     ts.append(secs)
  
  tsn=np.array(ts)


  # create the basemap 
  lon0=lon.min()
  lon1=lon.max()
  lat0=lat.min()
  lat1=lat.max()
    
  m = Basemap(projection='cyl',llcrnrlat=lat0,urcrnrlat=lat1,\
                 llcrnrlon=lon0,urcrnrlon=lon1,resolution='i')
    # grid
  xx,yy=np.meshgrid(lon,lat)


  # TIFF attributes
  dlat=lat[1]-lat[0]
  dlon=lon[1]-lon[0]
  geo=(lon.min(),dlon,0,lat.max(),0, -dlat)  
  #nodata=maxheight_ma.fill_value
  nband=1
    
  # read original full bathymetry
  deb0=Dep.read(path0+'{}/{}/{:02d}/med.dep'.format(t_ref.month,t_ref.day,t_ref.hour),grid.shape)
  b0=deb0.val[:-1,:-1]
    

  ## read places file ###
  points=pandas.read_csv(RUNPATH+'places.txt',delimiter='\t')
  plat,plon=points['lat'][:],points['long'][:]
    

  #@ FINAL FOLDER 
  if not date:
  
  ############# MAXHEIGHT ##################################################
  
    # store the max of all times
    maxheight=np.amax(ha,axis=0)
    maxheight_ma=np.ma.masked_where(w==True,maxheight)
    omin=maxheight_ma.min()
    omax=maxheight_ma.max()
    
    #fig=plt.figure(figsize=(btem.NCOLS/100.,btem.NROWS/100.)),dpi=100)
    fig=plt.figure()
    
    ma2=np.ma.masked_array(maxheight_ma,maxheight_ma>3.) # exclude outliers from bathymetry
    
    ma3=np.ma.masked_array(ma2,ma2<0.1)  # low threshold
    
    #clevs=np.linspace(-2.,3.5,12,endpoint=True)
    clevs=np.linspace(omin,omax,12,endpoint=True)
    #clevs=12
    
    CS2 = m.contourf(xx,yy,ma3,clevs,cmap=plt.cm.RdBu_r,animated=True)
    #cb = m.colorbar(CS2,"right", size="5%", pad="2%", ticks=clevs)
    m.drawcoastlines(linewidth=1.5)
    
    plt.tight_layout()
    
    plt.savefig(SAVEPATH+'final/P1_MAXHEIGHT_END.jpg',bbox_inches='tight', pad_inches=0)
    
    plt.show(block=False)
    
    plt.figure()
    
    CS3 = m.contourf(xx,yy,ma3,clevs,cmap=plt.cm.RdBu_r,animated=True)
    
    plt.tight_layout()
    
    plt.savefig(SAVEPATH+'final/P1_MAXHEIGHT_END.png',transparent=True,bbox_inches='tight', pad_inches=0)
    
    ###### LOCATIONS FILE ###################################################
    
    
    locfile(points,lat,lon,plat,plon,b,ha,tsn,tend,t_ref,SAVEPATH+'final/')
    
    # WRITE GTIFF
    
    # read default tif file
    tem=getmap(RUNPATH+'TIF_MAXHEIGHT_END.tif')
    
    # update values and save
    var=np.flipud(ma3.filled(tem.nan)) # save with the fill_value of tem
    filename=SAVEPATH+'final/TIF_MAXHEIGHT_END.tif'
    puttif(filename,var,geo,tem.Proj,nband,tem.nan)
    
    # read default tif file
    btem=getmap(RUNPATH+'bathymetry.tif')
    
    # update values and save
    var=np.flipud(b0)
    filename=SAVEPATH+'final/bathymetry.tif'
    puttif(filename,var,geo,btem.Proj,nband,btem.nan)
    
    
    # WRITE TO NETCDF
    
    writenc(SAVEPATH+'final/NETCDF_H.nc',lat,lon,tsn,ha)
    writenc(SAVEPATH+'final/NETCDF_U.nc',lat,lon,tsn,ua)
    writenc(SAVEPATH+'final/NETCDF_V.nc',lat,lon,tsn,va)

  #@ CALC_FOLDER 
  #### Folder of last computation #############################################
  
  fname='calc_{}'.format(datetime.datetime.strftime(tend,'%Y%m%d.%H'))
  if not os.path.exists(SAVEPATH+fname):
      os.makedirs(SAVEPATH+fname)
  
  # read default tif file
  tem=getmap(RUNPATH+'TIF_MAXHEIGHT_END.tif')
  
  hpr=ha[-73:,:,:]
  upr=ua[-73:,:,:]
  vpr=va[-73:,:,:]
  tpr=tsn[-73:]
  
  mh=np.amax(hpr,axis=0)
  
  mh0=hpr[0,:,:]
  mh0_ma=np.ma.masked_where(w==True,mh0)
  mh24=hpr[24,:,:]
  mh24_ma=np.ma.masked_where(w==True,mh24)
  mh48=hpr[48,:,:]
  mh48_ma=np.ma.masked_where(w==True,mh48)
  mh72=hpr[72,:,:]
  mh72_ma=np.ma.masked_where(w==True,mh72)
  
  mu0=upr[0,:,:]
  mu0_ma=np.ma.masked_where(w==True,mu0)
  mu0_ma[mu0_ma==0]=tem.nan
  mu24=upr[24,:,:]
  mu24_ma=np.ma.masked_where(w==True,mu24)
  mu24_ma[mu24_ma==0]=tem.nan
  mu48=upr[48,:,:]
  mu48_ma=np.ma.masked_where(w==True,mu48)
  mu48_ma[mu48_ma==0]=tem.nan
  mu72=upr[72,:,:]
  mu72_ma=np.ma.masked_where(w==True,mu72)
  mu72_ma[mu72_ma==0]=tem.nan
  
  mv0=vpr[0,:,:]
  mv0_ma=np.ma.masked_where(w==True,mv0)
  mv0_ma[mv0_ma==0]=tem.nan
  mv24=vpr[24,:,:]
  mv24_ma=np.ma.masked_where(w==True,mv24)
  mv24_ma[mv24_ma==0]=tem.nan
  mv48=vpr[48,:,:]
  mv48_ma=np.ma.masked_where(w==True,mv48)
  mv48_ma[mv48_ma==0]=tem.nan
  mv72=vpr[72,:,:]
  mv72_ma=np.ma.masked_where(w==True,mv72)
  mv72_ma[mv72_ma==0]=tem.nan
  
  
  plt.figure()
  
  maxh_ma=np.ma.masked_where(w==True,mh)
  maxh_ma2=np.ma.masked_array(maxh_ma,maxh_ma>3.)
  maxh_ma3=np.ma.masked_array(maxh_ma2,maxh_ma2<0.1)
  omin=maxh_ma.min()
  omax=maxh_ma.max()
  clevs=np.linspace(omin,omax,12,endpoint=True)
  CS3 = m.contourf(xx,yy,maxh_ma3,clevs,cmap=plt.cm.RdBu_r,animated=True)
  
  plt.tight_layout()
  
  plt.savefig(SAVEPATH+'{}/P1_MAXHEIGHT_END.png'.format(fname),transparent=True,bbox_inches='tight', pad_inches=0)
  
  var=np.flipud(maxh_ma3.filled(tem.nan))
  filename=SAVEPATH+'{}/TIF_MAXHEIGHT_END.tif'.format(fname)
  puttif(filename,var,geo,tem.Proj,nband,tem.nan)
  
  filename=SAVEPATH+'{}/TIF_H_{:08d}.tif'.format(fname,tpr[0].astype(int))
  puttif(filename,np.flipud(mh0_ma),geo,tem.Proj,nband,tem.nan)
  filename=SAVEPATH+'{}/TIF_H_{:08d}.tif'.format(fname,tpr[24].astype(int))
  puttif(filename,np.flipud(mh24_ma),geo,tem.Proj,nband,tem.nan)
  filename=SAVEPATH+'{}/TIF_H_{:08d}.tif'.format(fname,tpr[48].astype(int))
  puttif(filename,np.flipud(mh48_ma),geo,tem.Proj,nband,tem.nan)
  filename=SAVEPATH+'{}/TIF_H_{:08d}.tif'.format(fname,tpr[72].astype(int))
  puttif(filename,np.flipud(mh72_ma),geo,tem.Proj,nband,tem.nan)
  
  filename=SAVEPATH+'{}/u10x{:010d}.tif'.format(fname,tpr[0].astype(int))
  puttif(filename,np.flipud(mu0_ma),geo,tem.Proj,nband,tem.nan)
  filename=SAVEPATH+'{}/u10x{:010d}.tif'.format(fname,tpr[24].astype(int))
  puttif(filename,np.flipud(mu24_ma),geo,tem.Proj,nband,tem.nan)
  filename=SAVEPATH+'{}/u10x{:010d}.tif'.format(fname,tpr[48].astype(int))
  puttif(filename,np.flipud(mu48_ma),geo,tem.Proj,nband,tem.nan)
  filename=SAVEPATH+'{}/u10x{:010d}.tif'.format(fname,tpr[72].astype(int))
  puttif(filename,np.flipud(mu72_ma),geo,tem.Proj,nband,tem.nan)
  
  filename=SAVEPATH+'{}/u10y{:010d}.tif'.format(fname,tpr[0].astype(int))
  puttif(filename,np.flipud(mv0_ma),geo,tem.Proj,nband,tem.nan)
  filename=SAVEPATH+'{}/u10y{:010d}.tif'.format(fname,tpr[24].astype(int))
  puttif(filename,np.flipud(mv24_ma),geo,tem.Proj,nband,tem.nan)
  filename=SAVEPATH+'{}/u10y{:010d}.tif'.format(fname,tpr[48].astype(int))
  puttif(filename,np.flipud(mv48_ma),geo,tem.Proj,nband,tem.nan)
  filename=SAVEPATH+'{}/u10y{:010d}.tif'.format(fname,tpr[72].astype(int))
  puttif(filename,np.flipud(mv72_ma),geo,tem.Proj,nband,tem.nan)
  
  
  btem=getmap(RUNPATH+'bathymetry.tif')
  filename=SAVEPATH+'{}/bathymetry.tif'.format(fname)
  var=np.flipud(b0)
  puttif(filename,var,geo,btem.Proj,nband,btem.nan)
  
  
  writenc(SAVEPATH+'{}/NETCDF_H.nc'.format(fname),lat,lon,tpr,hpr)
  writenc(SAVEPATH+'{}/NETCDF_U.nc'.format(fname),lat,lon,tpr,upr)
  writenc(SAVEPATH+'{}/NETCDF_V.nc'.format(fname),lat,lon,tpr,vpr)
  
  
  # create locations file
  
  locfile(points,lat,lon,plat,plon,b,hpr,tpr,tend,t_ref,SAVEPATH+'{}/'.format(fname))
  
  #### Calc_input_deck file ################################################
  
  intime=(tend-t_ref).total_seconds()/3600.
  DateTsunami=datetime.datetime.strftime(t_ref,'%d %b %Y %H:%S')
  idCalc=datetime.datetime.strftime(tend,'%Y%m%d.%H')
  
  f=open(SAVEPATH+'{}/Calc_input_deck.txt'.format(fname),'wb')
  
  f.write('*********************************************\n')
  f.write('*       Tropical Cyclone input deck\n')
  f.write('*********************************************\n')
  f.write('* General data\n')
  f.write('Title=FIX_MED_SEA Storm Surge Calculation\n')
  f.write('*\n')
  f.write('DateTsunami={}   * bulDate is the date of the bulletin on which is based the time, i.e. date of fromBul\n'.format(DateTsunami))
  f.write('cyName=FIX_MED_SEA\n')
  f.write('idCyclone=FIX_MED_SEA\n')
  f.write('idCalc={}\n'.format(idCalc))
  f.write('*\n')
  f.write('*  Last position at time  *\n')
  f.write('Lat=               {}      * degree\n'.format(-1))
  f.write('Lon=               {}      * degree\n'.format(-1))
  f.write('Mag=               {}       * Cyclone s\n'.format(0))
  f.write('*  Calculation parameters  *\n')
  f.write('InTime=       {:.0f}.         * h\n'.format(intime))
  f.write('FinTime=      {:.0f}.        * h\n'.format(intime+72))
  f.write('Tsave=        {:.0f}.       * min\n'.format(60))
  f.write('fanning=0.015\n')
  f.write('*\n')
  f.write('faultform=-1    *  modifica     variab\n')
  f.write('*\n')
  f.write('*  Grid parameters  *\n')
  f.write('batgrid=      {:.0f}.        * min\n'.format(4))
  f.write('lonmin={:.4g}\n'.format(lon0))
  f.write('lonmax={:.4g}\n'.format(lon1))
  f.write('latmin={:.4g}\n'.format(lat0))
  f.write('latmax={:.4g}\n'.format(lat1))
  f.write('bathymetry=GEBCO30\n')
  
  f.close()
  
  return                  
