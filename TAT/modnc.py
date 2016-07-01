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


from mpl_toolkits.basemap import Basemap, shiftgrid
import lxml.etree as et
from pmap import *


path0='/mnt/web/brey/2016H/'
path0='/mnt/web/brey/venice/2016/'
basename='med'
SAVEPATH='D3/'
COPYPATH='/mnt/ECMWF/processed/2016/FIX_MED_SEA_D3/'

t_ref=datetime.datetime.strptime('2016060100','%Y%m%d%H')

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

tend=datetime.datetime.strptime(rstamp,'%Y%m%d%H')
#tend=tend-datetime.timedelta(days=1)
print tend

tstart=tend-datetime.timedelta(days=10)

dt=(tend-tstart).total_seconds()
ndt=dt/(3600*12)
ndt=np.int(ndt)+1


grid = Grid.fromfile(path0+'{}/{}/{:02d}/med.grd'.format(t_ref.month,t_ref.day,t_ref.hour))
lon=grid.x[0,:].data
lat=grid.y[:,0].data

####### SEA LEVEL #########
deb=Dep.read(path0+'{}/{}/{:02d}/med.dep'.format(t_ref.month,t_ref.day,t_ref.hour),grid.shape)
b=deb.val[:-1,:-1]
w=b<0.




hcombined=[]
ucombined=[]
vcombined=[]
tw=[]

for it in range(ndt): 
    idate=tstart+datetime.timedelta(hours=12*it)
    path=path0+'{}/{}/{:02d}/'.format(idate.month,idate.day,idate.hour)

    print path
    d = Dataset(path+'trim-'+basename+'.nc')
    h=d.variables['S1'][:12,:,:]
    u=d.variables['WINDU'][:12,:,:]
    v=d.variables['WINDV'][:12,:,:]
    time=d.variables['time'][:12]

    tm=(idate-tstart).total_seconds()+(time-time[0])
    tstamp=[]
    for l in tm : tstamp.append(tstart+datetime.timedelta(0,int(l)))

    tw.append(tstamp)
    hcombined.append(h)
    ucombined.append(u)
    vcombined.append(v)

tcw=np.array(tw).flatten()
cw=np.array(hcombined)
uw=np.array(ucombined)
vw=np.array(vcombined)

hf=d.variables['S1'][12:,:,:]
uf=d.variables['WINDU'][12:,:,:]
vf=d.variables['WINDV'][12:,:,:]
timef=d.variables['time'][12:]

tm=(idate-tstart).total_seconds()+(timef-time[0])
tstamp=[]
for l in tm : tstamp.append(tstart+datetime.timedelta(0,int(l)))

tcw=np.append(tcw,tstamp)
cw=np.append(cw,hf)
uw=np.append(uw,uf)
vw=np.append(vw,vf)


cw=cw.reshape(tcw.size,lon.size,lat.size)
uw=uw.reshape(tcw.size,lon.size,lat.size)
vw=vw.reshape(tcw.size,lon.size,lat.size)

ha=np.zeros((tcw.size,lat.size,lon.size))
ua=np.zeros((tcw.size,lat.size,lon.size))
va=np.zeros((tcw.size,lat.size,lon.size))

ts=[]

bh=np.zeros(ha.shape)
for i in range(tcw.size):
   ha[i,:,:]=cw[i,:,:].T
   ua[i,:,:]=uw[i,:,:].T
   va[i,:,:]=vw[i,:,:].T
   secs=(tcw[i]-t_ref).total_seconds()
   ts.append(secs)
   bh[i,:,:]=np.ma.masked_where(w==True,ha[i,:,:])

tsn=np.array(ts)


# Max Height graph

lon0=lon.min()
lon1=lon.max()
lat0=lat.min()
lat1=lat.max()

m = Basemap(projection='cyl',llcrnrlat=lat0,urcrnrlat=lat1,\
             llcrnrlon=lon0,urcrnrlon=lon1,resolution='i')


maxheight=np.amax(bh,axis=0)
maxheight_ma=np.ma.masked_where(w==True,maxheight)
omin=maxheight_ma.min()
omax=maxheight_ma.max()

#############################################################################################################
#interpolate on HyFlux grid
#############################################################################################################
#btem=getmap('bathymetry.tif')

#h1=np.zeros((tcw.size,btem.NROWS,btem.NCOLS))
#u1=np.zeros((tcw.size,btem.NROWS,btem.NCOLS))
#v1=np.zeros((tcw.size,btem.NROWS,btem.NCOLS))
#
#for i in range(tcw.size):
#  h1[i,:,:],x1,y1 = \
#    m.transform_scalar(ha[i,:,:],lon,lat,btem.NCOLS,btem.NROWS,returnxy=True,masked=True)
#  u1[i,:,:],x1,y1 = \
#    m.transform_scalar(ua[i,:,:],lon,lat,btem.NCOLS,btem.NROWS,returnxy=True,masked=True)
#  v1[i,:,:],x1,y1 = \
#    m.transform_scalar(va[i,:,:],lon,lat,btem.NCOLS,btem.NROWS,returnxy=True,masked=True)
#
#
#lon=x1[0,:]
#lat=y1[:,0]
#ha=h1
#ua=u1
#va=v1
#b=np.flipud(btem.data)
#w=b<0.
#############################################################################################################


maxheight=np.amax(ha,axis=0)
maxheight_ma=np.ma.masked_where(w==True,maxheight)

#clevs=np.linspace(-2.,3.5,12,endpoint=True)
clevs=np.linspace(omin,omax,12,endpoint=True)

xx,yy=np.meshgrid(lon,lat)

#fig=plt.figure(figsize=(btem.NCOLS/100.,btem.NROWS/100.)),dpi=100)
fig=plt.figure()

ma2=np.ma.masked_array(maxheight_ma,maxheight_ma>3.2)

ma3=np.ma.masked_array(ma2,ma2<0.1)

#CS2 = m.contourf(xx,yy,maxheight_ma,12,cmap=plt.cm.RdBu_r,animated=True)
CS2 = m.contourf(xx,yy,ma3,clevs,cmap=plt.cm.RdBu_r,animated=True)
#cb = m.colorbar(CS2,"right", size="5%", pad="2%", ticks=clevs)
m.drawcoastlines(linewidth=1.5)

plt.tight_layout()

plt.savefig(SAVEPATH+'final/P1_MAXHEIGHT_END.jpg',bbox_inches='tight', pad_inches=0)

plt.show(block=False)

plt.figure()

CS3 = m.contourf(xx,yy,ma3,clevs,cmap=plt.cm.RdBu_r,animated=True)
#CS3 = m.contourf(xx,yy,maxheight_ma,12,cmap=plt.cm.RdBu_r,animated=True)

plt.tight_layout()

plt.savefig(SAVEPATH+'final/P1_MAXHEIGHT_END.png',transparent=True,bbox_inches='tight', pad_inches=0)

## locations.xml ###
points=pandas.read_csv('places.txt',delimiter='\t')
plat,plon=points['lat'][:],points['long'][:]

# read xml templates

tree0=et.ElementTree(file='loc0.xml')
root=tree0.getroot() #rss
channel=root[0] ##channel

channel[2][0].text=et.CDATA(channel[2][0].text)

for elem in list(tree0.getiterator()):
   if elem.tag == 'pubDate': elem.text=datetime.datetime.strftime(t_ref,'%d %b %Y %H:%S')

tree=et.ElementTree(file='loc.xml')
litem=tree.getroot()[0]

#print et.tostring(litem,pretty_print=True)


def modtxt(tag,val):
   for elem in nitem.findall(tag): elem.text=val

# APPEND LOCATIONS


for l in range(points.shape[0]):
   if (lat.min() < points['lat'][l] < lat.max()) & (lon.min() < points['long'][l] < lon.max()):
  #    print points.values[l]
       nitem=copy.deepcopy(litem)

       i=np.abs(lon-plon[l]).argmin()
       j=np.abs(lat-plat[l]).argmin()

 #     print l, i,j

       bd=b[i-1:i+2,j-1:j+2]

       bda=np.where(bd>20.)

       if np.size(bda) > 0 : 

          bloc=bd[bda].min()

          pi,pj=np.where(bd==bd[bda].min())

          ip,jp=i-1+pi.ravel()[0],j-1+pj.ravel()[0]

          tseries=ha[:,ip,jp]

          maxh=tseries.max()
          print maxh

          if maxh-tseries[0] < 0.05 : continue

          loc=np.argwhere(tseries==maxh)
          print loc

          dh=tseries-tseries[0]
          try: 
             arr=np.argwhere(dh>0.05)[0]
          except:
             arr=loc

          tMax=tsn[loc.ravel()[0]]/3600.
          tArr=tsn[arr.ravel()[0]]/3600.

#      print points.values[l]

          for elem in nitem.findall('title'): 
            elem.text=et.CDATA('{}: {} ({} m)'.format(points['$country'][l],points['$name'][l],maxh))

          for elem in nitem.findall('description'): 
           elem.text=et.CDATA('Country: {}\n \
Location: {}\n \
Time (hh:mm): {}\n \
Maximum Height: {} m\n \
Time Max (hh:mm): {}'.format(points['$country'][l],points['$name'][l],'{}:00'.format(tArr.astype(int)),maxh,'{}:00'.format(tMax.astype(int))))

          modtxt('pubDate',datetime.datetime.strftime(tend,'%d %b %Y %H:00:00'))

          modtxt('cityName',points['$name'][l])

          modtxt('country',points['$country'][l])

          modtxt('maxHeight',maxh.astype(str))

          modtxt('ID',points['id'][l].astype(str))

          modtxt('timeMaxH','{}:00'.format(tMax.astype(int)))

          modtxt('timeMaxH_value',tMax.astype(str))

          modtxt('timeArrival','{}:00'.format(tArr.astype(int)))

          modtxt('timeArrival_value',tArr.astype(str))

          modtxt('cityClass',points['cityclass'][l].astype(str))

          if not pandas.isnull(points['popest'][l]) :
              modtxt('popEst',points['popest'][l])
          else:
              modtxt('popEst','')

          modtxt('{http://www.w3.org/2003/01/geo/wgs84_pos#}long',lon[ip].astype(str))
          modtxt('{http://www.w3.org/2003/01/geo/wgs84_pos#}lat',lat[jp].astype(str))

          modtxt('depth',bloc.astype(str))

          channel.append(nitem)


tree0.write(SAVEPATH+'final/locations.xml',xml_declaration=True,encoding='utf-8',method='xml')

# WRITE GTIFF

dlat=lat[1]-lat[0]
dlon=lon[1]-lon[0]
geo=(lon.min(),dlon,0,lat.max(),0, -dlat)  
#nodata=maxheight_ma.fill_value
nband=1

tem=getmap('TIF_MAXHEIGHT_END.tif')

mplot=np.amax(bh,axis=0)
mz=mplot>3.5
mplot[mz]=tem.nan

mz=mplot<0.1
mplot[mz]=tem.nan

var=np.flipud(mplot)
filename=SAVEPATH+'final/TIF_MAXHEIGHT_END.tif'
puttif(filename,var,geo,tem.Proj,nband,tem.nan)

btem=getmap('bathymetry.tif')

var=np.flipud(b)
filename=SAVEPATH+'final/bathymetry.tif'
puttif(filename,var,geo,btem.Proj,nband,btem.nan)




# WRITE TO NETCDF

def writenc(filename,lat,lon,tstamp,val):

 ni=lon.shape[0]
 nj=lat.shape[0]

 rootgrp = Dataset(filename, 'w', format='NETCDF3_64BIT')
 lats = rootgrp.createDimension('LAT', nj)
 lons = rootgrp.createDimension('LON', ni)
 time = rootgrp.createDimension('TIME', None)


 longitudes = rootgrp.createVariable('LON','f8',('LON',))
 latitudes = rootgrp.createVariable('LAT','f8',('LAT',))
 times = rootgrp.createVariable('TIME','f8',('TIME',))
 levels = rootgrp.createVariable('HA','f8',('TIME','LAT','LON'))

 rootgrp.description = ''
 rootgrp.history = 'DELFT3D - JRC Ispra European Commission'
 rootgrp.source = 'netCDF4 python module tutorial'
 latitudes.units = 'degrees_north'
 latitudes.point_spacing = 'even'
 longitudes.units = 'degrees_east'
 longitudes.point_spacing = 'even'
 levels.units = 'm'
 times.units = 'seconds'

 
 levels[:]=val
 times[:]=tstamp
 latitudes[:]=lat
 longitudes[:]=lon

 rootgrp.close()



writenc(SAVEPATH+'final/NETCDF_H.nc',lat,lon,tsn,bh)
writenc(SAVEPATH+'final/NETCDF_U.nc',lat,lon,tsn,ua)
writenc(SAVEPATH+'final/NETCDF_V.nc',lat,lon,tsn,va)

# last computation

fname='calc_{}'.format(datetime.datetime.strftime(tend,'%Y%m%d.%H'))
if not os.path.exists(SAVEPATH+fname):
    os.makedirs(SAVEPATH+fname)


hpr=bh[-73:,:,:]
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




def out(var,outname):

  plt.figure()

  #CS2 = m.contourf(xx,yy,var,12,cmap=plt.cm.RdBu_r,animated=True)
  CS2 = m.contourf(xx,yy,var,clevs,cmap=plt.cm.RdBu_r,animated=True)
  #cb = m.colorbar(CS2,"right", size="5%", pad="2%", ticks=clevs)
  m.drawcoastlines(linewidth=1.5)

  plt.tight_layout()

  plt.savefig(SAVEPATH+'{}/{}.jpg'.format(fname,outname),bbox_inches='tight', pad_inches=0)



plt.figure()

maxh_ma=np.ma.masked_where(w==True,mh)
maxh_ma2=np.ma.masked_array(maxh_ma,maxh_ma>3.2)
maxh_ma3=np.ma.masked_array(maxh_ma2,maxh_ma2<0.1)
CS3 = m.contourf(xx,yy,maxh_ma3,clevs,cmap=plt.cm.RdBu_r,animated=True)
  #CS3 = m.contourf(xx,yy,var,12,cmap=plt.cm.RdBu_r,animated=True)

plt.tight_layout()

plt.savefig(SAVEPATH+'{}/P1_MAXHEIGHT_END.png'.format(fname),transparent=True,bbox_inches='tight', pad_inches=0)

mplot=np.amax(hpr,axis=0)
mplot_ma=np.ma.masked_where(w==True,mplot)
mz=mplot_ma>3.5
mplot[mz]=tem.nan
mz=mplot_ma<0.1
mplot[mz]=tem.nan

var=np.flipud(mplot)
filename=SAVEPATH+'{}/TIF_MAXHEIGHT_END.tif'.format(fname)
puttif(filename,var,geo,tem.Proj,nband,tem.nan)

filename=SAVEPATH+'{}/TIF_H_{}.tif'.format(fname,tpr[0].astype(int))
puttif(filename,np.flipud(mh0_ma),geo,tem.Proj,nband,tem.nan)
filename=SAVEPATH+'{}/TIF_H_{}.tif'.format(fname,tpr[24].astype(int))
puttif(filename,np.flipud(mh24_ma),geo,tem.Proj,nband,tem.nan)
filename=SAVEPATH+'{}/TIF_H_{}.tif'.format(fname,tpr[48].astype(int))
puttif(filename,np.flipud(mh48_ma),geo,tem.Proj,nband,tem.nan)
filename=SAVEPATH+'{}/TIF_H_{}.tif'.format(fname,tpr[72].astype(int))
puttif(filename,np.flipud(mh72),geo,tem.Proj,nband,tem.nan)

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

#out(mh0,'OUT_TIF_H_00')
#out(mh0,'OUT_TIF_H_00')
#out(mh0,'OUT_TIF_H_00')
#out(mh24,'OUT_TIF_H_24')
#out(mh48,'OUT_TIF_H_48')
#out(mh72,'OUT_TIF_H_72')

btem=getmap('bathymetry.tif')
filename=SAVEPATH+'{}/bathymetry.tif'.format(fname)
puttif(filename,var,geo,btem.Proj,nband,btem.nan)


writenc(SAVEPATH+'{}/NETCDF_H.nc'.format(fname),lat,lon,tpr,hpr)
writenc(SAVEPATH+'{}/NETCDF_U.nc'.format(fname),lat,lon,tpr,upr)
writenc(SAVEPATH+'{}/NETCDF_V.nc'.format(fname),lat,lon,tpr,vpr)




## locations.xml ###
points=pandas.read_csv('places.txt',delimiter='\t')
plat,plon=points['lat'][:],points['long'][:]

# read xml templates

tree0=et.ElementTree(file='loc0.xml')
root=tree0.getroot() #rss
channel=root[0] ##channel

channel[2][0].text=et.CDATA(channel[2][0].text)

for elem in list(tree0.getiterator()):
   if elem.tag == 'pubDate': elem.text='01 Jan 2016 00:00'

tree=et.ElementTree(file='loc.xml')
litem=tree.getroot()[0]

#print et.tostring(litem,pretty_print=True)


def modtxt(tag,val):
   for elem in nitem.findall(tag): elem.text=val

# APPEND LOCATIONS


for l in range(points.shape[0]):
   if (lat.min() < points['lat'][l] < lat.max()) & (lon.min() < points['long'][l] < lon.max()):
  #    print points.values[l]
       nitem=copy.deepcopy(litem)

       i=np.abs(lon-plon[l]).argmin()
       j=np.abs(lat-plat[l]).argmin()

 #     print l, i,j

       bd=b[i-1:i+2,j-1:j+2]

       bda=np.where(bd>20.)

       if np.size(bda) > 0 : 

          bloc=bd[bda].min()

          pi,pj=np.where(bd==bd[bda].min())

          ip,jp=i-1+pi.ravel()[0],j-1+pj.ravel()[0]

          tseries=hpr[:,ip,jp]

          maxh=tseries.max()
          print maxh

#         if maxh-tseries[0] < 0.05 : continue

          loc=np.argwhere(tseries==maxh)
          print loc

          dh=tseries-tseries[0]
          try: 
             arr=np.argwhere(dh>0.05)[0]
          except:
             arr=loc


          tMax=tpr[loc.ravel()[0]]/3600.
          tArr=tpr[arr.ravel()[0]]/3600.

#      print points.values[l]

          for elem in nitem.findall('title'): 
            elem.text=et.CDATA('{}: {} ({} m)'.format(points['$country'][l],points['$name'][l],maxh))

          for elem in nitem.findall('description'): 
           elem.text=et.CDATA('Country: {}\n \
Location: {}\n \
Time (hh:mm): {}\n \
Maximum Height: {} m\n \
Time Max (hh:mm): {}'.format(points['$country'][l],points['$name'][l],'{}:00'.format(tArr.astype(int)),maxh,'{}:00'.format(tMax.astype(int))))

          modtxt('pubDate',datetime.datetime.strftime(tend,'%d %b %Y %H:00:00'))

          modtxt('cityName',points['$name'][l])

          modtxt('country',points['$country'][l])

          modtxt('maxHeight',maxh.astype(str))

          modtxt('ID',points['id'][l].astype(str))

          modtxt('timeMaxH','{}:00'.format(tMax.astype(int)))

          modtxt('timeMaxH_value',tMax.astype(str))

          modtxt('timeArrival','{}:00'.format(tArr.astype(int)))

          modtxt('timeArrival_value',tArr.astype(str))

          modtxt('cityClass',points['cityclass'][l].astype(str))

          if not pandas.isnull(points['popest'][l]) :
              modtxt('popEst',points['popest'][l])
          else:
              modtxt('popEst','')

          modtxt('{http://www.w3.org/2003/01/geo/wgs84_pos#}long',lon[ip].astype(str))
          modtxt('{http://www.w3.org/2003/01/geo/wgs84_pos#}lat',lat[jp].astype(str))

          modtxt('depth',bloc.astype(str))

          channel.append(nitem)


tree0.write(SAVEPATH+'{}/locations.xml'.format(fname),xml_declaration=True,encoding='utf-8',method='xml')


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


