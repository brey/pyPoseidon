import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from mpl_toolkits.basemap import Basemap, shiftgrid
import scipy.interpolate
import sys
import grid

#PATH='../BATHYMETRY/'

def pltm(minlat,maxlat,minlon,maxlon,lons,lats,topo,title=None):

# Create map
 m = Basemap(projection='cyl', llcrnrlat=minlat,urcrnrlat=maxlat,llcrnrlon=minlon, urcrnrlon=maxlon,resolution='l')
 fig = plt.figure(figsize=(10,8))
 cs = m.contourf(lons,lats,topo,cmap=plt.cm.jet)
 m.drawcoastlines()
 m.drawmapboundary()
 plt.title(title)
 cbar = plt.colorbar(orientation='horizontal', extend='both')
 cbar.ax.set_xlabel('meters')
 
# Save figure (without 'white' borders)
#plt.savefig('topo.png', bbox_inches='tight')
 plt.show(block=False)

def readem(minlat,maxlat,minlon,maxlon,filename,grid_x=None,grid_y=None,plot=False,interpolate=False):
#file=PATH+'GLOBAL/GEBCO_2014_2D.nc'
# open NetCDF data in 
 nc = netCDF4.Dataset(filename)
 ncv = nc.variables
#print ncv.keys()
 nv=ncv.keys()
 if nv[0] == 'elevation' :
          n3,n2,n1=nv[0],nv[1],nv[2]
 else:
          n1,n2,n3=nv[0],nv[1],nv[2]

 lon = ncv[n1][:]
 lat = ncv[n2][:]

 if 'topo30.grd' not in filename:

  if maxlon > 180:

    lon=lon+180.

    i1=np.abs(lon-minlon).argmin()
    if lon[i1] > minlon: i1=i1-1
    i2=np.abs(lon-maxlon).argmin()
    if lon[i2] < maxlon: i2=i2+1

    j1=np.abs(lat-minlat).argmin()
    if lat[j1] > minlat: j1=j1-1
    j2=np.abs(lat-maxlat).argmin()
    if lat[j2] < maxlat: j2=j2+1

    lons, lats = np.meshgrid(lon[i1:i2],lat[j1:j2])

    zlon=lon.shape[0]

    topo = ncv[n3][j1:j2,zlon/2+i1:]
    topo = np.hstack([topo,ncv[n3][j1:j2,:i2-zlon/2]])

  else:

    i1=np.abs(lon-minlon).argmin()
    if lon[i1] > minlon: i1=i1-1
    i2=np.abs(lon-maxlon).argmin()
    if lon[i2] < maxlon: i2=i2+1

    j1=np.abs(lat-minlat).argmin()
    if lat[j1] > minlat: j1=j1-1
    j2=np.abs(lat-maxlat).argmin()
    if lat[j2] < maxlat: j2=j2+1

    lons, lats = np.meshgrid(lon[i1:i2],lat[j1:j2])
    topo = ncv[n3][j1:j2,i1:i2]

 else:

  if minlon < 0:

    lon=lon-180.

    i1=np.abs(lon-minlon).argmin()
    if lon[i1] > minlon: i1=i1-1
    i2=np.abs(lon-maxlon).argmin()
    if lon[i2] < maxlon: i2=i2+1

    j1=np.abs(lat-minlat).argmin()
    if lat[j1] > minlat: j1=j1-1
    j2=np.abs(lat-maxlat).argmin()
    if lat[j2] < maxlat: j2=j2+1

    lons, lats = np.meshgrid(lon[i1:i2],lat[j1:j2])

    zlon=lon.shape[0]

    topo = ncv[n3][j1:j2,zlon/2+i1:]
    topo = np.hstack([topo,ncv[n3][j1:j2,:i2-zlon/2]])

  else:

    i1=np.abs(lon-minlon).argmin()
    if lon[i1] > minlon: i1=i1-1
    i2=np.abs(lon-maxlon).argmin()
    if lon[i2] < maxlon: i2=i2+1

    j1=np.abs(lat-minlat).argmin()
    if lat[j1] > minlat: j1=j1-1
    j2=np.abs(lat-maxlat).argmin()
    if lat[j2] < maxlat: j2=j2+1

    lons, lats = np.meshgrid(lon[i1:i2],lat[j1:j2])
    topo = ncv[n3][j1:j2,i1:i2]

 title=filename.split('/')[-1]
 if plot : pltm(minlat,maxlat,minlon,maxlon,lons,lats,topo,title=title +'- READ')

 if interpolate :
# interpolate on the given grid
  #flip on lat to make it increasing for RectBivariateSpline
  ilon=lons[0,:]
  ilat=lats[:,0]
  sol=scipy.interpolate.RectBivariateSpline(ilon,ilat,topo.T)

  itopo=[]
  for x,y in zip(grid_x.ravel(),grid_y.ravel()):
      itopo.append(sol(x,y).ravel()[0])

  itopo=np.array(itopo)
  itopo=itopo.reshape(grid_x.shape)
  if plot : pltm(minlat,maxlat,minlon,maxlon,grid_x,grid_y,itopo,title=title+'- interpolated')
  return itopo
 else:
  return lons,lats,topo



#############################################################################################################
### MAIN


# Definine the domain of interest
if __name__ == "__main__":
 try:
    minlon=sys.argv[1]
    maxlon=sys.argv[2]
    minlat=sys.argv[3]
    maxlat=sys.argv[4]
    filename=sys.argv[5]
 except:
    minlat = 28.5
    maxlat = 47.5
    minlon = -5.5
    maxlon = 43.

 l1,l2,ba = readem(float(minlat),float(maxlat),float(minlon),float(maxlon),filename)
 m1,m2,bg = readem(float(minlat),float(maxlat),float(minlon),float(maxlon),'../BATHYMETRY/GLOBAL/topo30.grd')

 topo=ba-bg
 pltm(float(minlat),float(maxlat),float(minlon),float(maxlon),m1,m2,topo,title='GEBCO_Deltares-GEBCO (DIFF)')
 plt.show(block=False)
