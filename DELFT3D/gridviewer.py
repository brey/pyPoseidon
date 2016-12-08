import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
from grid import *
from dep import *
#from nefis2 import *
import mdf
import sys
from netCDF4 import Dataset


def view(path,basename):

  inp, inp_order = mdf.read(path+basename+'.mdf')

  grdf=inp['Filcco'] # filename of grid file 
  depf=inp['Fildep'] # filename of dep file

  
  grid = Grid.fromfile(path+grdf)
  lon=grid.x
  lat=grid.y

  #lon=d.getdata('map-const','XZ')
  #lat=d.getdata('map-const','YZ')

  llcrnrlon=lon.min()
  llcrnrlat=lat.min()
  urcrnrlon=lon.max()
  urcrnrlat=lat.max()

  lon_0=(lon.max()-lon.min())/2.

  fig1 = plt.figure(figsize=(10,6))
# ax = fig1.add_axes([0.1,0.1,0.8,0.8])
  ax=fig1.add_subplot(111)

  m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,
            projection='cyl',lat_1=llcrnrlat,lat_2=urcrnrlat,lon_0=lon_0,
            resolution ='i',area_thresh=1000.)

  dep  = Dep.read(path+depf,grid.shape)

  x, y = m(lon, lat)

  bt=dep.val[:-1,:-1]

  #BATH = m.contour(x,y,bt,10,linewidths=0.5,colors='k',animated=True)
  #plt.clabel(BATH, fmt = '%03g', colors = 'k', fontsize=14)
  BATH = m.contourf(x,y,bt,20,cmap=plt.cm.RdBu_r,animated=True)
  m.colorbar()
# GP=m.scatter(x,y,s=2)

  parallels = np.arange(-90.,90,5.)
  meridians = np.arange(0.,360.,5.)


  m.drawcoastlines(linewidth=1.5)
  m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
  m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

  fig1.tight_layout(pad=1.08, h_pad=None, w_pad=None,rect=(0.05, 0, 1, 1))
 #fig1.tight_layout()

#  ax.set_title('bathymetry  ')

  try:
     __IPYTHON__ 
     plt.show(block=False)
  except:
     plt.show()



if __name__ == "__main__":
    path=sys.argv[1]
    basename=sys.argv[2]
    view(path,basename)

