import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
from grid import *
from dep import *
#from nefis2 import *
import mdf
import sys
from netCDF4 import Dataset


def view(path1,path2,basename):

  inp, inp_order = mdf.read(path1+basename+'.mdf')

  grdf=inp['Filcco'] # filename of grid file 
  depf=inp['Fildep'] # filename of dep file

  
  grid1 = Grid.fromfile(path1+grdf)
  lon1=grid1.x
  lat1=grid1.y

  grid2 = Grid.fromfile(path2+grdf)
  lon2=grid2.x
  lat2=grid2.y

  #lon=d.getdata('map-const','XZ')
  #lat=d.getdata('map-const','YZ')

  llcrnrlon=min(lon1.min(),lon2.min())
  llcrnrlat=min(lat1.min(),lat2.min())
  urcrnrlon=max(lon1.max(),lon2.max())
  urcrnrlat=max(lat1.max(),lat2.max())

  lon_0=(lon1.max()-lon1.min())/2.

  fig1 = plt.figure(figsize=(10,8))
  ax = fig1.add_axes([0.1,0.1,0.8,0.8])

  m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,
            projection='cyl',lat_1=llcrnrlat,lat_2=urcrnrlat,lon_0=lon_0,
            resolution ='i',area_thresh=1000.)

  x1, y1 = m(lon1, lat1)
  x2, y2 = m(lon2, lat2)

  m.scatter(x1,y1,s=20)
  m.scatter(x2,y2,s=40,marker='+')

  parallels = np.arange(-90.,90,5.)
  meridians = np.arange(0.,360.,5.)


  m.drawcoastlines(linewidth=1.5)
  m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
  m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

  ax.set_title('grid')

  try:
     __IPYTHON__ 
     plt.show(block=False)
  except:
     plt.show()



if __name__ == "__main__":
    path1=sys.argv[1]
    path2=sys.argv[2]
    basename=sys.argv[3]
    view(path1,path2,basename)

