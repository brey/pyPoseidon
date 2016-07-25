import numpy as np
import matplotlib.pyplot as plt
import datetime

from nefis2 import *

from grid import *

from mpl_toolkits.basemap import Basemap, shiftgrid


ni=200
nj=100
lon0=-10.
lat0=28.
lon1=48.
lat1=48.

date=datetime.date.today()
time=0 

d = Nefis('../test/trim-med.dat','../test/trim-med.def',access='r')

h = d.getdata('map-series','S1')
u = d.getdata('map-series','WINDU')
v = d.getdata('map-series','WINDV')
#lon = d.getdata('map-series','XZ')
#lat = d.getdata('map-series','YZ')

h = np.array(h)
u = np.array(u)
v = np.array(v)



m = Basemap(projection='cyl',llcrnrlat=lat0,urcrnrlat=lat1,\
             llcrnrlon=lon0,urcrnrlon=lon1,resolution='l')

# define parallels and meridians to draw.
parallels = np.arange(-90.,90,20.)
meridians = np.arange(0.,360.,20.)


grd = Grid.fromfile('../test/med.grd')


uu=u[0,:,:].T
vv=v[0,:,:].T
hh=h[0,:,:].T

vel=np.sqrt(uu**2+vv**2)

fig1 = plt.figure(figsize=(8,10))
ax = fig1.add_axes([0.1,0.1,0.8,0.8])


CS1 = m.contour(grd.x,grd.y,vel,10,linewidths=0.5,colors='k',animated=True)
CS2 = m.contourf(grd.x,grd.y,vel,10,cmap=plt.cm.RdBu_r,animated=True)

Q = m.quiver(grd.x,grd.y,uu,vv,scale=1500)
# make quiver key.
qk = plt.quiverkey(Q, 0.1, 0.1, 20, '20 m/s', labelpos='W')
# draw coastlines, parallels, meridians.
m.drawcoastlines(linewidth=1.5)
m.drawparallels(parallels)
m.drawmeridians(meridians)

cb = m.colorbar(CS2,"right", size="5%", pad="2%")
cb.set_label('m/s')
# set plot title
ax.set_title('Wind Speed at '+ datetime.datetime.strftime(date+datetime.timedelta(hours=time), '%a %b %d  %H:%M:%S %Z %Y'))


fig2 = plt.figure(figsize=(8,10))
ax = fig2.add_axes([0.1,0.1,0.8,0.8])

H = m.contourf(grd.x,grd.y,hh,10,cmap=plt.cm.RdBu_r,animated=True)
# draw coastlines, parallels, meridians.
m.drawcoastlines(linewidth=1.5)
m.drawparallels(parallels)
m.drawmeridians(meridians)
cb = m.colorbar(CS2,"right", size="5%", pad="2%")
cb.set_label('m')
# set plot title
ax.set_title('Water Level at '+ datetime.datetime.strftime(date+datetime.timedelta(hours=time), '%a %b %d  %H:%M:%S %Z %Y'))
    

plt.show(block=False)
