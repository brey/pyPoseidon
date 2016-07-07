import numpy as np
import matplotlib.pyplot as plt
import datetime

from netCDF4 import Dataset


from grid import *
from dep import *
import mdf

from mpl_toolkits.basemap import Basemap, shiftgrid

import matplotlib.animation as animation


####################################################################
# Specify ffmpeg path
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'

FFwriter = animation.FFMpegWriter()
####################################################################


#date=datetime.datetime(2016,2,15,12,0)


PATH='/mnt/web/brey/2016H/6/14/12/'
PATH='/mnt/web/brey/venice/2016/6/14/12/'
PATH='/mnt/web/brey/NEPARTAK/2016/7/6/00/'
tag='nepartak'

inp, ord = mdf.read(PATH+tag+'.mdf')

date=datetime.datetime.strptime(inp['Itdate'],'%Y-%m-%d')
date=date+datetime.timedelta(minutes=inp['Tstart'][0])

d = Dataset(PATH+'trim-'+tag+'.nc')

h = d.variables['S1'][:]
u = d.variables['WINDU'][:]
v = d.variables['WINDV'][:]
x = d.variables['XCOR'][:]
y = d.variables['YCOR'][:]


lon0=x.min()
lon1=x.max()
lat0=y.min()
lat1=y.max()

#READ GRID /BATHYMETRY
grd=Grid.fromfile(PATH+tag+'.grd')
deb=Dep.read(PATH+tag+'.dep',grd.shape)
b=deb.val[:-1,:-1]
w=b<0.
w=w.T



m = Basemap(projection='cyl',llcrnrlat=lat0,urcrnrlat=lat1,\
             llcrnrlon=lon0,urcrnrlon=lon1,resolution='l')

# define parallels and meridians to draw.
parallels = np.arange(-90.,90,20.)
meridians = np.arange(0.,360.,20.)


uu=u[-3,:,:]
vv=v[-3,:,:]
hh=h[-3,:,:]

nt=h.shape[0]

vel=np.sqrt(uu**2+vv**2)

###MASK
#um=np.ma.masked_array(uu,uu==-999.)
#vm=np.ma.masked_array(vv,vv==-999.)
vel=np.sqrt(uu**2+vv**2)


time=nt-1 

####### WATER VELOCITY ############

fig1 = plt.figure(figsize=(10,8))
ax = fig1.add_axes([0.1,0.1,0.8,0.8])

clevs=10


CS1 = m.contour(x,y,vel,clevs,linewidths=0.5,colors='k',animated=True)
CS2 = m.contourf(x,y,vel,clevs,cmap=plt.cm.RdBu_r,animated=True)

ni,nj=uu.shape
stepi=ni/50
stepj=nj/30

Q = m.quiver(x[::stepi,::stepj],y[::stepi,::stepj],uu[::stepi,::stepj],vv[::stepi,::stepj], units='x')#,scale=500)
# make quiver key.
qk = plt.quiverkey(Q, 0.1, 0.1, 20, '20 m/s', labelpos='W')
# draw coastlines, parallels, meridians.
m.drawcoastlines(linewidth=1.5)
m.drawparallels(parallels)
m.drawmeridians(meridians)

cb = m.colorbar(CS2,"right", size="5%", pad="2%")
cb.set_label('m/s')
# set plot title
ax.set_title('Winds at '+ datetime.datetime.strftime(date+datetime.timedelta(hours=time), '%a %b %d  %H:%M:%S %Z %Y'))


####### SEA LEVEL #########
bh=np.ma.masked_where(w==True,hh)

fig2 = plt.figure(figsize=(10,8))
ax = fig2.add_axes([0.1,0.1,0.8,0.8])

[[i,j]]=np.argwhere(bh==bh.max())
H = m.contourf(x,y,bh,clevs,cmap=plt.cm.RdBu_r,animated=True)
H1 =  m.plot(x[i,j],y[i,j],'kx',markersize=20)
fig2.text(.01,.05,'iter= '+np.str(nt-1),horizontalalignment='left',size='small')
# draw coastlines, parallels, meridians.
m.drawcoastlines(linewidth=1.5)
m.drawparallels(parallels)
m.drawmeridians(meridians)
cb = m.colorbar(H,"right", size="5%", pad="2%")
cb.set_label('m')
# set plot title
ax.set_title('Water Level at '+ datetime.datetime.strftime(date+datetime.timedelta(hours=time), '%a %b %d  %H:%M:%S %Z %Y'))
    

plt.show(block=False)



# check convergence

#for i in range(h.shape[0]-1):
#   wl=np.abs(h[i+1,:,:]-h[i,:,:])
#   print wl.max()

iframes=h.shape[0]-1


def updatefig1(nt):
    global CS1,CS2
    bh=np.ma.masked_where(w==True,h[nt,:,:])
    for c in CS1.collections: c.remove()
    CS1 = m.contour(x,y,bh,clevs,linewidths=0.5,colors='k')
    for c in CS2.collections: c.remove()
    CS2 = m.contourf(x,y,bh,clevs,cmap=plt.cm.jet)

def updatefig2(nt):
    global H,H1,bh
    bh=np.ma.masked_where(w==True,h[nt,:,:])
    for c in H.collections: c.remove()
#   if fig2.get_default_bbox_extra_artists():
#      for c in fig2.get_default_bbox_extra_artists(): c.pop()
    fig2.delaxes(fig2.axes[1]) 
    fig2.texts.pop()
    H = m.contourf(x,y,bh,clevs,cmap=plt.cm.RdBu_r,animated=True)
#   [[i,j]]=np.argwhere(bh==bh.max())
#   H1 = m.plot(x[i,j],y[i,j],'kx',markersize=20)
    cb = m.colorbar(H,"right", size="5%", pad="2%")
    plt.title('Water Level at '+ datetime.datetime.strftime(date+datetime.timedelta(hours=nt), '%a %b %d  %H:%M:%S %Z %Y'))
    fig2.text(.01,.05,'iter= '+np.str(nt),horizontalalignment='left',size='small')

ani = animation.FuncAnimation(fig2, updatefig2, frames=iframes, repeat=False)

    
ani.save('sea_level.mp4',writer = FFwriter)


def anim():

 plt.ion(); plt.figure(3);

 for i in range(h.shape[0]):
    plt.clf()
    bh=np.ma.masked_where(w==True,h[i,:,:])
    H = m.contourf(x,y,bh,clevs,cmap=plt.cm.RdBu_r,animated=True)
# draw coastlines, parallels, meridians.
    m.drawcoastlines(linewidth=1.5)
    m.drawparallels(parallels)
    m.drawmeridians(meridians)
    cb = m.colorbar(H,"right", size="5%", pad="2%")
    cb.set_label('m')
# set plot title
    plt.title('Water Level at '+ datetime.datetime.strftime(date+datetime.timedelta(hours=i), '%a %b %d  %H:%M:%S %Z %Y'))
    plt.text(.01,.05,'iter= '+np.str(i),horizontalalignment='left',size='large')
    plt.draw();

#anim()
