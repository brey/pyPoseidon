
# coding: utf-8

# In[1]:

from matplotlib import animation,rc
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from grid import *
from dep import *
import glob
import numpy as np
import sys

from IPython.display import HTML


# In[2]:

rc('animation',html='html5')


# In[3]:

PATH='/home/critechuser/TC/THOMA/20101031.12/'
path0='/home/critechuser/TC/THOMA/'


# In[4]:

basename='TOMAS'


# In[5]:

grid=Grid.fromfile(PATH+basename+'.grd')
lon=grid.x[0,:].data
lat=grid.y[:,0].data


# In[6]:

deb=Dep.read(PATH+basename+'.dep',grid.shape)
b=deb.val
w=np.isnan(b)


# ### Merging the 12h forecasts for longer time span

# In[7]:

tstart=datetime.datetime(2010,10,21,0)


# In[8]:

tend=datetime.datetime(2010,10,31,12)


# In[9]:

dt=(tend-tstart).total_seconds()
ndt=dt/(3600*12)
ndt=np.int(ndt)+1

nl=np.arange(ndt)
n=2
im=0

spl=[nl[i:i + n] for i in xrange(0, len(nl), n)]

for k in spl:

 im += 1
# In[9]:

# In[10]:

 hcombined=[]
 tw=[]


# In[11]:

 for it in k:
    idate=tstart+datetime.timedelta(hours=12*it)
    folder=datetime.datetime.strftime(idate,'%Y%m%d.%H')
    path=path0+'{}/'.format(folder)
    d=Dataset(path+'trim-'+basename+'.nc')
    hcombined.append(d['S1'][:12,:,:])
    time=d['time'][:12]
    tm =(idate-tstart).total_seconds()+(time-time[0])
    for l in tm : tw.append(tstart+datetime.timedelta(0,int(l)))
    


# In[12]:

 tcw=np.array(tw).flatten()


# In[13]:

 hh=np.array(hcombined)


# In[14]:

 hh=hh.reshape(tcw.size,w.shape[1],w.shape[0])


# In[15]:

 ha=np.transpose(hh,axes=(0,2,1))


# In[16]:

 iframes=ha.shape[0]
 fig = plt.figure(figsize=(10,8))
 ax = fig.add_subplot(111)
 bh=np.ma.masked_where(w==True,ha[0,:,:])
 #H=plt.imshow(np.flipud(h[0,:,:].T),animated=True, vmin=h.min(), vmax=h.max())
 H=ax.imshow(np.flipud(bh),animated=True, vmin=-.5, vmax=.5)
 plt.colorbar(H, orientation='horizontal')
 date_text=ax.text(0.02,0.95, datetime.datetime.strftime(tcw[0],'%a %b %d  %H:%M:%S %Z %Y' ))


# In[17]:

 def animate(i):
    bh=np.ma.masked_where(w==True,ha[i,:,:])
#    H.set_array(np.flipud(ha[i,:,:]))
    H.set_array(np.flipud(bh))
    date_text.set_text(datetime.datetime.strftime(tcw[i],'%a %b %d  %H:%M:%S %Z %Y' ))
    return H


# In[18]:

# call the animator.  blit=True means only re-draw the parts that have changed.
 anim = animation.FuncAnimation(fig, animate, frames=iframes, interval=200, blit=False, repeat=False)


# In[19]:

#  display the animation
#anim


# In[ ]:

#save animation
 anim.save('tmp/THOMAS/thomas{}.mp4'.format(im), fps=10, extra_args=['-vcodec','libx264','-pix_fmt', 'yuv420p'])
 print 'movie {} saved'.format(im)
 
 anim=[]
 hh=[]
 ha=[]
 bh=[]
 H=[]


