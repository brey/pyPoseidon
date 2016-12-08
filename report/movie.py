
# coding: utf-8

# In[1]:

from matplotlib import animation,rc
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from grid import *
from dep import *
import glob
import numpy as np

from IPython.display import HTML


# In[2]:

rc('animation',html='html5')


# In[3]:

path0='/home/critechuser/TC/THOMA/'
PATH=path0+'20101021.00/'
basename='TOMAS'


# ### Merging the 12h forecasts for longer time span

# In[6]:

tstart=datetime.datetime(2010,10,21,0)


# In[7]:

tend=datetime.datetime(2010,10,31,12)


d=Dataset(PATH+'trim-'+basename+'.nc')
xz=d.variables['XZ'][:]
yz=d.variables['YZ'][:]


# In[8]:

dt=(tend-tstart).total_seconds()
ndt=dt/(3600*12)
ndt=np.int(ndt)

nl=np.arange(ndt)
n=2
im=0

spl=[nl[i:i + n] for i in xrange(0, len(nl), n)]

for k in spl:

 im += 1
# In[9]:

 hcombined=[]
 tw=[]

# In[10]:

 for it in k:
    idate=tstart+datetime.timedelta(hours=12*it)
    folder=datetime.datetime.strftime(idate,'%Y%m%d.%H')
    path=path0+'{}/'.format(folder)
    d=Dataset(path+'trim-'+basename+'.nc')
    hcombined.append(d['S1'][:12,:,:])
    time=d['time'][:12]
    tm =(idate-tstart).total_seconds()+(time-time[0])
    for l in tm : tw.append(tstart+datetime.timedelta(0,int(l)))
    


# In[11]:

 tcw=np.array(tw).flatten()


# In[12]:

 hh=np.array(hcombined)


# In[15]:

 hh=hh.reshape(tcw.size,xz.shape[1],xz.shape[0])


# In[16]:

 ha=np.transpose(hh,axes=(0,2,1))


# In[17]:

 iframes=ha.shape[0]
 fig = plt.figure()
 ax = fig.add_subplot(111)
 bh=np.ma.masked_where(xz==0,ha[0,:,:])
#H=plt.imshow(np.flipud(h[0,:,:].T),animated=True, vmin=h.min(), vmax=h.max())
 H=ax.imshow(bh,animated=True, vmin=-.5, vmax=2.)
#H=ax.contourf(xz,yz,bh,animated=True, vmin=-.5, vmax=2.)
# H=ax.contourf(xz,yz,bh,animated=True, vmin=-.5, vmax=2.)
 plt.colorbar(H, orientation='horizontal')
 date_text=ax.text(0.02,0.95, datetime.datetime.strftime(tcw[0],'%a %b %d  %H:%M:%S %Z %Y' ))


# In[18]:

 def animate(i):
    bh=np.ma.masked_where(xz==0,ha[i,:,:])
#    H.set_array(np.flipud(ha[i,:,:]))
    H.set_array(ha[i,:,:])
    date_text.set_text(datetime.datetime.strftime(tcw[i],'%a %b %d  %H:%M:%S %Z %Y' ))
    return H


# In[19]:


# call the animator.  blit=True means only re-draw the parts that have changed.
 anim = animation.FuncAnimation(fig, animate, frames=iframes, interval=200, blit=False, repeat=False)


# In[ ]:

#  display the animation
#anim


# In[20]:

#save animation
 anim.save('tmp/THOMAS/thomas{}.mp4'.format(im), fps=10, extra_args=['-vcodec','libx264'])


# In[ ]:


