#python modules
import numpy as np
import glob
import pickle
from netCDF4 import Dataset
import folium
import matplotlib.pyplot as plt
import pandas as pd
import datetime
import base64
from folium.element import IFrame
from folium import plugins
import cmocean
from mpl_toolkits.basemap import Basemap
from folium import colormap as cm
from folium.map import *


# In[3]:

get_ipython().magic('load_ext pypath_magic')


# In[4]:

# add temporarily to the python the required local paths. Make sure you clean up afterwards.
# It won't run if ~/pypath_magic.pth is not empty to begin with.
get_ipython().magic('pypath -a /home/critechuser/REPOS/StormS/src/')
get_ipython().magic('pypath -a /home/critechuser/REPOS/StormS/utils/')


# In[5]:

#local modules
from grid import *
from dep import *
#from point_history import getmes
#from get_point_map import get
#from get_point_history import pget


# In[6]:

# clean up 
get_ipython().magic('pypath -d /home/critechuser/REPOS/StormS/src/')
get_ipython().magic('pypath -d /home/critechuser/REPOS/StormS/utils/')


# In[32]:

import itertools


# In[7]:

get_ipython().magic('matplotlib notebook')


# In[8]:

path='tmp/'


# In[9]:

# check contents of folder
tfolder = glob.glob(path+'*')
tfolder


# In[18]:

# exclude pkl file from contents of folder
tfolder = [t for t in tfolder if '.pkl' not in t]
tfolder
tfolder = [t for t in tfolder if 'png' not in t]
tfolder
tfolder = [t for t in tfolder if 'html' not in t]
tfolder


# In[19]:

tl = [t.split('/')[1] for t in tfolder]


# In[20]:

tl


# In[21]:

calc_dir=path+tl[0]+'/' # define project folder, usually the first time stamp folder


# In[22]:

pfile = glob.glob(path+'*.pkl')[0]
pfile


# In[23]:

# read pref file
with open(pfile, 'r') as f:
    dic=pickle.load(f)


# In[24]:

basename=dic['bname']


# In[25]:

#read obs points dictionary
with open(calc_dir+basename+'.pkl', 'r') as f:
    ptr=pickle.load(f)


# ### read grid 

# In[26]:

d=Dataset(calc_dir+'trim-'+basename+'.nc')


# In[27]:

xg=d['XCOR'][:] #grid points
yg=d['YCOR'][:]


# In[28]:

xz=d['XZ'][:] #
yz=d['YZ'][:]


# ## Define map coordinates

# In[29]:

#center of lat/lon window
plat = yz.mean()
plon = xz.mean()


# In[109]:

get_ipython().run_cell_magic('skip', '', 'a=np.arange(xg.shape[0])\nb=np.arange(xg.shape[1])\nglist=[]\nfor r in itertools.product(a,b): glist.append([r[0],r[1]])\n#glist')


# In[364]:

grid=FeatureGroup(name='Grid')


# In[365]:

#draw horizontal lines of grid
for k in xrange(xg.shape[0]):
    try:
        xgrid=zip(yg[k,:],xg[k,:])
        ygrid=zip(yg[k,:],xg[k,:])
        folium.PolyLine(xgrid,weight=.5, color='black').add_to(grid)
        folium.PolyLine(ygrid,weight=.5, color='black').add_to(grid)
    except Exception as e:
        print e
        pass


# In[366]:

#draw horizontal lines of grid
for k in xrange(1,xz.shape[0]-1):
    try:
        xgrid=zip(yz[k,1:-1],xz[k,1:-1])
        ygrid=zip(yz[k,1:-1],xz[k,1:-1])
        folium.PolyLine(xgrid,weight=.3, color='blue').add_to(grid)
        folium.PolyLine(ygrid,weight=.3, color='blue').add_to(grid)
    except Exception as e:
        print e
        pass


# In[367]:

#draw vertical
for k in xrange(xg.shape[1]):
    try:
        xgrid=zip(yg[:,k],xg[:,k])
        ygrid=zip(yg[:,k],xg[:,k])
        folium.PolyLine(xgrid,weight=.5, color='black').add_to(grid)
        folium.PolyLine(ygrid,weight=.5, color='black').add_to(grid)
    except Exception as e:
        print e
        pass


# In[368]:

#draw vertical
for k in xrange(1,xz.shape[1]-1):
    try:
        xgrid=zip(yz[1:-1,k],xz[1:-1,k])
        ygrid=zip(yz[1:-1,k],xz[1:-1,k])
        folium.PolyLine(xgrid,weight=.3, color='blue').add_to(grid)
        folium.PolyLine(ygrid,weight=.3, color='blue').add_to(grid)
    except Exception as e:
        print e
        pass


# In[114]:

get_ipython().run_cell_magic('skip', '', "for k,l in glist:\n# draw grid points        \n  folium.CircleMarker([yg[k,l],xg[k,l]], popup='{},{}'.format(k,l), radius=50).add_to(mapa)\n# draw staggered pressure points\n  folium.CircleMarker([yz[k,l],xz[k,l]], popup='{},{}'.format(k,l), radius=50, fill_color='red').add_to(mapa)")


# In[380]:

mapa = folium.Map(location=[plat, plon], zoom_start=4)


# In[381]:

#%%skip
folium.LatLngPopup().add_to(mapa) # click to show lat lon


# In[382]:

mapa.add_children(grid)
folium.LayerControl().add_to(mapa)


# In[383]:

mapa


# ## obs points

# In[384]:

df = pd.read_csv('../src/SeaLevelBuoys2.csv', encoding = 'utf8')


# In[385]:

df.head()


# In[386]:

df['jref'] = [ptr[i] if i in ptr.keys() else np.nan for i in df.ID.values]


# In[387]:

df.head()


# In[388]:

#drop the not relevant locations
df = df.dropna(subset=['jref'])
df.jref = df.jref.apply(pd.to_numeric) #convert to integer


# In[389]:

df = df.reset_index(drop=True) #reset the index


# In[390]:

mapa = folium.Map(location=[plat, plon], zoom_start=4) # reset map


# In[34]:

## put on the map


# In[35]:

mc = folium.MarkerCluster().add_to(mapa)


# In[36]:

for idx, v in df.iterrows():
    folium.Marker([v.lat,v.lon], popup=v['name']).add_to(mc)


# In[37]:

mapa


# ## overlay obs point data

# In[43]:

#read netcdf
h=Dataset(calc_dir+'trih-'+basename+'.nc')


# In[44]:

ha = h['ZWL'][:]
t = h['time'][:]


# In[45]:

idate=datetime.datetime.strptime(tl[0],'%Y%m%d.%H' )


# In[46]:

#construct datetime
tw=[]
for it in t:
        tw.append(idate+datetime.timedelta(seconds=np.int(it)))
ttw=[(item-t[0])/60. for item in t]


# In[47]:

resolution, width, height = 75, 8, 4


# In[146]:

mapa = folium.Map(location=[plat, plon], zoom_start=4) # reset map


# In[147]:

marker_cluster = folium.MarkerCluster().add_to(mapa)


# In[148]:

plt.ioff()
for iobs in range(5):#  for testing see full output from scipt
  
    gdic={'time':tw,'ha':ha[:,df.iloc[iobs].jref.astype(int)]}

    gf = pd.DataFrame(gdic)

    gf.time = gf.time.apply(pd.to_datetime) #convert to datetime
    gf = gf.set_index('time')

#create a graph, choose a station
    station = df.iloc[iobs]['name']
    print station
    
    fig, ax = plt.subplots(figsize=(width, height))
    ax = gf.plot(ax = ax, legend=False)
    ax.set_ylabel('Sea surface height (m)')
    ax.set_title(station)
    png = 'tmp/png/{}.png'.format(station.encode('utf-8').strip())
    fig.savefig(png, dpi=resolution)
    
    encoded = base64.b64encode(open(png, 'rb').read())
    html = '<img src="data:image/png;base64,{}">'.format
    iframe = IFrame(html(encoded), width=(width*resolution)+20, height=(height*resolution)+20)
    popup = folium.Popup(iframe, max_width=2650)

    icon = folium.Icon(color="red", icon="ok")
    marker = folium.Marker(location=[df.iloc[iobs].lat, df.iloc[iobs].lon], popup=popup, icon=icon)
    mapa.add_children(marker);


# In[149]:

mapa


# ## save html

# In[144]:

mapa.save('tmp/viewer.html')


# ## add storm surge to map 

# In[30]:

#read values
hz=d['S1'][:]


# In[31]:

plt.ion()


# In[154]:

#READ GRID /BATHYMETRY
grd=Grid.fromfile(calc_dir+basename+'.grd')
deb=Dep.read(calc_dir+basename+'.dep',grd.shape)
b=deb.val[:,:]
w=np.isnan(b)
w=w.T


# In[92]:

lons=xz[1:-1,1:-1]
lats=yz[1:-1,1:-1]


# In[93]:

lon0, lon1, lat0, lat1 = lons.min(),lons.max(),lats.min(),lats.max()


# In[94]:

#define projection Mercator
m = Basemap(projection='merc',llcrnrlat=lat0,urcrnrlat=lat1,             llcrnrlon=lon0,urcrnrlon=lon1,resolution='h')


# In[95]:

# define parallels and meridians to draw.
parallels = np.arange(-90.,90,20.)
meridians = np.arange(0.,360.,20.)


# In[96]:

xx,yy=m(lons,lats)


# In[155]:

#mask values
w.shape,hz[-1,1:-1,1:-1].shape


# In[305]:

pz =np.ma.masked_where(w[1:-1,1:-1]==True, hz[-1,1:-1,1:-1])


# In[98]:

xx.shape,yy.shape,pz.shape


# In[356]:

# create png image
plt.figure()
m.contourf(xx,yy,pz,cmap=plt.get_cmap('YlGn'))
# draw coastlines, parallels, meridians.
#m.drawcoastlines(linewidth=1.5)
#m.drawparallels(parallels)
#m.drawmeridians(meridians)


plt.savefig('tmp/sspng',transparent=True,bbox_inches='tight', pad_inches=0)


# In[192]:

from scipy.ndimage import imread


# In[352]:

# read in png file to numpy array
data = imread('tmp/sspng.png')


# ## WITH numpy manipulation 

# In[333]:

from PIL import Image, ImageEnhance


# In[334]:

im = Image.fromarray(plt.get_cmap('YlGn')(np.flipud(pz.T), bytes=True))


# In[347]:

plt.imshow(im)


# In[341]:

im2=ImageEnhance.Sharpness(im)


# In[348]:

plt.imshow(im2.enhance(2))


# In[316]:

v = np.asarray(im)


# In[331]:

v.shape


# ## test image overlay

# In[321]:

min_lon=json.loads(str(lon0))
max_lon=json.loads(str(lon1))
min_lat=json.loads(str(lat0))
max_lat=json.loads(str(lat1))


# In[360]:

mapa = folium.Map(location=[plat, plon], zoom_start=4) # reset map


# In[361]:

Surge = plugins.ImageOverlay(data, opacity=0.8,         bounds =[[min_lat, min_lon], [max_lat, max_lon]],        attr = 'Storm Surge')#,mercator_project=True)

# add a name to the layer
Surge.layer_name = 'Storm Surge'

# Overlay the image
mapa.add_children(Surge)

# add the layer control
folium.LayerControl().add_to(mapa)


# In[362]:

colormap = cm.linear.YlGn.scale(0, 3)#.to_step(10)
colormap.caption = 'Storm Surge'
mapa.add_child(colormap)


# ## test heatmap plugin

# In[258]:

mapa = folium.Map(location=[plat, plon], zoom_start=4) # reset map


# In[259]:

ss = [[a,b,c] for (a,b,c) in zip(lats.flatten(),lons.flatten(),hz[-1,1:-1,1:-1].flatten())]


# In[260]:

from folium.plugins import HeatMap

HeatMap(ss).add_to(mapa)

mapa


# In[240]:

mapa.save('tmp/viewer.html')

