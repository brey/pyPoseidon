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


#local modules
from grid import *
from dep import *
from cmap2gdf import *
from float_image import FloatImage
from Bind import BindColormap
#from point_history import getmes
#from get_point_map import get
#from get_point_history import pget


path='../../../test/'

# check contents of folder
tfolder = glob.glob(path+'*')

# exclude pkl file from contents of folder
tfolder = [t for t in tfolder if '.pkl' not in t]
tfolder = [t for t in tfolder if 'png' not in t]
tfolder = [t for t in tfolder if 'html' not in t]

tl = [t.split('/')[-1] for t in tfolder]

calc_dir=path+tl[0]+'/' # define project folder, usually the first time stamp folder

pfile = glob.glob(path+'*.pkl')[0]

# read pref file
with open(pfile, 'r') as f:
    dic=pickle.load(f)

basename=dic['bname']

#read obs points dictionary
with open(calc_dir+basename+'.pkl', 'r') as f:
    ptr=pickle.load(f)

# ### read grid 
d=Dataset(calc_dir+'trim-'+basename+'.nc')

xg=d['XCOR'][:] #grid points
yg=d['YCOR'][:]

xz=d['XZ'][:] #
yz=d['YZ'][:]

# ## Define map coordinates

#center of lat/lon window
plat = yz.mean()
plon = xz.mean()

grid=FeatureGroup(name='Grid')

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

mapa = folium.Map(location=[plat, plon], zoom_start=4)

folium.LatLngPopup().add_to(mapa) # click to show lat lon

mapa.add_children(grid)


# ## obs points

# In[384]:

df = pd.read_csv('../src/SeaLevelBuoys2.csv', encoding = 'utf8')

df['jref'] = [ptr[i] if i in ptr.keys() else np.nan for i in df.ID.values]

#drop the not relevant locations
df = df.dropna(subset=['jref'])
df.jref = df.jref.apply(pd.to_numeric) #convert to integer


df = df.reset_index(drop=True) #reset the index

obs=FeatureGroup(name='Observation points')
## put on the map

mc = folium.MarkerCluster().add_to(obs)


for idx, v in df.iterrows():
    folium.Marker([v.lat,v.lon], popup=v['name']).add_to(mc)

mapa.add_children(obs)

# ## overlay obs point data

comp=FeatureGroup(name='Observations')
#read netcdf
h=Dataset(calc_dir+'trih-'+basename+'.nc')

ha = h['ZWL'][:]
t = h['time'][:]

idate=datetime.datetime.strptime(tl[0],'%Y%m%d.%H' )

#construct datetime
tw=[]
for it in t:
        tw.append(idate+datetime.timedelta(seconds=np.int(it)))
ttw=[(item-t[0])/60. for item in t]

resolution, width, height = 75, 8, 4

marker_cluster = folium.MarkerCluster().add_to(comp)

plt.ioff()
for iobs in range(df.shape[0]):#  for testing see full output from scipt
  
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
    png = path+'png/{}.png'.format(station.encode('utf-8').strip())
    fig.savefig(png, dpi=resolution)
    
    encoded = base64.b64encode(open(png, 'rb').read())
    html = '<img src="data:image/png;base64,{}">'.format
    iframe = IFrame(html(encoded), width=(width*resolution)+20, height=(height*resolution)+20)
    popup = folium.Popup(iframe, max_width=2650)

    icon = folium.Icon(color="red", icon="ok")
    marker = folium.Marker(location=[df.iloc[iobs].lat, df.iloc[iobs].lon], popup=popup, icon=icon)
    marker_cluster.add_children(marker);


mapa.add_children(comp)
# ## add storm surge to map 


#read values
hz=d['S1'][:]

#READ GRID /BATHYMETRY
grd=Grid.fromfile(calc_dir+basename+'.grd')
deb=Dep.read(calc_dir+basename+'.dep',grd.shape)
b=deb.val[:,:]
w=np.isnan(b)
w=w.T

lons=xz[1:-1,1:-1]
lats=yz[1:-1,1:-1]

lon0, lon1, lat0, lat1 = lons.min(),lons.max(),lats.min(),lats.max()

#mask values
w.shape,hz[-1,1:-1,1:-1].shape

pz =np.ma.masked_where(w[1:-1,1:-1]==True, hz[-1,1:-1,1:-1])

# contour plot

plt.figure()
nb_class = 20 
collec_poly = plt.contourf(lons,lats,pz, nb_class, cmap=plt.get_cmap('YlGn'), alpha=0.5)

gdf = collec_to_gdf(collec_poly) # From link above
gdf.to_json()
colors = [p.get_facecolor().tolist()[0] for p in collec_poly.collections]
gdf['RGBA'] = colors
gdf['RGBA'] = gdf['RGBA'].apply(convert_to_hex)

colors = []
Contours = folium.GeoJson(
    gdf,
    style_function=lambda feature: {
        'fillColor': feature['properties']['RGBA'],
        'color' : feature['properties']['RGBA'],
        'weight' : 1,
        'fillOpacity' : 0.5,
        }
    )

t = folium.FeatureGroup(name='Storm Surge')
t.add_children(Contours)

mapa.add_children(t)

# #### float an image

ec = 'https://ec.europa.eu/jrc/sites/jrcsh/themes/jrc_multisite_subtheme/logo.png'

FloatImage(ec, bottom=0, left=2).add_to(mapa)

ecmwf = 'http://www.ecmwf.int/sites/default/files/ECMWF_Master_Logo_RGB_nostrap.png'

FloatImage(ecmwf, bottom=95, left=5, height=4).add_to(mapa)

delft3d='https://oss.deltares.nl/image/image_gallery?uuid=13baea0c-1a8c-44c4-b4fa-89c7d4d68f3c&groupId=183920&t=1352473709104'

FloatImage(delft3d, bottom=85, left=5, height=8).add_to(mapa)

cm1 = cm.linear.YlGn.scale(0, 3)#.to_step(10)
cm1.caption = 'Storm Surge' 

mapa.add_child(cm1)
mapa.add_child(BindColormap(t,cm1))

mapa.add_child(folium.map.LayerControl())

mapa.save(path+'viewer.html')


