import numpy as np
import pandas
import collections
import pickle
from itertools import product

from dep import *
from grid import *

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid



def newpoint(i,j,l1,l2,m,grd,bath):
       pool=[]
       loc=[]
       for ia,ja in product(range(-1,2),range(-1,2)):
                if  bath.val[i+ia,j+ja] > 0. : 
                            pi=i+ia
                            pj=j+ja
                            xx=grd.x.T[pi,pj] 
                            yy=grd.y.T[pi,pj]
                            m.plot(xx,yy,'gx',markersize=10,linewidth=6)
                            pool.append([pi,pj]) 
                            loc.append(np.sqrt((xx-l1)**2+(yy-l2)**2)) 
                else: 
                            pi=i 
                            pj=j

       if not loc : print l1,l2
       if loc:
           idx = np.array(loc).argmin()
           pi,pj = pool[idx]
           xx=grd.x.T[pi,pj] 
           yy=grd.y.T[pi,pj]
           m.scatter(xx,yy,s=60,facecolors ='none', edgecolors='b')
       return pi+1,pj+1 # the fortran/python index issue ??




def createf(path,basename,lat0,lat1,lon0,lon1,grd,bath):

 bath.val=bath.val.T

#  BUOYS WEBCRITECH
 dat=pandas.read_csv('SeaLevelBuoys2.csv')
 print dat.columns

 ID=dat[dat.columns[0]]
 lon=dat['lon']
 lat=dat['lat']
 names=dat['name']

 fig1 = plt.figure(figsize=(10,8))
 ax = fig1.add_axes([0.1,0.1,0.8,0.8])

 m = Basemap(projection='cyl',llcrnrlat=lat0,urcrnrlat=lat1,\
             llcrnrlon=lon0,urcrnrlon=lon1,resolution='l')

# define parallels and meridians to draw.
 parallels = np.arange(-90.,90,20.)
 meridians = np.arange(0.,360.,20.)

 iobs=[]
 loci=[]
 k=0
 pdic={}
 for l1,l2 in zip(lon,lat):
  if (lon0 <= l1 <= lon1) & (lat0 <= l2 <= lat1) :
    ih=np.abs(grd.x-l1).argmin()
    jh=np.abs((grd.y-l2).T).argmin()
# test
    xx=grd.x.T[ih,jh] 
    yy=grd.y.T[ih,jh]
    m.plot(xx,yy,'rx',markersize=10,linewidth=6)
    pi=ih+1 # the fortran/python index issue ??
    pj=jh+1 # the fortran/python index issue ??

    if bath.val[ih,jh] < 0. :  pi,pj=newpoint(ih,jh,l1,l2,m,grd,bath) # choose new wet point 
                 
    m1=np.argwhere(lon==l1) 
    m2=np.argwhere(lat==l2) 
    if m1.all()==m2.all():
     if (pi,pj) not in iobs:
      for m0 in m1.flatten():
        pdic[ID[m0]]=k
      iobs.append((pi,pj))
      loci.append((l1,l2))
      k=k+1
     else:
      g=np.where(np.array(iobs)==[pi,pj])
      c=collections.Counter(g[0])
      rk = [value for value, count in c.items() if count > 1][0]
      for m0 in m1.flatten():
        pdic[ID[m0]]=rk
    else: print m1,m2
    
 with open(path+basename+'.pkl', 'w') as f:
   pickle.dump(pdic,f)

 f=open(path+basename+'.obs','w')
 for p in list(iobs): 
  f.write('{:<12}{:>12}{:>12}\n'.format(p,p[0],p[1]))
   
 f.close()   

#PLOT
 count=-1
 for p in loci:
  m.plot(p[0],p[1],'ro')   
  count=count+1
  plt.annotate(count,xy=(p[0],p[1]), xytext=(-15,15), textcoords='offset points', size='large', color='r', label='ACTUAL')

 count=-1
 for l1,l2 in iobs:
  count=count+1
  xx=grd.x.T[int(l1)-1,int(l2)-1] # fortran/python conversion
  yy=grd.y.T[int(l1)-1,int(l2)-1]
  m.plot(xx,yy,'k+',markersize=10,linewidth=6, label='OBS')
# plt.annotate(count,xy=(xx,yy), xytext=(-5,5), textcoords='offset points', size='large', color='b')

 m.drawcoastlines(linewidth=1.5)
 m.drawparallels(parallels)
 m.drawmeridians(meridians)

 plt.show(block=False)



if __name__ == "__main__":
# read grd file
    grd=Grid.fromfile('/DATA/HWRF/tmp/hwrf.grd')
# read dep file
    bath=Dep.read('/DATA/HWRF/tmp/hwrf.dep',grd.shape)
    bath.val=bath.val.T

    lon0=grd.x.min()
    lon1=grd.x.max()
    lat0=grd.y.min()
    lat1=grd.y.max()

    createf(path,basename,lat0,lat1,lon0,lon1,grd,bath)

