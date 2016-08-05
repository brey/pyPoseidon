import numpy as np
import pandas
import collections
import pickle
from itertools import product

from dep import *
from grid import *

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid


def createf(path,basename,lat0,lat1,lon0,lon1,grd,bath):

# bathymetry
 bath.val=bath.val.T

# grid
 dx=grd.x[0,1]-grd.x[0,0]
 dy=grd.y[1,0]-grd.y[0,0]

#  BUOYS WEBCRITECH
 dat=pandas.read_csv('SeaLevelBuoys2.csv')
 
 dat.columns=np.append(['ID'],dat.columns.get_values()[1:]) # fix ID header
#print dat.columns

 ID=dat['ID']
 lon=dat['lon']
 lat=dat['lat']
 names=dat[['ID','name']]
 names=names.set_index(['ID'])
 dnames=names.T.to_dict('list')
 
#------------------------------------------------------------------
# mapping
 fig1 = plt.figure(figsize=(10,8))
 ax = fig1.add_axes([0.1,0.1,0.8,0.8])

 m = Basemap(projection='cyl',llcrnrlat=lat0,urcrnrlat=lat1,\
             llcrnrlon=lon0,urcrnrlon=lon1,resolution='l')

# define parallels and meridians to draw.
 parallels = np.arange(-90.,90,20.)
 meridians = np.arange(0.,360.,20.)
#------------------------------------------------------------------

 iobs=[] # initialize
 loci=[]
 inames=[]
 k=0
 pdic={}
 for l1,l2 in zip(lon,lat): # all points in database
  if (lon0 <= l1 <= lon1) & (lat0 <= l2 <= lat1) :
    ih=np.abs(grd.x-l1).argmin()
    jh=np.abs((grd.y-l2).T).argmin()
# plot first guess
    xx=grd.x.T[ih,jh]-dx/2. 
    yy=grd.y.T[ih,jh]-dy/2.
    m.plot(xx,yy,'gx',markersize=10,linewidth=6)
    pi=ih+1 # the fortran/python index issue ??
    pj=jh+1 # the fortran/python index issue ??

# define the nearby points
    lon8=grd.x.T[ih-1:ih+2,jh-1:jh+2]
    lat8=grd.y.T[ih-1:ih+2,jh-1:jh+2]
    val8=bath.val[ih-1:ih+2,jh-1:jh+2]
# choose maximum depth
    l=np.argwhere(np.array(val8)==np.nanmax(np.array(val8)))
# if l exists
    if len(l) > 0:

     [i,j]=l.ravel()
    
     x=lon8[i,j]
     y=lat8[i,j]
     ih=np.abs(grd.x-x).argmin()
     jh=np.abs((grd.y-y).T).argmin()
     pi=ih+1 # the fortran/python index issue ??
     pj=jh+1 # the fortran/python index issue ??

    else:
   #  print l,l1,l2,val8
      continue
                 
    m1=np.argwhere(lon==l1) 
    m2=np.argwhere(lat==l2) 
    if m1.all()==m2.all():
     if (pi,pj) not in iobs:
     #print '--->',pi,pj
      for m0 in m1.flatten():
      # print 'm0=', m0
        pdic[ID[m0]]=k
        if ID[m0] not in inames and (pi,pj) not in iobs:
           iobs.append((pi,pj))
           loci.append((l1,l2))
 #         inames.append(names.ix[ID[m0],'name'])
           inames.append(ID[m0])
        else:
           k=k-1
           pdic[ID[m0]]=k
        k=k+1
     else:
    # print pi,pj
      g=np.where(np.array(iobs)==[pi,pj])
      c=collections.Counter(g[0])
      rk = [value for value, count in c.items() if count > 1][0]
      for m0 in m1.flatten():
        pdic[ID[m0]]=rk
    else: print 'problem with {}, {}'.format(m1,m2)
    
 with open(path+basename+'.pkl', 'w') as f:
   pickle.dump(pdic,f)

 f=open(path+basename+'.obs','w')
 for p,n in zip(list(iobs),list(inames)): 
  f.write('{:<12}{:>12}{:>12}\n'.format(n,p[0],p[1]))
   
 f.close()   
 
#cnames = np.array([ '({},)'.format(l) for l in np.array(inames).ravel() ])
#inames=[ name[:15] for name in inames] 

#op=pandas.DataFrame(iobs,index=inames)
#op['C']=''
#opp= op[['C',0,1]]
#opp.to_csv('test',header=0,sep='\t')

#PLOT
 cs = m.contourf(grd.x,grd.y,bath.val[:-1,:-1].T,cmap=plt.cm.jet)
 count=-1
 for p in loci:
  m.plot(p[0],p[1],'ro')   
  count=count+1
  plt.annotate(count,xy=(p[0],p[1]), xytext=(-15,15), textcoords='offset points', size='large', color='r', label='ACTUAL')

 count=-1
 for l1,l2 in iobs:
  count=count+1
  xx=grd.x.T[int(l1)-1,int(l2)-1]-dx/2. # fortran/python conversion
  yy=grd.y.T[int(l1)-1,int(l2)-1]-dy/2.
  m.plot(xx,yy,'rx',markersize=15,linewidth=30, label='OBS')
# plt.annotate(count,xy=(xx,yy), xytext=(-5,5), textcoords='offset points', size='large', color='b')

 m.drawcoastlines(linewidth=1.5)
 m.drawparallels(parallels)
 m.drawmeridians(meridians)

 plt.show(block=False)

#return iobs,inames

if __name__ == "__main__":
    path='/DATA/critechuser/tmp2/'
    basename='med'
# read grd file
    grd=Grid.fromfile(path+basename+'.grd')
# read dep file
    bath=Dep.read(path+basename+'.dep',grd.shape)

    lon0=grd.x.min()
    lon1=grd.x.max()
    lat0=grd.y.min()
    lat1=grd.y.max()

    createf(path,basename,lat0,lat1,lon0,lon1,grd,bath)

