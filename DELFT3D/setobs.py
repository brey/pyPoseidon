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
 bval=bath.val.T[:-1,:-1]

# grid
 dx=grd.x[0,1]-grd.x[0,0]
 dy=grd.y[1,0]-grd.y[0,0]
 gx=grd.x[0,:]-dx/2. 
 gy=grd.y[:,0]-dy/2.

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
    ih=np.abs(gx-l1).argmin()
    jh=np.abs(gy-l2).argmin()
# plot first guess
    x0=grd.x.T[ih,jh]
    y0=grd.y.T[ih,jh]
    m.plot(x0-dx/2.,y0-dy/2.,'go',markersize=10,markerfacecoloralt='green',fillstyle=None)
    pi=ih+1 # the fortran/python index issue ??
    pj=jh+1 # the fortran/python index issue ??

    print dat[dat['lon']==l1]['ID'].values, ih,jh,bval[ih,jh]

    l=[]

   #if np.isnan(bval[ih,jh]) :
    if True :

  # retrieve the 4 nearest diagonal points of the i,j
      lon4=grd.x.T[ih-1:ih+3:2,jh-1:jh+3:2]
      lat4=grd.y.T[ih-1:ih+3:2,jh-1:jh+3:2]
      val4=bval[ih-1:ih+3:2,jh-1:jh+3:2]
#     print zip(lon4.flatten(),lat4.flatten(),val4.flatten())
# print '==================='
      lon8=grd.x.T[ih-1:ih+2,jh-1:jh+2]
      lat8=grd.y.T[ih-1:ih+2,jh-1:jh+2]
      val8=bval[ih-1:ih+2,jh-1:jh+2]
#     print zip(lon8.flatten(),lat8.flatten(),val8.flatten())
# print '==================='
  # define the quandrant
      A=l1-x0
      B=l2-y0
      for xx,yy in zip(lon4.flatten(),lat4.flatten()):
#        print xx,yy
         C=np.sign([xx-l1,A])
         D=np.sign([yy-l2,B])
         if C[0]==C[-1] and D[0] == D[-1] : 
           corner=[xx,yy]
           indx=[np.abs(gx+dx/2-xx).argmin(),np.abs(gy+dy/2-yy).argmin()]

      print indx
      print [ih,jh,x0,y0,bval[ih,jh]]
      print [ih,indx[1],x0,corner[1],bval[ih,indx[1]]]
      print [indx[0],indx[1],corner[0],corner[1],bval[indx[0],indx[1]]]
      print [indx[0],jh,corner[0],y0,bval[indx[0],jh]]
  
      bath4=np.array([bval[ih,jh],bval[ih,indx[1]],bval[indx[0],indx[1]],bval[indx[0],jh]])
      lb4=np.array([[ih,jh],[ih,indx[1]],[indx[0],indx[1]],[indx[0],jh]])

      print 'finite values', np.isfinite(bath4).sum()

      if np.isfinite(bath4).sum() > 1: 
          
   #      print l1,l2,bath4

          ih=np.abs(gx-l1).argmin()
          jh=np.abs(gy-l2).argmin()
          pi=ih+1 # the fortran/python index issue ??
          pj=jh+1 # the fortran/python index issue ??
        

      elif np.isfinite(bath4).sum() == 1:
# choose minimum depth
     #  l=np.argwhere(np.array(val8)==np.nanmin(np.array(val8)))
          knan=np.argwhere(np.isfinite(bath4)).flatten()
          print k
          [k1,k2]=lb4[knan[0]]
          print k1,k2

          i1=np.sign(k1-ih)
          i2=np.sign(k2-jh)

          ih=ih+i1  # opposite diagonal 
          jh=jh+i2

          print 'DIAGO', ih,jh
          
          pi=ih+1 # the fortran/python index issue ??
          pj=jh+1 # the fortran/python index issue ??

    else:
# recompute close i,j considering the staggered grid

          ih=np.abs(gx-l1).argmin()
          jh=np.abs(gy-l2).argmin()
          pi=ih+1 # the fortran/python index issue ??
          pj=jh+1 # the fortran/python index issue ??
        

# if l exists
    if len(l) > 0:

     [i,j]=l.ravel()
    
     x=lon8[i,j]-dx/2.
     y=lat8[i,j]-dy/2.
   # print x,y
     m.plot(x,y,'ko',markersize=5)
     ih=np.abs(gx-x).argmin()
     jh=np.abs((gy-y).T).argmin()
     pi=ih+1 # the fortran/python index issue ??
     pj=jh+1 # the fortran/python index issue ??
  #  print bval[ih,jh]

#   else:
#     continue
    print 'FINAL', dat[dat['lon']==l1]['ID'].values, ih,jh,bval[ih,jh]
                 
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
 #         print ID[m0], pi,pj,l1,l2,bval[pi-1,pj-1]
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
 cs = m.contourf(grd.x,grd.y,bval.T,cmap=plt.cm.jet)
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
    path='../../../tmp/20160101.00/'
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

