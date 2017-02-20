import numpy as np
import datetime
import sys
from netCDF4 import Dataset
from grid import *
from dep import *
import itertools


#path0='/home/critechuser/2016/'

def get(t0,t1,path0,basename,plat,plon):

  plat=float(plat)
  plon=float(plon)

# print plon,plat

  tstart=datetime.datetime.strptime(t0,'%Y%m%d.%H')
  tend=datetime.datetime.strptime(t1,'%Y%m%d.%H')
  dt=(tend-tstart).total_seconds()
  ndt=dt/(3600*12)
  ndt=np.int(ndt)+1


  path=path0+'{}/'.format(t0)

  # read computational grid
  grid = Grid.fromfile(path+basename+'.grd')
  lon=grid.x[0,:].data
  lat=grid.y[:,0].data
  bath=Dep.read(path+basename+'.dep',grid.shape)
  bath.val=bath.val.T[:-1,:-1]
  dx=lon[1]-lon[0]
  dy=lat[1]-lat[0]
  x=lon-dx/2.
  y=lat-dy/2.


  i=np.abs(lon-plon).argmin()
  j=np.abs(lat-plat).argmin()
  i0=i
  j0=j
# print i0,j0
  x0=grid.x.T[i,j]
  y0=grid.y.T[i,j]
# print i,j,x0,y0,bath.val[i,j]
  # retrieve the 4 nearest diagonal points of the i,j
  lon4=grid.x.T[i-1:i+3:2,j-1:j+3:2]
  lat4=grid.y.T[i-1:i+3:2,j-1:j+3:2]
  val4=bath.val[i-1:i+3:2,j-1:j+3:2]
# print zip(lon4.flatten(),lat4.flatten(),val4.flatten())
# print '==================='
  lon8=grid.x.T[i-1:i+2,j-1:j+2]
  lat8=grid.y.T[i-1:i+2,j-1:j+2]
  val8=bath.val[i-1:i+2,j-1:j+2]
# print zip(lon8.flatten(),lat8.flatten(),val8.flatten())
# print '==================='
  # define the quandrant
  A=plon-x0
  B=plat-y0
  for xx,yy in zip(lon4.flatten(),lat4.flatten()):
#        print xx,yy
         C=np.sign([xx-plon,A])
         D=np.sign([yy-plat,B])
         if C[0]==C[-1] and D[0] == D[-1] : 
           corner=[xx,yy]
           indx=[np.abs(lon-xx).argmin(),np.abs(lat-yy).argmin()]

# print indx
# print [i,j,x0,y0,bath.val[i,j]]
# print [i,indx[1],x0,corner[1],bath.val[i,indx[1]]]
# print [indx[0],indx[1],corner[0],corner[1],bath.val[indx[0],indx[1]]]
# print [indx[0],j,corner[0],y0,bath.val[indx[0],j]]

  bath4=[]
  for r in itertools.product([indx[0],i],[indx[1],j]): bath4.append(bath.val[r[0],r[1]])
# print bath4
  if np.isfinite(bath4).any() : 

      i=np.abs(x-plon).argmin()
      j=np.abs(y-plat).argmin()
 #    print i,j

  else:

     try:
  
        l=np.argwhere(np.array(val8)==np.nanmin(np.array(val8)))
        [ii,jj]=l.ravel()
    
        wx=lon8[ii,jj]
        wy=lat8[ii,jj]
        i=np.abs(x-wx).argmin()
        j=np.abs((y-wy).T).argmin()
  
# [i,j]=indx

     except:
        i,j=i0,j0

  combined=[]
  tw=[]

  for it in range(ndt-1): 
    idate=tstart+datetime.timedelta(hours=12*it)
    idt=datetime.datetime.strftime(idate,'%Y%m%d.%H')
    path=path0+'{}/'.format(idt)

 #  print path
    d = Dataset(path+'trim-'+basename+'.nc')
   #h=d.variables['S1'][:12,i:i+2,j:j+2]
    h=d.variables['S1'][:12,i,j]
    

    time=d.variables['time'][:12]

    tm=(idate-tstart).total_seconds()+(time-time[0])

    tstamp=[]
    for l in tm : tstamp.append(tstart+datetime.timedelta(0,int(l)))

    tw.append(tstamp)
    combined.append(h)#.mean(axis=(1,2)))

  tcw=np.array(tw).flatten()
  cw=np.array(combined).flatten()
# Add the last one with forecasting

  idate=tend
  idt=datetime.datetime.strftime(idate,'%Y%m%d.%H')
  path=path0+'{}/'.format(idt)
# print path

  d = Dataset(path+'trim-'+basename+'.nc')
 #h=d.variables['S1'][:,i:i+2,j:j+2]
  h=d.variables['S1'][:,i,j]
    

  time=d.variables['time'][:]

  tm=(idate-tstart).total_seconds()+(time-time[0])

  tstamp=[]
  for l in tm : tstamp.append(tstart+datetime.timedelta(0,int(l)))

  tcw=np.append(tcw,np.array(tstamp))
  cw=np.append(cw,np.array(h))#.mean(axis=(1,2))))

# w1=tcw<datetime.datetime.strptime('20160228.00','%Y%m%d.%H')
# cwm = cw[w1].mean()
  
  return tcw,cw,lat[j]-dy/2.,lon[i]-dx/2.,j,i

if __name__ == "__main__":
    tstart=sys.argv[1]
    tend=sys.argv[2]
    path=sys.argv[3]
    basename=sys.argv[4]
    plat =sys.argv[5]
    plon =sys.argv[6]
    t,ha,mlat,mlon,jm,im=get(tstart,tend,path,basename,plat,plon)

