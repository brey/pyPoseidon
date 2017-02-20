import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from scipy.interpolate import interp1d, splrep, splev
from family import family
from xmlreader import xmlr


def vtcalc(time,lat,lon,vmax,inpfile,outdirplot):
  """ 
  Input from bulletin: time array, latitude, longitute, maximum velocity for each time
  Output : 

  """
  radius=6378388 #137  # earth's radius according to WGS 84
  deg2m=np.pi*radius/180 # ds on cicle equals ds=r*dth <-- dth=pi/180 
  
  if lon == None :
  
     #read $inpfile1 tab
     try:
          inpdat=xmlr(inpfile)
          time=inpdat.time
          lat=inpdat.lat
          lon=inpdat.lon
          vmax=inpdat.vmax
          ne64=inpdat.ne64
          se64=inpdat.se64
          sw64=inpdat.sw64
          nw64=inpdat.nw50
          ne50=inpdat.ne50
          se50=inpdat.se50
          sw50=inpdat.sw50
          nw50=inpdat.nw50
          ne34=inpdat.ne34
          se34=inpdat.se34
          sw34=inpdat.sw34
          nw34=inpdat.nw34
     except:
          txtfile=os.path.splitext(inpfile)[0]+'.txt'
          inpdat=readtxt(txtfile)
          time=inpdat.time
          lat=inpdat.lat
          lon=inpdat.long
          vmax=inpdat.vmax
          ne64=inpdat['64ne']
          se64=inpdat['64se']
          sw64=inpdat['64sw']
          nw64=inpdat['50nw']
          ne50=inpdat['50ne']
          se50=inpdat['50se']
          sw50=inpdat['50sw']
          nw50=inpdat['50nw']
          ne34=inpdat['34ne']
          se34=inpdat['34se']
          sw34=inpdat['34sw']
          nw34=inpdat['34nw']
#    time,lat,lon,vmax,ne64,se64,sw64,nw64,ne50,se50,sw50,nw50,ne34,se34,sw34,nw34,notes=np.loadtxt(inpfile,skiprows=1).T
   #exec longSign.pr
     with open ('../longSign.pr', 'r') as myfile:
        data=myfile.read()
     sig1=np.int(data[4::])
  
     sig=sign(lon)
  
     m=sig != sig1
     if np.sum(m) != 0 : 
 
        if (sig1 > 0): lon[m] = lon[m]+360
        if (sig1 < 0): lon[m] = lon[m]-360
  
# tsec=np.array(time)*3600  #translate time from hours to sec
  tsec=time*3600  #translate time from hours to sec
  nv=np.size(tsec)
  
  t1000=np.linspace(tsec[0],tsec[-1],1000)  # create an array of equally spaced points that spans the time [0:1000:120]
  
  x=lon
  y=lat
 
#------------------------------------------------------------------------------------------------------
# linear (if we have less or equal that 2 points, or spline interpolation on the array created above
#------------------------------------------------------------------------------------------------------
  
  if (nv <= 2):
        fx=interp1d(tsec,x)
        fy=interp1d(tsec,y)
        x1000 = fx(t1000)
        y1000 = fy(t1000)
  elif nv==3 :
        fx=interp1d(tsec,x,kind='quadratic')
        fy=interp1d(tsec,y,kind='quadratic')
        x1000 = fx(t1000)
        y1000 = fy(t1000)
  else: 
        tck=splrep(tsec,x,s=0)
        x1000 = splev(t1000,tck,der=0)
        tck=splrep(tsec,y,s=0)
        y1000 = splev(t1000,tck,der=0)
  
# compute the arc length s for all the points (time), that is the trajectory in time

  s1000 = np.cumsum(np.sqrt(np.ediff1d(x1000)**2+np.ediff1d(y1000)**2))
  s1000 = np.insert(s1000,0,0)  # adding 0 for the first point
  
# plot the initial points and the interpolated tranjectory

  plt.figure()
  plt.plot(x,y,'x',x1000,y1000,'-')
  plt.xlabel('lon')
  plt.ylabel('lat')
  plt.savefig(outdirplot+'/track.png')
# show(block=False)
  
  dt=tsec[1]-tsec[0]  # time increment, depends on the number of points, default 12 hours
  t1=np.maximum(tsec-dt,tsec[0]) # extrapolate backwards (or zero for bulNo 1) for a dt
  t2=np.minimum(tsec+dt,tsec[-1]) #extrapolate forward for dt
  
  
  #define  a linear interpolation function s(t) 
  [fs]=family(t1000,s1000)
  s1=fs(t1)  # compute values for points t1 
  s2=fs(t2)  # compute values for points t2
  
  #define  a linear interpolation function x(t) 
  [fx]=family(t1000,x1000)
  x1=fx(t1)
  x2=fx(t2)
  
  
  #define  a linear interpolation function y(t) 
  [fy]=family(t1000,y1000)
  y1=fy(t1)
  y2=fy(t2)
  
  t1h=t1/3600 #translate to hours
  t2h=t2/3600
  
  cincl=(x2-x1)/(s2-s1) # compute cos(phi)=dx/ds on the tranjectory
  sincl=(y2-y1)/(s2-s1) # compute sin(phi)=dy/ds on the tranjectory
  vt=(s2-s1)/(t2-t1)  # compute translational velocity ds/dt
  
  vtrx = vt * cincl *deg2m*np.cos(np.radians(lat)) #compute translational vx=vt*cos(phi) adjusted for latitude
  vtry = vt * sincl *deg2m
  
  vtr = np.sqrt(vtrx**2+vtry**2)
  
  time1=time
  
#------------------------------------------------------------------------------------------------------
#  save to txt
#------------------------------------------------------------------------------------------------------
  var=np.column_stack([time1,vtr,cincl,sincl,vtrx,vtry])
  svar=['time1','vtr','cincl','sincl','vtrx','vtry']
  he="\t\t".join(svar)
  fmt="\t".join(["%13.6g"]*(np.shape(var)[1]))
#  print repr(svar[i]).rjust(25)
# val=np.c_[time1,vtr,cincl,sincl,vtrx,vtry]
# for row in val:
#   print '%s' % (' '.join('%10.8s' % i for i in row))
  np.savetxt('vtcal.txt',var,header=he, fmt=fmt, comments='\t')
  
  t1=(tsec[:-1]+tsec[1:])/2
  
# compute .....
  
  p1=np.ediff1d(x)/np.ediff1d(tsec)*deg2m*np.cos(np.radians(y[1:]-np.ediff1d(y)/2)) 
  p2=np.ediff1d(y)/np.ediff1d(tsec)*deg2m 

#------------------------------------------------------------------------------------------------------
  # plot  

  plt.figure()
  colors=['r','b','g','m']
  ax=plt.subplot(111)
  l1=plt.plot(tsec,vtrx,'r',label='vt * cincl *deg2m*cos(lat)')
  l2=plt.plot(t1,p1,'b',label='delta(x)/delta(tsec)*deg2m*cos(y(2:)-delta(y)/2)')
  l3=plt.plot(tsec,vtry,'g',label='vt * sincl *deg2m')
  l4=plt.plot(t1,p2,'m',label='delta(y)/delta(tsec)*deg2m')
  leg=plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=1, mode="expand", borderaxespad=0.)
  for color,text in zip(colors,leg.get_texts()):
       text.set_color(color)
  
  #### getting the x axis label in sci notation
  xma=np.ma.masked_equal(tsec,0)
  np.ma.set_fill_value(xma, 999999.)
  sc=np.floor(np.log10(np.abs(xma))).astype(int)
  scale_pow=-np.max(sc)
  def my_formatter_fun(x, p):
      return "%.2f" % (x * (10 ** scale_pow))
  
  ax.get_xaxis().set_major_formatter(ticker.FuncFormatter(my_formatter_fun))
  ax.set_xlabel('x 10^{0}'.format(-scale_pow))
  ax.xaxis.set_label_coords(1.05, -0.25)
  
  
  plt.tight_layout(rect=[.02,.02,.98,.8])
  plt.savefig(outdirplot+'/vt.png')
# plt.show(block=False)

#------------------------------------------------------------------------------------------------------
  return vt,vtr,vtrx,vtry,cincl,sincl
  
  
if __name__ == '__main__':

    path='../storage/tmp/1000130/1/'
    time,lat,lon,vmax=np.loadtxt(path+'inpData.txt',skiprows=1,usecols=(0,1,2,3)).T
    vtcalc(time,lat,lon,vmax)
    
