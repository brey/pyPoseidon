import sys
import numpy as np
import matplotlib.pyplot as plt
from xmldic import bxml
from xmlreader import xmlr
from bunch import *
from family import family
from vtcalc import vtcalc
from read import readtxt
#from holsq import lsq
#from slsqp import slsqp
#from hroot import hroot
#from visvel import visvel
import os



# -------------------------------------------------------------------------
# Const
# -------------------------------------------------------------------------
nm2m=1852. # 1 nautical mile to meters
kt2ms=nm2m/3600.  # knots to m/s
omega=2*np.pi/(3600.*24.) # angular speed omega=2pi*f(=frequency of earth : 1 cycle per day) 2pi* 1 / day in seconds
rhoa=1.15 #air density  Kg/m^3
radius=6378388 #137. # earth's radius according to WGS 84
deg2m=np.pi*radius/180.  # ds on cicle equals ds=r*dth - dth=pi/180
pn=101000.  # Atmospheric pressure [N/m^2] (101KPa - enviromental pressure)

tetaNE=45. #mean angle [degrees] of quadrants
tetaNW=135. #        "
tetaSW=225. #        "
tetaSE=315. #        "

npmin=2

kmin=0  # low limit of parameter k (=xn-.5-> k=0-> x=0.5)
kmax=0.15 # upper limit for k (->xn=.65)  WHY?

dpmin=10.e2  # minimum value of  pressure drop P_central - P_env(=101kPa).
dpmax=200.e2   # maximum value of  pressure drop P_central - P_env(=101kPa).
rvmaxmin=10.e3  # default minimum value of Rmax[m] 

bmin=0.8 # minimum value of holland parameter b
#bmax=2.5
bmax=1.8  # maximum value of holland parameter b
b0=1.2  # initial estimation of holland parameter b

rmax0=20.e3  # intial estimation for radius of maximum wind [m] (20km)
maxR=500.e3  # maximum radius of TC [m] (500Km)


def expf(t,A,K,C):
	return A*np.exp(K*t)+C



def wr2h(infofile,inpfile,outfile):

# -------------------------------------------------------------------------
# Input/Output File
# -------------------------------------------------------------------------

  info=bxml(infofile)

  try:
     inpdat=xmlr(inpfile)
     time=inpdat.time
     lat=inpdat.lat
     lon=inpdat.lon
     vmax=inpdat.vmax
     ne64=inpdat.ne64
     se64=inpdat.se64
     sw64=inpdat.sw64
     nw64=inpdat.nw64
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
     nw64=inpdat['64nw']
     ne50=inpdat['50ne']
     se50=inpdat['50se']
     sw50=inpdat['50sw']
     nw50=inpdat['50nw']
     ne34=inpdat['34ne']
     se34=inpdat['34se']
     sw34=inpdat['34sw']
     nw34=inpdat['34nw']
#time,lat,lon,vmax,ne64,se64,sw64,nw64,ne50,se50,sw50,nw50,ne34,se34,sw34,nw34,notes=np.loadtxt(inpfile,skiprows=1).T


#-------------------------------------------------------------------------

  bulNo=np.int(info.setexp.bulNo)

#nb=$n

  nb=10000  #size of random numbers used

  fk=0.92  # coefficient for going from 1m to 10m in velocities ????????????????????????
  fvtr=1.0

#fk=1.0
#fvtr=0.01

  stitle=info.setexp.source
  hurName=info.setexp.hurName
  bulDate=info.setexp.bulDate

# -------------------------------------------------------------------------
# Read input File
# -------------------------------------------------------------------------


#sig1=np.loadtxt('../longSign.pr')

  sig=np.sign(lon)  # get the sign of longitute 
  sig1=sig[0] # uses the first sign to set the adjustment below. 

  m=sig != sig1 # map the values of lon that have sign different from the first lon

  if sum(m)>0:
# adjust the lon values going from -180:180
	if sig1 > 0:
		lon[m]=lon[m]+360.
	elif sig1 < 0:
		lon[m]=lon[m]-360.

# storing the input quantities

  TIME =  time  
  LAT = lat
  LON = lon
  VMAX0 = vmax
  ntime=np.size(TIME) # number of time steps from bulletin


# initialize arrays
  deltaptot = np.zeros(ntime)
  btot = np.zeros(ntime)
  ktot=np.zeros(ntime)
  rmaxtot = np.zeros(ntime)
  vmax0tot = np.zeros(ntime)
  vmax0ktot = np.zeros(ntime)
  vmax1tot = np.zeros(ntime)
  vtrtot = np.zeros(ntime)
  vtrxtot = np.zeros(ntime)
  vtrytot = np.zeros(ntime)
  rmsetot=np.zeros(ntime)
  biastot=np.zeros(ntime)



#==========================================================================
# vt
#==========================================================================
  try:
    if ntime > 1:   # if we have a bulletin with time array ...
     try:
        vt,vtr,vtrx,vtry,cincl,sincl=vtcalc(time,lat,lon,vmax,inpfile,'./') #get translational velocity
     except: 
        print 'VTCALC FAILED'
        sys.exit()

#  why use a coef fvtr ???????????????
     vtr=vtr*fvtr
     vtrx=vtrx*fvtr
     vtry=vtry*fvtr

# store variables
     SINFI = sincl
     COSFI = cincl
     VTR = vtr
     VTRX = vtrx
     VTRY = vtry

    else: 

# extrapolate from previous bulletin

	time1,sincl,cincl,vtr,vtrx,vtry = np.load.txt('../'+np.str(bulNo-1)+'/vtcal.txt').T 
 	[sincl,cincl, vtr, vtrx, vtry] =family(time1,sincl,cincl,vtr,vtrx,vtry)

	SINFI = sincl(time+6)
	COSFI = cincl(time+6)
	VTR = vtr(time+6)
	VTRX = vtrx(time+6)
	VTRY = vtry(time+6)



  except:   
        print 'setting  default values for vtcalc'
	SINFI = 1
	COSFI = 0
	VTR = 0
	VTRX = 0
	VTRY = 0

#==========================================================================
# loop over all times in bulletin
#==========================================================================

  done=0
  for t in range(ntime):

   time = TIME[t]
   lat = LAT[t]
   lon = LON[t]
   vmax0 = VMAX0[t]
   vmax0k = vmax0*fk
   V0=np.array([64, 64, 64, 64, 50, 50, 50, 50, 34, 34, 34, 34])*kt2ms*fk #translate knots to m/s
   R0=np.array([ne64[t], se64[t], sw64[t], nw64[t], ne50[t], se50[t], sw50[t], nw50[t], ne34[t], se34[t], sw34[t], nw34[t]])*nm2m #translate from nautical miles to m
   an=np.array([tetaNE, tetaSE, tetaSW, tetaNW,tetaNE, tetaSE, tetaSW, tetaNW,tetaNE, tetaSE, tetaSW, tetaNW])
   vtr = VTR[t]
   vtrx = VTRX[t]
   vtry = VTRY[t]
   sinfi = SINFI[t]
   cosfi = COSFI[t]  

#--------------------------------------------------
# Choose valid velocity 
#--------------------------------------------------
   M=(R0 != 0) & (R0 < maxR)  # find which radii are not zero and less than  Rmax(=500km)
   if  np.sum(M) < npmin : # if there are less than 2 

   	if (done) :  # if not the time[0] use the previous values
   		rmax = rmaxtot[t-1]
  		b = btot[t-1]
   		k = ktot[t-1]

   	else:
   		rmax = rmax0
  		b = b0
   		k = kmin

  	bias=None
   	rmse=None

   	dp = vmax0**2*rhoa*np.exp(1)/b
   	deltaptot[t] = dp
   	btot[t] = b
   	rmaxtot[t] = rmax 
   	vmax0tot[t] = vmax0
   	vmax0ktot[t] = vmax0k
   	vmax1tot[t] = np.max([vmax0k-vtr, vmax0k/2])
   	vtrtot[t] = vtr
   	vtrxtot[t] = vtrx
   	vtrytot[t] = vtry
   	ktot[t] =k
   	biastot[t]=bias
   	rmsetot[t] = rmse
        continue

   V0 = V0[M]
   R = R0[M]
   an = an[M]
   sinan = np.sin(np.radians(an+90))  # an +90 = angle of tangential wind
   cosan=np.cos(np.radians(an+90))

   npv=np.size(V0)

#--------------------------------------------------
# calculate V wind radii (V - translational velocity)
#--------------------------------------------------
   RATIO = (rmax0/R)**b0    # assume exponential decay eqs (13) from JRC report
   EXPRATIO = np.exp(-RATIO)  #                       "

   VT=vtr*(cosfi * cosan + sinfi * sinan)*(1-EXPRATIO)   # Eq (15) from JRC report
#   VT=vtr*(cosfi * cosan + sinfi * sinan)*(RATIO)

   if (lat<0): VT=-VT   # reverse for south hemishpere
   VV = V0-VT   # substract translational velocity from TC velocity

   # calculate  f WR

   deltalatWR=R/deg2m*np.sin(np.radians(an))
   latWR=lat+deltalatWR

   fWR=2*omega*np.abs(np.sin(np.radians(latWR))) # Coriolis parameter f=2*Omega*sin(lat)
   Vnco=((VV+R*fWR/2)**2-(R*fWR/2)**2)**0.5
   V=Vnco


#--------------------------------------------------
# Choose valid velocities
#--------------------------------------------------
   M=V>0

   if np.sum(M) == 0 : continue

   V0 = V0[M]
   VV = VV[M]
   V=V[M]
   R=R[M]

#--------------------------------------------------
#  calculate  vmax = vmax0 - vt
#--------------------------------------------------
   vmax0vt=np.max([vmax0k-vtr,np.max(V)])


#--------------------------------------------------   
#  calculate coriolis for vmax
#--------------------------------------------------
   if (lat>0): 
      #sinfivmax=np.sin(fi-90)=-np.cos(fi)
      sinfivmax=-cosfi
   else:
      #sinfivmax=np.sin(fi+90)=np.cos(fi)
      sinfivmax=cosfi


#  fitf1,p1=lsq(R,V,maxR,rhoa,b0,rvmaxmin,dpmin,kmin,vmax0vt)
#  fitf2,p2,p3=slsqp(R,V,maxR,rhoa,b0,rvmaxmin,dpmin,kmin,vmax0vt)
#  hroot(R,V,maxR,rhoa,b0,rvmaxmin,dpmin,kmin,vmax0vt)
#  visvel(R,V,maxR,rhoa,b0,rvmaxmin,dpmin,kmin,vmax0vt)
#===========================================================================
# Monte Carlo method for estimating Rmax
#===========================================================================
   # k

   if (npv>npmin):
     K=kmin+(kmax-kmin)*np.random.rand(nb)
   elif(done): 
     K = np.ones(nb)*ktot[t-1]
   else:
     K = np.ones(nb)*kmin

   #  DP
   DP=dpmin*(dpmax/dpmin)**np.random.rand(nb)   

   #  Rmax
   rvmaxmin_=np.min([rvmaxmin,np.min(R)*0.5])  # update the minimum  value for Rmax with the R.min/2 from input
   RMAX=rvmaxmin_*(np.min(R)*0.99/rvmaxmin_)**np.random.rand(nb) # range  min(10000,Rmin/2)<Rmax<.99*Rmin (scaled)

#--------------------------------------------------
# calculate vmax1 = v max0k -vt - Coriolis effect (function of RMAX)
#--------------------------------------------------
   deltalatvmax=RMAX/deg2m*sinfivmax  # for each Rmax we compute the lat deviation for the velocity
   latvmax=lat+deltalatvmax

   fvmax=2*omega*np.abs(np.sin(np.radians(latvmax))) # Coriolis coef f  

   fvmax2=RMAX*fvmax/2
   vmax1=((vmax0vt+fvmax2)**2-fvmax2**2)**0.5
   mask=vmax1<np.max(V)
   np.copyto(vmax1,np.max(V),where=mask)
 
#----------------------------------
# use the random values of vmax,dp above we compute b (from Holland 2010 - eqs (7))

   B=(rhoa*np.exp(1)/DP)*vmax1**2  


   m=(B >= bmin) & ( B <= bmax) & (lat*latvmax > 0)  # mask B that fits all 3 criteria
   nb1 = np.sum(m) #number of 'True' values

#  mask arrays accordingly 
   K=K[m]
   DP=DP[m]
   RMAX=RMAX[m]
   B = B[m]

   nval = np.size(V)  # number of V > 0 
   Vcalc = []
   RMS = np.zeros(nb1)

   try:
      	r 
   except NameError:
        pass
   else: 
	r=None


# check values for all V
   for i in range(nval):
  	r = R[i]
  	ratio=(r-RMAX)/(maxR-RMAX)
  	X=0.5 + np.min([np.max(ratio,0),1])*K   #compute x using random k  & Rmax
  	Vcalc=np.append(Vcalc,((B/rhoa) * DP* (RMAX/r)**B * np.exp(-(RMAX/r)**B))**X)  # compute & store V

   for i in range (nb1):
  	RMS[i]=np.sqrt(np.average((Vcalc[i::nb1]-V)**2))  # compute deviation from estimated and given values

   value=nb1
   totvalue=nb

# -------------------------------------------------------------------------
# select final velocities
# -------------------------------------------------------------------------
   m=RMS == np.min(RMS)  #find minimum RMS

# select the minimizing quantities
   rmse=RMS[m][0]
   dp=DP[m][0] 
   b=B[m][0]
   rmax=RMAX[m][0]
   k=K[m][0]

   vmax1 = np.sqrt(b*dp/(rhoa*np.exp(1)))  # compute estimated vmax


# print on screen
   var=[t,time,rmse,dp,b,rmax,k ,vmax0,vmax0k,vmax0vt,np.max(VV),np.max(V),vmax1]
   varn=['t','time','rmse','dp','b','rmax','k' ,'vmax0','vmax0k','vmax0vt','np.max(VV)','np.max(V)','vmax1']
#  for i in range(np.size(var)):
#	   print repr(varn[i]).rjust(25), repr(var[i]) 


#******************************************************************************
# PLOT Vg
#******************************************************************************
   Rmax=rmax

   rx=np.arange(100)
   ex=np.log(max(R)*1.2/1000.)/100.
   r=expf(rx,1000.,ex,0.) 

   ratio=(r-rmax)/(maxR-rmax)
   
   mask=ratio<0.
   np.copyto(ratio,0.,where=mask)
   mask=ratio>1.
   np.copyto(ratio,1.,where=mask)
  
   x=0.5 + ratio*k

   Vg=( (b/rhoa) * dp* (Rmax/r)**b * np.exp(-(Rmax/r)**b) )**x 

   [fVg]=family(r,Vg)

   bias=np.average(fVg(R)-V)
   rmse=np.sqrt(np.average((fVg(R)-V)**2))

#  READ DATA FROM VPL RUN FOR COMPARISON
   name = '%2.0f' % time
#  Rvpl, Vvpl=np.loadtxt('../storage/tmp/1000130/1/d_R'+name+'.',skiprows=1).T
#  rvpl, Vgvpl=np.loadtxt('../storage/tmp/1000130/1/d_r'+name+'.',skiprows=1).T

#  plt.figure()
#  plt.plot(R/1000,V,'kx',r/1000,Vg,'k--')#,r/1000,fitf1(p1,r),'c--',r/1000,fitf2(p2,r),'g--',r/1000,fitf2(p3,r),'b--',Rvpl/1000,Vvpl,'ro',rvpl/1000,Vgvpl,'r-')
#  plt.axis([r.min()/1000,500.,0.,60.])
#  plt.xlabel('Radius [km]')
#  plt.ylabel('Velocity (m/s)')
#  plt.title(stitle)
#  plt.figtext(.5,.05,'TIME='+np.str(time)+'     Size of random sample:'+np.str(nb),horizontalalignment='center')
#  plt.figtext(.5,.03,'B='+'%.4f' % b+' Rmax='+'%.4f' % (Rmax/1000)+' [km]  DeltaP='+'%.4f' % (dp/100)+' [mBar] k='+'%.4f' % k+' Vmax_adv='+'%.4f' % vmax0+' [m/sec]',horizontalalignment='center',size='small') 
#  plt.figtext(.5,.01,'rmse='+'%.4f' % rmse+' bias='+'%.4f' % bias+' value(\%)='+'%.4f' %(value/totvalue*100)+' Vt='+'%.4f' % vtr+' [m/sec] Vmax_Hol='+'%.4f' % vmax1+' [m/sec]',horizontalalignment='center',size='small') 
#  plt.subplots_adjust(bottom=.15)
#  plt.savefig(outdirplot+'/vg_'+np.str(time)+'.png')

#******************************************************************************
   done = 1  # first time is done
# -------------------------------------------------------------------------
#  store variables
   deltaptot[t] = dp
   btot[t] = b
   rmaxtot[t] = rmax 
   vmax0tot[t] = vmax0
   vmax0ktot[t] = vmax0k
   vmax1tot[t] = vmax1
   vtrtot[t] = vtr
   vtrxtot[t] = vtrx
   vtrytot[t] = vtry
   ktot[t] =k
   biastot[t]=bias
   rmsetot[t] = rmse

# -------------------------------------------------------------------------
#  loop in time ends here
# -------------------------------------------------------------------------



  l=btot != -1  # map non negative values of b


# map arrays accordingly
  time=TIME[l]
  xhc=LON[l]
  yhc=LAT[l]
  deltap=deltaptot[l]
  rmax=rmaxtot[l]
  b = btot[l]
  k = ktot[l]
  vmax0 = vmax0tot[l]
  vmax0k = vmax0ktot[l]
  vmax1 = vmax1tot[l]
  vtr = vtrtot[l]
  vtrx = vtrxtot[l]
  vtry = vtrytot[l]
  bias=biastot[l]
  rmse = rmsetot[l]


#******************************************************************************
# OUTPUT file
#******************************************************************************

  vmax=vmax0  #  original vmax as in inpData.txt
  vmax0=vmax1   # attenzione, in bulAdr.pr si aspetta che vmax0 sia vmax a 10 min average ,  meno coriolis e vel translazione

  var=np.column_stack([time,xhc,yhc,b,k,rmax,deltap,vmax,vmax0,vtr,vtrx,vtry,bias,rmse])
  svar=['time','xhc','yhc','b','k','rmax','deltap','vmax','vmax0','vtr','vtrx','vtry','bias','rmse']
  he="\t\t".join(svar)
  fmt="\t".join(["%13.6g"]*(np.shape(var)[1]))

  np.savetxt(outfile,var,header=he, fmt=fmt, comments='\t')


if __name__ == '__main__':
#===========================================================================
# INPUT DATA
#===========================================================================
 infofile=sys.argv[1]
 inpfile=sys.argv[2]
 outfile=sys.argv[3]
 print infofile,inpfile,outfile
 wr2h(infofile,inpfile,outfile)
