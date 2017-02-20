#######################################################
# Create Grid
#######################################################
# 	print "Syntax: numb_intervals cellsize_min HurName GDACS_ID bul_no "

#******************************************************
# IMPORT Python lib ...
#******************************************************
# matplotlib notebook
import pandas
import numpy as np
from scipy.interpolate import griddata
import sys
import os
import pp
import time

def create(minxc, maxxc, minyc, maxyc, ngrid, hurName,namefile,workdir,dcell):
	#******************************************************
	# SET FUNCTIONS
	#******************************************************
	# ----------------------------------------
	# f coriolis
	# ----------------------------------------
	def fcor(lat,omega):
		return 2*omega*np.sin(np.deg2rad(lat))
		
	# ----------------------------------------
	# Calculate x 
	# (!!!! da sistemare)
	# ----------------------------------------
	def calcx(r,rmax,maxR,k):   
		m=np.shape(r)	
		x=np.zeros(m)
		for i in range(0, np.shape(r)[0], 1):   
			for j in range(0, np.shape(r)[1], 1):
				ratio=(r[i][j]-rmax)/(maxR-rmax)
				ratio2=min(max(ratio,0),1)
				xi=0.5+ratio2*k          
				x[i][j]=xi           
		return x

	# ----------------------------------------
	# Inflow angle correction
	#  (not yet implemented in the time loop)
	# ----------------------------------------
	def calcbeta(r,rmax):

		mtot=len(r[1]) # da sistemare
	   
		beta=np.zeros(np.shape(r))
		for i in range(0, np.shape(r)[0], 1):   
			for j in range(0, np.shape(r)[1], 1):   
				if r[i][j]>(rmax*1.2):
					betai=25    
				elif r[i][j]>rmax:
					betai=85*(r[i][j]/rmax)-65    
				else: 
					betai=10*r[i][j]/rmax
				beta[i][j]=betai    
		return beta

	# ----------------------------------------
	# Pressure
	# ----------------------------------------
	def pres(r, dp, pn, rmax, b):
		pc=pn-dp
		return pc+dp*np.exp(-(rmax/r)**b)

	# r = radius
	# pc = Central Pressure (Pc)
	# dp = Pressure Drop=(Pn-Pc)
	# rmax = Radius of max winds (Rmax)
	# b = Holland's parameter (B)

	# ----------------------------------------
	# Velocity
	# ----------------------------------------
	# Velocity without vt, without Coriolis
	def vel(r,b,rhoa,rmax,maxR,k,dp,x):
		return (b/rhoa*(rmax/r)**b*dp*np.exp(-(rmax/r)**b))**x

	# Velocity without vt, with Coriolis
	def velcor(r,b,rhoa,rmax,maxR,k,dp,x,f):
		return (b/rhoa*(rmax/r)**b*dp*np.exp(-(rmax/r)**b)+(r*(np.absolute(f))/2)**(1/x))**x-r*(np.absolute(f))/2

	# Velocity with vt, without Coriolis
	def vel2(r,sinteta,costeta,b,rhoa,rmax,maxR,k,dp,x,vtx,vty):

		vr=vel(r,b,rhoa,rmax,maxR,k,dp,x)
		
		vx=-vr*sinteta
		vy=vr*costeta
		
		vxvt=vx+vtx*(1.-np.exp(-(rmax0/r)**b0))
		vyvt=vy+vty*(1.-np.exp(-(rmax0/r)**b0))
		
		return vxvt, vyvt

	# Velocity with vt, with Coriolis
	def vel2cor(r,sinteta,costeta,b,rhoa,rmax,maxR,k,dp,x,vtx,vty,f):
		
		vr=velcor(r,b,rhoa,rmax,maxR,k,dp,x,f)
		
		vx=-vr*sinteta
		vy=vr*costeta
		
		vxvt=vx+vtx*(1.-np.exp(-(rmax0/r)**b0))
		vyvt=vy+vty*(1.-np.exp(-(rmax0/r)**b0))
		
		return vxvt, vyvt
	
	def exportValues(x,y,vals,filename):
		print x.shape[0]
		print x.shape[1]
		out_file = open(filename,"w")

		for i in range(x.shape[0]):
			for j in range(x.shape[1]):
				out_file.write(str(x[i,j])+","+str(y[i,j])+","+str(vals[i,j])+"\n")
		out_file.close()
	#######################################################
	# Holland model
	# reading GDACS data
	#######################################################
	import numpy as np
	import pandas
	from scipy.interpolate import griddata

	print "Entering CreateHolland: " + str(minxc)+" "+str(maxxc) +" "+ str(minyc)+" "+str(maxyc) 
	#******************************************************
	# SET PARAMETERS
	#******************************************************
	nm2m=1852. # 1 nautical mile to meters
	kt2ms=nm2m/3600.  # knots to m/s
	omega=2*np.pi/(3600.*24.) # angular speed omega=2pi*f(=frequency of earth : 1 cycle per day) 2pi* 1 / day in seconds
	rhoa=1.15 #air density  Kg/m^3
	radius=6378388 #137. # earth's radius according to WGS 84
	deg2m=np.pi*radius/180.  # ds on cicle equals ds=r*dth - dth=pi/180
	pn=101000.  # Atmospheric pressure [N/m^2] (1010hPa - enviromental pressure)

	tetaNE=45. #mean angle [degrees] of North Eastern quadrant
	tetaNW=135. #        "              North Western
	tetaSW=225. #        "              South West
	tetaSE=315. #        "              South East

	maxR=500.e3  # maximum radius of TC [m] (500Km)

	# deg to rad (NEW)
	deg2rad=np.pi/180 # np.deg2rad
	deg2m=deg2rad*radius

	# set Holland initial 
	rmax0=20000
	b0=1.2

	

	# read input data
	
	print('Reading input data....')
	#namefile().split('\t')
	#totcols=14 --> outdata
	totcols=16 # --> hollanddata
	datH=pandas.read_csv(namefile,header=0,delimiter='\t',usecols=np.arange(0,totcols))

	# store input data
	timeH=datH['time'].values
	latH=datH['yhc'].values
	lonH=datH['xhc'].values
	bH=datH['b'].values
	kH=datH['k'].values
	rmaxH=datH['rmax'].values
	dpH=datH['deltap'].values
	vmaxmsH=datH['vmax'].values
	vmax0msH=datH['vmax0'].values
	vtrxH=datH['vtrx'].values
	vtryH=datH['vtry'].values


	print('------------------')

	#******************************************************
	# INTERPOLATE HOLLAND DATA
	#******************************************************
	# !!!!!!!!!!!! not finished
	print('Interpolating Holland data...')
	minint=15.
	#totH=120. #!!! last value of DatH
	print len(timeH)
	#print datH
	print timeH[len(timeH)-1]
	totH=timeH[(len(timeH)-1)]
	intH=totH*(60/minint)+1
	print 'min int: ', minint
	print 'intH: ', intH
	print 'totH: ', totH

	xvals = np.linspace(0., totH, intH)
	timeH1=xvals
	latH1=np.interp(xvals, timeH, latH)
	lonH1=np.interp(xvals, timeH, lonH)
	bH1=np.interp(xvals, timeH, bH)
	kH1=np.interp(xvals, timeH, kH)
	rmaxH1=np.interp(xvals, timeH, rmaxH)
	dpH1=np.interp(xvals, timeH, dpH)
	vmaxmsH1=np.interp(xvals, timeH, vmaxmsH)
	vmax0msH1=np.interp(xvals, timeH, vmax0msH)
	vtrxH1=np.interp(xvals, timeH, vtrxH)
	vtryH1=np.interp(xvals, timeH, vtryH)



	print ('zzGDACS/Holland time'), timeH
	#print ('New time: '), timeH1
	print('------------------') 

	#******************************************************
	# CREATE GRID
	#******************************************************

	#dcell=2./60. #(maxxc-minxc)/(2/60) # 2min
	maxxc=maxxc 
	nintx=(maxxc-minxc)/dcell+1
	ninty=(maxyc-minyc)/dcell+1
	
	print ('Creating grid...')
	print 'min lon: ', minxc
	print 'max lon: ', maxxc
	print 'min lat: ', minyc
	print 'max lat: ', maxyc
	print "dcell=",dcell
	max1=nintx*dcell+minxc
	print nintx,ninty, nintx,dcell, max1
	# create grid degrees
	#xcell0=np.linspace(minxc,maxxc-dcell,nintx-1)
	
	xcell0=np.linspace(minxc,maxxc,nintx)
	ycell0=np.linspace(minyc,maxyc,ninty)
	xcell,ycell=np.meshgrid(xcell0,ycell0)

	print 'xcell shape', xcell.shape
	print 'ycell shape', ycell.shape
	print('------------------')

	#******************************************************
	# LOOP: CREATE HOLLAND WIND FIELD
	# Calculate wind (and pressure) fields for all times
	#******************************************************
	print ('Creating WIND FIELD...LOOP....')
	print ('start....')
	done=0
	it=0
	itot=len(timeH1)
	#itot=10 #!!! test

	cvel=[]

	#cpres=[]
	for it in range(itot):
		lat=latH1[it]
		lon=lonH1[it]
		b=bH1[it]
		vmax=vmaxmsH1[it]
		vmax0=vmax0msH1[it]
		k=kH1[it]
		rmax=rmaxH1[it]
		dp=dpH1[it]
		pc=pn-dp
		vtx=vtrxH1[it]
		vty=vtryH1[it]
		#vt=vtrH1[it]

		# Define Lat and Lon of TC's center
		xhc=lon
		yhc=lat
		#print 'TC center (lat/lon): ', yhc, '/', xhc

		
		# calculate distance between TC center and point (cell)
		coscellm=np.cos(np.deg2rad(yhc-ycell))
		dycell=(ycell-yhc)*deg2m
		dxcell=(xcell-xhc)*deg2m*coscellm
		rcell=np.sqrt((dxcell)**2+(dycell)**2) # max(rcell,1000)
		rcell=np.maximum(rcell,1000)
		costeta=dxcell/rcell
		sinteta=dycell/rcell
		
		
		# Calculate x
		x=calcx(rcell,rmax,maxR,k)
		
		
		# calculate f
		f=fcor(ycell, omega)
		
		# calculate beta
		#beta=calcbeta(rcell,rmax)*deg2rad # (use np deg2rad)
		#cosb=np.cos(beta)
		#sinb=np.sin(beta)
		
		
		# calculate pressure field
		#prcell=pres(rcell, dp, pn, rmax, b)
	   
		
		# calculate v (with vt)
		
		vxvt, vyvt = vel2cor(rcell,sinteta,costeta,b,rhoa,rmax,maxR,k,dp,x,vtx,vty,f)
		vvt=np.sqrt(vxvt**2+vyvt**2)
		
		#vx, vy, vxvt, vyvt = vel2(rcell,sinteta,costeta,b,rhoa,rmax,maxR,k,dp,x,vtx,vty)
		#vx, vy, vxvt, vyvt = vel2cor(rcell,sinteta,costeta,b,rhoa,rmax,maxR,k,dp,x,vtx,vty,f)
		#v=np.sqrt(vx**2+vy**2)
		
		cvel.append(vvt)
	# print x.shape, f.shape
	# print "b,rhoa,rmax,maxR,k,dp,x,vtx,vty,f ",b,rhoa,rmax,maxR,k,dp,x,vtx,vty,f
	# print "vxvt=", vxvt[59,1],vxvt[60,1],vxvt[61,1],vxvt[62,1],vxvt[63,1]
	# print "vyvt=", vyvt[59,1],vyvt[60,1],vyvt[61,1],vyvt[62,1],vyvt[63,1]
	# print "rcell=", rcell[59,1],rcell[60,1],rcell[61,1],rcell[62,1],rcell[63,1]
	# print "dxcell=",dxcell[59,1],dxcell[60,1],dxcell[61,1],dxcell[62,1],dxcell[63,1]
	# print "dycell=",dycell[59,1],dycell[60,1],dycell[61,1],dycell[62,1]		,dycell[63,1]
	
	
	print ('....end loop')
	print('------------------')


	#******************************************************
	# CREATE MAX WIND FIELD
	#******************************************************
	print ('creating max wind field ....')
	cvel=np.array(cvel)
	#cpres=np.array(cpres)

	print 'cvel shape: ', cvel.shape
	print 'rcell shape: ', rcell.shape

	maxvel=np.amax(cvel,axis=0)
	#minpres=np.amin(cpres,axis=0)
	print('------------------')
	
	#filename='tmp_wind10m_'+hurName+"_"+str(ngrid)+'.txt'
	#print "exporting "+filename
	#exportValues(xcell,ycell,maxvel,filename)
	#return ngrid
	
	#******************************************************
	# CREATE MAP AND SAVE DATA
	#******************************************************
	# ----------------------------------------
	# create map
	# ----------------------------------------
	print ('creating map....')
	from osgeo import gdal,gdal_array
	import osr

	dataTypeformat={1:np.byte,2:np.int32,3:np.int32,4:np.float32,5:np.float32,6:np.byte}
	VSType={1:'VS_BOOLEAN',2:'VS_NOMINAL',3:'VS_ORDINAL',4:'VS_SCALAR',5:'VS_DIRECTION',6:'VS_LDD'}

	def putmap(filename,var,geo,TYPE,nodata):
	 driver=gdal.GetDriverByName('GTiff')
	 varw=var.astype(dataTypeformat[TYPE])
	 gtype=gdal_array.NumericTypeCodeToGDALTypeCode(varw.dtype)
	 NROWS,NCOLS = var.shape
	 VS='PCRASTER_VALUESCALE={}'.format(VSType[TYPE])
	 dst_ds=driver.Create(filename,NCOLS,NROWS,1,gtype,[VS])
	 proj=osr.SpatialReference()
	 proj.ImportFromEPSG(4326)
	 dst_ds.SetProjection(proj.ExportToWkt())
	 dst_ds.SetGeoTransform(geo)
	 dst_ds.GetRasterBand(1).WriteArray(varw)
	 dst_ds.GetRasterBand(1).SetNoDataValue(nodata)
	 dst_ds.FlushCache()
	 dst_ds=None
	 return

	#SAVEPATH='/mnt/pandora/Operations_Critech/Emergencies/2016/20160929_Caribbean_TC_MATTHEW/maps_7Oct/data/HWRF/'

	ih='test'
	print xcell.min(),xcell.max(),ycell.min(),ycell.max(),dcell
	TYPE=4     

	geo=(xcell.min(),dcell,0,ycell.max(),0, -dcell)  
	nodata=-9999.
	
	filenameTIF=workdir+'/tmp_wind10m_'+hurName+"_"+str(ngrid)+'.tif' #!!!!!H --> name bul_no
	print filenameTIF
	putmap(filenameTIF,np.flipud(maxvel),geo,TYPE,nodata)
	return ngrid
