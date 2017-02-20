#######################################################
# Create wind filed + classify pop
#######################################################
# 	print "Syntax: numb_intervals cellsize_min HurName GDACS_ID year"

#******************************************************
# IMPORT Python lib 
#******************************************************
import pandas
import numpy as np
from scipy.interpolate import griddata
import sys
import os
import pp
import time
import shutil
import datetime
from HollandH import create
from classifyH import classFile



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

#totcols=14 --> outdata
totcols=16 # --> hollanddata

ddeg=2 # addition of 2 deg to borders




#******************************************************
# SET INPUT/OUTPUT DIR/FILES
#******************************************************
GDACShome='/mnt/web/cycloneSurgeVM/'
basedir='.'
#basedir="/home/critechuser/cycloneSurge/Holland"
print 'Home dir'
print os.getcwd()
os.chdir (basedir)
print os.getcwd()


now = datetime.datetime.now()
#log=open("logfile.txt","a")
#log.write(str(now) + ": create_grid_ppH.py " +str(sys.argv)+"\n" )
#log.close()

print ("-------------------------------------------------------------------------------")
print('>>INPUT ARGUMENTS')
print ("-------------------------------------------------------------------------------")
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

if len(sys.argv)<5:
	print "Syntax: numb_intervals cellsize_min HurName GDACS_ID bul_no year "
	quit()
	
nint=int(sys.argv[1])
cellSize0=float(sys.argv[2])
cellSize=cellSize0/60
hurName=sys.argv[3]
GDACS_ID=sys.argv[4]
bul=int(sys.argv[5])
yearHur=int(sys.argv[6])

print 'numb_intervals: ' + str(nint )
print 'cell size: '+ str(cellSize)
print 'Hurricane name: ' + hurName 
print 'GDACS ID: ' + GDACS_ID 
print 'bulletin number: '+ str(bul)
print ("-------------------------------------------------------------------------------")


GDACSpath=GDACShome+GDACS_ID

if bul == 0: # da sistemare
	bul='final'
	GDACSdir=GDACSpath+'/'+bul
	namefile=GDACSdir+'/hollandData.txt'
	workdir="./"+str(yearHur)+"/"+hurName+"/"+str(cellSize0)+"min/"+bul
else:
	GDACSdir=GDACSpath+'/{}'.format(bul)
	namefile=GDACSdir+'/outData.txt'
	workdir="./"+str(yearHur)+"/"+hurName+"/"+str(cellSize0)+"min/"+str(bul)


if os.path.exists(workdir):
	shutil.rmtree(workdir)
 
os.makedirs(workdir)

# check if hollandData file exists in GDACSdir
if not os.path.exists(namefile):
	print "hollandData.txt file  does not exists: " + namefile
	exit()
	
print ("-------------------------------------------------------------------------------")
print('>>INPUT DATA')
print ("-------------------------------------------------------------------------------")
print 'Name:'+hurName
print 'Input GDACS Dir: '+GDACSpath
print 'Bulletin nr:'+(str(bul))
print 'namefile: '+ namefile+'\n'
print 'workdir: '+ workdir+'\n'
print 'Base dir:' + basedir




#******************************************************
# Reading input file
#******************************************************
# read input data
print ("-------------------------------------------------------------------------------")
print('>> Reading input data: '+ namefile)
print ("-------------------------------------------------------------------------------")
datH=pandas.read_csv(namefile,header=0,delimiter='\t',usecols=np.arange(0,totcols))

# store input data
timeH=datH['time'].values
latH=datH['yhc'].values
lonH=datH['xhc'].values


nrow=0 #float(nint/2) #!!!! da sistemare
ncol=nint #float(nint/2) #!!!! da sistemare
print 'nrow', nrow
print 'ncol', ncol

min_xhc=min(lonH)
max_xhc=max(lonH)
min_yhc=min(latH)
max_yhc=max(latH)


minx=round(min_xhc-ddeg,1)
maxx=round(max_xhc+ddeg,1)
miny=round(min_yhc-ddeg,1)
maxy=round(max_yhc+ddeg,1)

dxtot=round((maxx-minx),1)
dytot=round((maxy-miny),1)

if ncol==0:
  dx=maxx-minx
else:
  dx=round((maxx-minx)/ncol,2)

nn=int(dx/cellSize)
dx=(nn+1)*cellSize


dy=round((maxy-miny),1)  


print '--------------'
print '> GRID parameters'
print '--------------'
print 'min lon grid: ', minx
print 'max lon grid: ', maxx
print 'min lat grid: ', miny
print 'max lat grid: ', maxy
print 'dx tot: ', dxtot
print 'dy tot: ', dytot
print 'ninterval: ', nint
print 'dx: ', dx
print 'dy: ', dy


#minx0=minx
#miny0=miny
#minxi=minx0
#minyi=miny0


print ("-------------------------------------------------------------------------------")
print('>> Create jobserver + run hollandH.py')
print ("-------------------------------------------------------------------------------")

job_server = pp.Server()
job_server.set_ncpus(nint)
jobs = []
start_time=time.time()
print "Starting ", job_server.get_ncpus(), " workers"

print 'Starting HOLLAND' 

minxi=minx
minyi=miny
maxyi=maxy
ng=0 
for j in range(int(ncol)):
	ng=ng+1
	maxxi=minxi+dx
	print minxi,maxxi,minyi,maxyi
	print ("launching job "+str(ng))
	jobs.append(job_server.submit(create, (minxi,maxxi,minyi,maxyi,ng,hurName,namefile,workdir,cellSize)))
	minxi=maxxi 

  
part_sum1 = sum([job() for job in jobs])

# Print the partial sum
print "Partial sum is", part_sum1

print "Time elapsed: ", time.time() - start_time, "s"
job_server.print_stats()


print ("-------------------------------------------------------------------------------")
print('>> Merge files')
print ("-------------------------------------------------------------------------------")
# #outfile="/vmax_"+hurName+"_"+'{}'.format(bul)+"_"+'{}'.format(cellSize0)+"min_final.tif"

if bul == 'final':
	outfile="/vmax_"+hurName+"_"+bul+'{}'.format(cellSize0)+"min.tif"

else:
	outfile="/vmax_"+hurName+"_"+'{}'.format(bul)+"_"+'{}'.format(cellSize0)+"min.tif"

print outfile

var="time python gdal_merge.py "+workdir+"/tmp*.tif -o "+workdir+outfile
print var
os.system(var) 

#-----
#print 'copy to GDACS'
# var="cp "+workdir+outfile+" "+ GDACSdir
# print var
# os.system (var)
#-----



for j in range(int(ncol)):
	tmpFilename=workdir+'/tmp_wind10m_'+hurName+"_"+str(j+1)+'.tif'
	#os.remove(tmpFilename)
#print 'tmp tif files deleted'


print ("-------------------------------------------------------------------------------")
print('>> classFile')
print ("-------------------------------------------------------------------------------")
var2=workdir+"/vmax_"+hurName+"_"+'{}'.format(cellSize0)+"min.tif" #!!!
var3=workdir+"/ClassPop_"+hurName+"_"+'{}'.format(cellSize0)+"min.xml" #!!!
var4=GDACSdir+"/ClassPop_"+hurName+"_"+'{}'.format(cellSize0)+"min.xml" #!!!

print var2
print var3
print var4
if bul == 'final':
	desc="Cyclone: "+hurName+" BullNo: "+bul+" cellSize="+str(cellSize0)
else:
	desc="Cyclone: "+hurName+" BullNo: "+str(bul)+" cellSize="+str(cellSize0)

# classFile (var2,var3,workdir,desc)
# print'classFile completed'


# # print ("-------------------------------------------------------------------------------")
# # print('>> Copying xml file to GDACS dir')
# # print ("-------------------------------------------------------------------------------")

# # print "Copying xml file to GDACS dir:  shutil.copyfile("+var3+","+var4+")\n "
# # shutil.copyfile(var3,var4)

# # os.system("ls -la "+GDACSdir+"/vmax*")
# # os.system("ls -la "+GDACSdir+"/ClassPop_*")

# # print "Removing the files:"
# # rmfile1=workdir+"/popfile_clipped.tif"
# # rmfile2=workdir+"/u10res.tif"
# # print rmfile1
# # print rmfile2
# # #os.remove(rmfile1)
# # #os.remove(rmfile2)

# print "\n\n **** Job completed  **** :) **** " 






