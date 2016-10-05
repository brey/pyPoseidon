import numpy as np
import datetime
import sys
import os
from shutil import copy2

#from getd import wmap
#from getmedf import wmap
from meteo import wmap
from idelft3d import meteo2delft3d
import mdf
from modnc import fTAT

import subprocess
import logging

#Grid size (fixed)

ni=727
nj=285
lon0=-5.5
lat0=28.5
lon1=43.
lat1=47.5

path='/mnt/web/brey/2016H/'
SAVEPATH='/mnt/ECMWF/processed/2016/FIX_MED_SEA_D3/'
RUNPATH='/home/critechuser/DELFT3D/python/'

nt=72

files=['config_d_hydro.xml','med.mdf','med.grd','med.enc','med.obs','med.dep', 'med.pkl','run_flow2d3d.sh']

def go(rundate):
   logging.info(rundate)
   sys.stdout.write('RUNNING {}\n'.format(rundate))
   sys.stdout.flush()

# previous date for reference

   previous_day=rundate-datetime.timedelta(days=1)

# create the folder/run

   month_dir=path+'{}/'.format(rundate.month)  #month subfolder  

   if not os.path.exists(month_dir):
    os.makedirs(month_dir)

   day_dir=month_dir+'{}/'.format(rundate.day)  #day subfolder  

   folder=day_dir+'{:02d}/'.format(rundate.hour)  #daily run subfolder :  0 or 12

   if not os.path.exists(folder):
    os.makedirs(folder)


# copy necessary files
   ppath = day_dir+'00/' if rundate.hour == 12 else path+'{}/{}/{}/'.format(previous_day.month,previous_day.day,12)  # previous path
   
   logging.info(ppath)

   for filename in files:

      copy2(ppath+filename,folder+filename)

# copy restart file

   inresfile='tri-rst.med.'+datetime.datetime.strftime(rundate,'%Y%m%d.%H%M%M')

   outresfile='restart.'+datetime.datetime.strftime(rundate,'%Y%m%d.%H%M%M')

   copy2(ppath+inresfile,folder+'tri-rst.'+outresfile)

#get new meteo 

   yyyy=rundate.year
   mm=rundate.month
   dd=rundate.day
   hh=rundate.hour

   sys.stdout.write('process meteo\n')
   sys.stdout.flush()


 # p,u,v,lon,lat,bat = wmap(yyyy,mm,dd,hh,0,nt,lon0,lon1,lat0,lat1,ni,nj,save=False)
 # p,u,v,lon,lat = wmap(yyyy,mm,dd,hh,0,3*(nt+1),lon0,lon1,lat0,lat1,ni,nj)
   p,u,v,lat,lon = wmap(yyyy,mm,dd,hh,0,3*(nt+1),lon0,lon1,lat0,lat1)


   sys.stdout.write('\n')

#write u,v,p files 

   dlat=lat[1,0]-lat[0,0]
   dlon=lon[0,1]-lon[0,0]

   sys.stdout.write('create delft3d files\n')
   sys.stdout.flush()

   meteo2delft3d(p,u,v,lat0,lon0,dlat,dlon,rundate,nt,path=folder,curvi=False)

   sys.stdout.write('\n')

# modify mdf file

   inp, ord = mdf.read(folder+'med.mdf')

  # adjust iteration date
   tstart=rundate.hour*60
   inp['Itdate']=datetime.datetime.strftime(rundate,'%Y-%m-%d')
   inp['Tstart']=[tstart]
   inp['Tstop']=[4320+tstart]
   inp['Flmap']  = [tstart,60,4320+tstart]
   inp['Flhis']  = [tstart,1,4320+tstart]


   if 'Restid' not in ord : ord.append('Restid')
  # adjust restart file   
   inp['Restid']=outresfile

  # update mdf
   mdf.write(inp, folder+'med.mdf',selection=ord)
# run case

#  p=u=v=lon=lat=None

   sys.stdout.write('start run\n')
   sys.stdout.flush()

 
   os.chdir(folder)
  #subprocess.call(folder+'run_flow2d3d.sh',shell=True)
   os.system('./run_flow2d3d.sh')

   sys.stdout.write('\n')


   sys.stdout.write('start analysis for TAT\n')
   sys.stdout.flush()

  # TAT interface
   try:
     fTAT()
   except Exception as e:
     print e
     sys.exit()

   sys.stdout.write('\n')
   sys.stdout.write('analysis finished\n')
   sys.stdout.flush()
   sys.stdout.write('\n')

   fname='calc_{}'.format(datetime.datetime.strftime(rundate,'%Y%m%d.%H'))
   now=datetime.datetime.strftime(datetime.datetime.now(),'%d %b %d %H:%M:%S CET %Y')
   with open(SAVEPATH+'{}/completed.txt'.format(fname),'w') as f:
       f.write(now)
   with open(SAVEPATH+'final/completed.txt'.format(fname),'w') as f:
       f.write(now)


   with open(RUNPATH+'logMED.txt'.format(fname),'a') as f:
       f.write('{}\n'.format(datetime.datetime.strftime(rundate,'%Y%m%d.%H')))
   
   sys.stdout.write('completed\n')
   sys.stdout.flush()
   sys.stdout.write('\n')


if __name__ == "__main__":
    inp=sys.argv[1]
    rdate=datetime.datetime.strptime(inp,'%Y%m%d.%H')
    go(rdate)



