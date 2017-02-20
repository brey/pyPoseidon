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
#from modnc import fTAT

import subprocess
import logging


SAVEPATH='/mnt/ECMWF/processed/2016/FIX_MED_SEA_D3/'
RUNPATH='/home/critechuser/'


def go(rundate,path,dic,TAT=False):

   bname=dic['bname']
   nt=dic['nt']
   ni=dic['ni']
   nj=dic['nj']
   lon0=dic['lon0']
   lat0=dic['lat0']
   lon1=dic['lon1']
   lat1=dic['lat1']

   files=['config_d_hydro.xml',bname+'.mdf',bname+'.grd',bname+'.enc',bname+'.obs',bname+'.dep', bname+'.pkl','run_flow2d3d.sh']

   logging.info(rundate)
   sys.stdout.write('RUNNING {}\n'.format(rundate))
   sys.stdout.flush()

# previous date for reference

   previous_run=rundate-datetime.timedelta(hours=12)

# create the folder/run path

   folder=datetime.datetime.strftime(rundate,'%Y%m%d.%H' )

   rpath=path+'{}/'.format(folder)   

   if not os.path.exists(rpath):
    os.makedirs(rpath)

# copy necessary files
   pfolder=datetime.datetime.strftime(previous_run,'%Y%m%d.%H' )
   ppath =  path+'{}/'.format(pfolder)  # previous path
   
   logging.info(ppath)

   for filename in files:

      copy2(ppath+filename,rpath+filename)

# copy restart file

   inresfile='tri-rst.'+bname+'.'+datetime.datetime.strftime(rundate,'%Y%m%d.%H%M%M')

   outresfile='restart.'+datetime.datetime.strftime(rundate,'%Y%m%d.%H%M%M')

   copy2(ppath+inresfile,rpath+'tri-rst.'+outresfile)

#get new meteo 

   sys.stdout.write('process meteo\n')
   sys.stdout.flush()

   check=[os.path.exists(rpath+f) for f in ['u.amu','v.amv','p.amp']]   
   if np.any(check)==False :

       p,u,v,lat,lon = wmap(rundate,0,3*(nt+1),lon0,lon1,lat0,lat1)


       sys.stdout.write('\n')

#write u,v,p files 

       dlat=lat[1,0]-lat[0,0]
       dlon=lon[0,1]-lon[0,0]

       sys.stdout.write('create delft3d files\n')
       sys.stdout.flush()

       meteo2delft3d(p,u,v,lat0,lon0,dlat,dlon,rundate,nt,path=rpath,curvi=False)

   else:
       sys.stdout.write('meteo files present\n')

   sys.stdout.write('\n')

# modify mdf file

   inp, ord = mdf.read(rpath+bname+'.mdf')

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
   mdf.write(inp, rpath+bname+'.mdf',selection=ord)
# run case

#  p=u=v=lon=lat=None

   sys.stdout.write('start run\n')
   sys.stdout.flush()

 
   os.chdir(rpath)
  #subprocess.call(rpath+'run_flow2d3d.sh',shell=True)
   os.system('./run_flow2d3d.sh')

   sys.stdout.write('\n')


   if TAT :

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



