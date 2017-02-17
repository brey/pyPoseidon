import numpy as np
import datetime
import sys
import os
import logging

import mdf
from modnc import fTAT

path='/mnt/web/brey/2016H/'
SAVEPATH='/mnt/ECMWF/processed/2016/FIX_MED_SEA_D3/'
RUNPATH='/home/critechuser/DELFT3D/python/'

def go(rundate):

  # TAT interface
   try:
     fTAT(datetime.datetime.strftime(rundate,'%Y%m%d.%H'))
   except Exception as e:
     print e
     sys.exit()

   fname='calc_{}'.format(datetime.datetime.strftime(rundate,'%Y%m%d.%H'))
   now=datetime.datetime.strftime(datetime.datetime.now(),'%d %b %d %H:%M:%S CET %Y')
   with open(SAVEPATH+'{}/completed.txt'.format(fname),'w') as f:
       f.write(now)
   with open(SAVEPATH+'final/completed.txt'.format(fname),'w') as f:
       f.write(now)


   with open(RUNPATH+'logMED.txt'.format(fname),'a') as f:
       f.write('{}\n'.format(datetime.datetime.strftime(rundate,'%Y%m%d.%H')))
   

if __name__ == "__main__":
    inp=sys.argv[1]
    rdate=datetime.datetime.strptime(inp,'%Y%m%d.%H')
    go(rdate)



