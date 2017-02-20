import numpy as np
import os
import xml.etree.ElementTree as et
from bunch import bunchify
from datetime import datetime 

def xmlr(fname):

      tree=et.parse(fname)
      root=tree.getroot()

      vars=['.//advisory_datetime', './/longitude', './/latitude', './/wind_speed', './/windrad_nm_64kt_ne' , './/windrad_nm_64kt_se', './/windrad_nm_64kt_sw', './/windrad_nm_64kt_nw', './/windrad_nm_50kt_ne' , './/windrad_nm_50kt_se', './/windrad_nm_50kt_sw', './/windrad_nm_50kt_nw', './/windrad_nm_34kt_ne' , './/windrad_nm_34kt_se', './/windrad_nm_34kt_sw', './/windrad_nm_34kt_nw']
      qv=['time' ,  'lon'  ,  'lat'  , 'vmax' ,  'ne64'  , 'se64' ,  'sw64' ,  'nw64'  , 'ne50' ,  'se50' ,  'sw50'  , 'nw50' ,  'ne34' ,  'se34' ,  'sw34'  , 'nw34' ]

      q={k: [] for k in qv}
      
      for x, y in zip(vars,qv):
        for i in root.findall(x):
	  q[y]=np.append(q[y],i.text) 

      for key in q:
	if key != 'time' : 
                q[key]=q[key].astype(float)
        else:
                t0=datetime.strptime(q[key][0], '%d %b %Y %H:%M')#:%S %p')
                q[key][0]=0
                for l in range(1,q[key].size):
                   t=datetime.strptime(q[key][l], '%d %b %Y %H:%M')#:%S %p') 
                   dt=t-t0
                   q[key][l]=np.float((dt.total_seconds()/3600))

      q['time']=q['time'].astype(float)
      return bunchify(q)

if __name__ == "__main__":
    inp='/home/brey/cycloneSurge/storage/public/1000147/input/1/inpData.xml'
    xmlr(inp)
