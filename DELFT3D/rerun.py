import numpy as np
import datetime
import sys
from workflow import *

def comp(sd,fd):

    
  dt=(fd-sd).total_seconds()
  ndt=dt/(3600*12)
  ndt=np.int(ndt)+1

  for it in range(ndt): 
    idate=sd+datetime.timedelta(hours=12*it)
    print datetime.datetime.strftime(idate,'%Y%m%d.%H')
    go(idate)

  


if __name__ == "__main__":
    inp1=sys.argv[1]
    sdate=datetime.datetime.strptime(inp1,'%Y%m%d.%H')
    inp2=sys.argv[2]
    fdate=datetime.datetime.strptime(inp2,'%Y%m%d.%H')
    comp(sdate,fdate)

