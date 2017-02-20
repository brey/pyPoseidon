import numpy as np
import datetime
import sys
from workflowTOMAS import *
import logging
import pickle


def comp(sd,fd,path,TAT):

  case=path.split('/')[-2]

  with open(path+case+'.pkl', 'r') as f:
    att=pickle.load(f)


  TAT=False

  logging.basicConfig(filename=path+case+'.log',level=logging.INFO)
    
  dt=(fd-sd).total_seconds()
  ndt=dt/(3600*12)
  ndt=np.int(ndt)+1

  for it in range(ndt): 
    idate=sd+datetime.timedelta(hours=12*it)
    logging.info(datetime.datetime.strftime(idate,'%Y%m%d.%H'))
    go(idate,path,att,TAT)

  


if __name__ == "__main__":
    inp1=sys.argv[1]
    sdate=datetime.datetime.strptime(inp1,'%Y%m%d.%H')
    inp2=sys.argv[2]
    fdate=datetime.datetime.strptime(inp2,'%Y%m%d.%H')
    path=sys.argv[3]
    TAT=sys.argv[4]
    comp(sdate,fdate,path,TAT)

