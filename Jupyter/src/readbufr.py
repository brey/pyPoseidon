import numpy as np

import glob

from pybufr_ecmwf.bufr import BUFRReader
from pybufr_ecmwf.raw_bufr_file import RawBUFRFile
from pybufr_ecmwf.bufr_interface_ecmwf import BUFRInterfaceECMWF

from bunch import Bunch
import datetime


def read_bufr(stormname):
    # get a list of all bufr files in directory
    tr=glob.glob('/mnt/ECMWF/bufr/2016/*ECMF*tropical_cyclone_track_{}*'.format(stormname))
        
    trdat={} # create a dictionary for storing the data

    for name in tr:
        
        bob=BUFRReader(name,warn_about_bufr_size=False,verbose=False)
        try:
           bob.get_next_msg()
        except EOFError as e:
            print e

        nsubsets = bob.get_num_subsets()

        for subs in range(1, nsubsets+1):

            # add header strings
            (list_of_names, list_of_units) = bob.get_names_and_units(subs)

            # currently not used
            # list_of_unexp_descr = bob.bufr_obj.py_unexp_descr_list

            data = bob.get_subset_values(subs)

        hdata=zip(list_of_names,data,list_of_units)

        t=[b for a,b,c in hdata if a == 'TIME PERIOD OR DISPLACEMENT']
        t=np.array(t)
        
        p_msl=[b for a,b,c in hdata if a == 'PRESSURE REDUCED TO MEAN SEA LEVEL']
        p_msl=np.array(p_msl)

        lons=[b for a,b,c in hdata if a == 'LONGITUDE (COARSE ACCURACY)']
        lons=np.array(lons)

        lats=[b for a,b,c in hdata if a == 'LATITUDE (COARSE ACCURACY)']
        lats=np.array(lats)
        
        maxv=[b for a,b,c in hdata if a == 'WIND SPEED AT 10 M']
        maxv=np.array(maxv)


        m=np.abs(lats)<90  # screen out the nonvalues

        yyyy=[b for a,b,c in hdata if a == 'YEAR']
        mm=[b for a,b,c in hdata if a == 'MONTH']
        dd=[b for a,b,c in hdata if a == 'DAY']
        hh=[b for a,b,c in hdata if a == 'HOUR']
        ss=[b for a,b,c in hdata if a == 'MINUTE']

        #index=datetime.datetime.strptime('%04i%02i%02i%02i%02i' % (yyyy[0],mm[0],dd[0],hh[0],ss[0]),'%Y%m%d%H%S')
        index='%04i%02i%02i%02i%02i' % (yyyy[0],mm[0],dd[0],hh[0],ss[0])
        
        # set time zero
        
        t=np.append([0],t)
        
        # add to dictionary with the date data as index and the value as a bunched dictionary
        
        value=Bunch()
        value['center']=[lons[0],lats[0]]
        value['t']=t
        value['p']=p_msl
        value['plons']=lons[1::2]
        value['plats']=lats[1::2]
        value['maxv']=maxv
        value['ulons']=lons[2::2]
        value['ulats']=lats[2::2]

        
        trdat[index]=value 

        
    return trdat    


        

