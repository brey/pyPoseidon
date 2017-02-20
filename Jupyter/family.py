import numpy as np
from scipy.interpolate import interp1d


def family(x, *args):

        farg=[]
	for arg in args:
   		farg.append(interp1d(x,arg))
        
        return farg
