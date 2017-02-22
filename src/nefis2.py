# -*- coding: utf-8 -*-
"""
Neutral File System (NEFIS) Python Wrapper

@authors: julio werner <juliowerner@ufpr.br>; 
"""

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import ctypes as ct
from sys import byteorder
import numpy as np

## LOAD NEFIS LIBRARY
"""
LINUX:
libNefisS0.so can be found at 
src/utils_lgpl/nefis/packages/nefis/src_so/.libs/
after compiling.

WINDOWS:
nefis.dll
\src\utils_lgpl\nefis\lib\Release\dynamic
after compiling the release
"""

#libnefispath = "../libs/libnefis.so"
libnefispath = "../libs/libNefisSO.so.0.0.0"
libnefis = ct.CDLL(libnefispath) # LOAD NEFIS C LIBRARY

## SOME CONSTANTS
MAX_CEL_DIM        = 100 # Maximum number of elements (extracted from nef-def.h)
MAX_DESC           = 64
MAX_DIM            = 5
MAX_VARgroups      = 10
MAX_NAME           = 16
SIZE_ERROR_MESSAGE = 1024

##  TYPES DEFINITIONS
pint      = ct.POINTER(ct.c_int)
pfloat    = ct.POINTER(ct.c_float)
pchar     = ct.c_char_p
aint3     = ct.ARRAY(ct.c_int, 3)
aint5     = ct.ARRAY(ct.c_int, 5)
aint35    = ct.ARRAY(aint3, 5)
achartype = ct.ARRAY(ct.ARRAY(ct.c_char, MAX_NAME+1), MAX_CEL_DIM)

## FUNCTIONS ARGTYPES
libnefis.Crenef.argtypes = [pint, pchar, pchar, ct.c_char, ct.c_char]
libnefis.Inqfel.argtypes = libnefis.Inqnel.argtypes = [pint,
    pchar, pchar, pchar, pchar, pchar,
    pint, pint, pint, pint]
libnefis.Inqfcl.argtypes = libnefis.Inqncl.argtypes = [pint,
    pchar, pint, pint, achartype]
libnefis.Inqncl.argtypes = libnefis.Inqncl.argtypes = [pint,
    pchar, pint, pint, achartype]
libnefis.Inqfgr.argtypes = libnefis.Inqngr.argtypes = [pint,
    pchar, pchar, pint, pint, pint]
libnefis.Inqcel.argtypes = [pint, pchar,  pint, achartype]
libnefis.Inqelm.argtypes = [pint, pchar, pchar, pint, pchar, pchar,
    pchar, pint, pint]
libnefis.Inqmxi.argtypes = [pint, pchar, pint]

class Group:
    def __init__(self, grpname, cellname, ndm, dms, grpord):
        self.name     = grpname.value
        self.cellname = cellname.value
        self.ndm      = ndm.value  # number of dimensions
        self.dms      = list(dms)  # dimensions
        self.ord      = list(grpord)  # fast order to read
        self.nelems   = int()
        self.elements = list()
        self.maxi     = int()
        
    def __str__(self):
        output = "Group Name: "           + self.name      + "\n"
        output += "Cell Name: "           + self.cellname  + "\n"
        output += "N. of Dimensions: "    + str(self.ndm)  + "\n"
        output += "Max Number of index: " + str(self.maxi) + "\n"
        output += "Group shape: "         + str(self.dms)  + "\n"
        output += "Faster Order: "        + str(self.ord)
        return output
        
  
class Element:
    def __init__(self, name, tp, nbyt, desc, ndim, shape, elbytes):
        self.name   = name.strip()
        self.type   = tp.strip()
        self.nbytes = nbyt
        self.desc   = desc
        self.ndim   = ndim
        self.shape  = shape
        self.bytes  = elbytes

    def __str__(self):
        output = "Element Name: "           + self.name        + "\n"
        output += "Description: "           + str(self.desc)   + "\n"
        output += "Element Type: "          + self.type        + "\n"
        output += "Bytes per data: "        + str(self.bytes)  + "\n"
        output += "Total Amount of Bytes: " + str(self.nbytes) + "\n"
        output += "N. of Dimensions: "      + str(self.ndim)   + "\n"
        output += "Shape: "                 + str(self.shape)  + "\n"
        return output
        
class Nefis:
    """
    Neutral File System (NEFIS)
    """
    def __init__(self, datpath, defpath, access='r'):
        self._datpath     = datpath
        self._defpath     = defpath
        self._access      = access  # 'r'=readonly; 'c' = create, 'u' = update
        self.groups       = dict()
        self._fid         = ct.byref(ct.c_int())  # file identification
        self._byteorder   = byteorder[0]  # 'L'=Little endian, 'B'=Big endian
        self._allelements = dict()
        self._update()
        
    def _update(self):
        dat_file = ct.create_string_buffer(self._datpath)
        def_file = ct.create_string_buffer(self._defpath)
        coding   = ct.c_char(self._byteorder)
        rdwr     = ct.c_char(self._access)   
        error    = libnefis.Crenef(self._fid, dat_file, def_file, coding, rdwr)
        if error:
           libnefis.Neferr(error, ct.c_char_p)
        self._updategroupinfo()

    def _updategroupinfo(self):
        error = False
        first_time = True
        while not error:
            grpnam = ct.create_string_buffer(MAX_NAME)
            celnam = ct.create_string_buffer(MAX_NAME)
            grpndm = ct.c_int(5)
            grpdms = (ct.c_int*5)()
            grpord = (ct.c_int*5)()
            if first_time:
                error = libnefis.Inqfgr(self._fid, grpnam, celnam,
                                        grpndm, grpdms, grpord)
                first_time = False
            else:
                error = libnefis.Inqngr(self._fid, grpnam, celnam,
                                        grpndm, grpdms, grpord)
            if error: # Must to improve error handle.
                break
            else:
                tmp = Group(grpnam, celnam, grpndm, grpdms, grpord)
                maxindex = ct.c_int()
                libnefis.Inqmxi(self._fid, grpnam, maxindex)
                tmp.maxi = maxindex.value
                self.groups[grpnam.value.strip()] = tmp 
        self._getallelement()
        for group in self.groups.values():
            nelems = ct.c_int(MAX_CEL_DIM)
            celnam = ct.create_string_buffer(group.cellname, 16)
            elmnms = ((ct.c_char * (MAX_NAME+1)) * MAX_CEL_DIM)()
            error = libnefis.Inqcel(self._fid, celnam, nelems, elmnms)
            group.nelems = nelems.value
            for name in elmnms[:group.nelems]:
                group.elements.append(self._allelements[name.value.strip()])
                
    def _getallelement(self):
        error = False
        first_time = True
        while not error:
            elname = ct.create_string_buffer(MAX_NAME)
            elmtyp = ct.create_string_buffer(8)
            nbytes = ct.c_int()
            nbytsg = ct.c_int()
            elmqty = ct.create_string_buffer(MAX_NAME)
            elmunt = ct.create_string_buffer(MAX_NAME)
            elmdes = ct.create_string_buffer(MAX_DESC)
            elmndm = ct.c_int()
            elmdms = (ct.c_int*5)()
            if first_time:
                error = libnefis.Inqfel(self._fid, elname, elmtyp, elmqty,
                                        elmunt, elmdes, nbytes, nbytsg,
                                        elmndm, elmdms)        
                first_time = False
            else:
                error = libnefis.Inqnel(self._fid, elname, elmtyp, elmqty,
                                        elmunt, elmdes, nbytes, nbytsg,
                                        elmndm, elmdms)        
            if error:
                break
            self._allelements[elname.value.strip()] = Element(elname.value,
                elmtyp.value, nbytsg.value, elmdes.value, elmndm.value,
                list(elmdms)[:elmndm.value], nbytes.value)

    def getdata(self, groupname, elementname, grprange = ''):    
        grpnam  = ct.create_string_buffer(groupname  , MAX_NAME)
        elmnam  = ct.create_string_buffer(elementname, MAX_NAME)     
        group   = self.groups[groupname]
        element = self._allelements[elementname]
        shape   = tuple(element.shape)[::-1]
        if element.type.lower() == 'real':
            dt = np.dtype('f4')
            if element.ndim == 1:
                narray = element.shape[0]
                if narray == 0:
                    libnefis.Getelt.argtypes = [pint, pchar, pchar, aint35,
                                                aint5, pint, pfloat]
                    dtbuffer = ct.c_float()
                elif narray > 0:
                    typec = np.ctypeslib.ndpointer(dtype=dt, ndim=1,
                                                   shape=shape)
                    libnefis.Getelt.argtypes = [pint, pchar, pchar, aint35,
                                               aint5, pint, typec]
                    dtbuffer = np.empty(shape, dtype=dt)
            elif element.ndim >= 2:
                typec = np.ctypeslib.ndpointer(dtype=dt,
                                               ndim=element.ndim,
                                               shape=shape)
                dtbuffer = np.empty(shape, dtype=dt)
                libnefis.Getelt.argtypes = [pint, pchar, pchar, aint35,
                                               aint5, pint, typec]
        elif element.type.lower() == 'integer':
            dt = np.dtype('i4')
            if element.ndim == 1:
                narray = element.shape[0]
                if narray == 0:
                    libnefis.Getelt.argtypes = [pint, pchar, pchar, aint35,
                                                aint5, pint, pint]
                    dtbuffer = ct.c_int()
                elif narray > 0:
                    typec = np.ctypeslib.ndpointer(dtype=dt, ndim=1,
                                                   shape=shape)
                    libnefis.Getelt.argtypes = [pint, pchar, pchar, aint35,
                                               aint5, pint, typec]
                    dtbuffer = np.empty(shape, dtype=dt)
            elif element.ndim >= 2:
                typec = np.ctypeslib.ndpointer(dtype=dt,
                                               ndim=element.ndim,
                                               shape=shape)
                dtbuffer = np.empty(shape, dtype=dt)
                libnefis.Getelt.argtypes = [pint, pchar, pchar, aint35,
                                               aint5, pint, typec]
        output = list()
        if grprange == '':
            grprange = range(group.maxi)
        for it in grprange:
            # Get all data
            x = aint3(1+it, 1+it, 1)
            y = aint3(0, 0, 0)
            z = aint3(0, 0, 0)
            t = aint3(0, 0, 0) 
            w = aint3(0, 0, 0) 
            uindex = aint35(x, y, z, t, w)
            usrord = aint5(1, 2, 3 , 4 , 5)  # dummy parameter
            buflen = ct.c_int(element.nbytes)
            error = libnefis.Getelt(self._fid, grpnam, elmnam, uindex,
                                    usrord, buflen, dtbuffer)
            if error:
                error_string = ct.create_string_buffer(SIZE_ERROR_MESSAGE)
                libnefis.Neferr(error, error_string)
            output.append(dtbuffer.copy())
        if  len(output) == 1:
            return dtbuffer
        else:
            return output
