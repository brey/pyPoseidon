Storm Surge estimation  software
==============================

This is a development project for utilizing DELFT3D open source code for considering storm surge. The purpose is to create, run and analyze a storm surge computation based on Python scripts and thus avoid the Matlab based GUI that deltares use. The procedure is outlined in Jupyter notebooks in order to have clarity and transparency.  

## Getting Started

A number of data are inclusive in the project folders but there are some required outside sources. The development is done in Linux environment. Some (maybe many) modifications need to be made for other platforms. 

There are some system dependencies that need to be specified. These are defined in appropriate env variables.

* Related to DELFT3D

export D3D=path/to/your/delft3d/installation. e.g. export D3D=/opt/DELFT3D/5740/bin

export ARCH=platform e.g. export ARCH=lnx64

export LD_LIBRARY_PATH=required libraries i.e. netcdf,mpich, etc.  e.g export LD_LIBRARY_PATH=/usr/local/lib

* Related to input

export ECMWF=path/to/ECMWF/data e.g. export ECMWF=/opt/ECMWF # it is assumed that in this folder the structure is year/month/day/yyyymmdd.hh/

In order to keep it compact the rest of the paths are kept relative to the /src folder. One can use symbolic links to link there the required data. See corresponding folders for more details. 

* Related to ipython

In the Notebooks pypath_magic is used to load local modules. In order to use the scipts however add the /src/ and /util/ folder to PYTHONPATH


* Related to Cython

For compiling the cython module redtoreg.pyx it might be required to set the path for the numpy headers like

export CFLAGS='-I/usr/local/lib/python2.7/site-packages/numpy/core/include:$CFLAGS'

This is done only once. 

### Prerequisities

DELFT3D needs to be compiled for your system. You can download it from http://oss.deltares.nl/web/delft3d/source-code. Follow the instruction therein.  

A number of Python modules are required. A complete list is available in a file named requirements.txt.


## Tests

No tests are available at the moment.

## Authors

* **George Breyiannis** 


## Acknowledgments

* A number of scripts have been used (some with some modifications) from https://sourceforge.net/projects/openearthtools/. 
* All the people who spend time and energy to provide solutions online.  

## License
* Since a few scripts from OpenEarthTools have been used, the project is released under the corresponding GPL v3 license. 

  This library is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

