#!/bin/bash
    #
    # This script is an example for running Delft3D-FLOW
    # Adapt and use it for your own purpose
    #
    # adri.mourits@deltares.nl
    # 27 Dec 2010
    # 
    #
    # This script starts a single-domain Delft3D-FLOW computation on Linux
    #


    #
    # Set the config file here
    # 
#argfile=config_flow2d3d.xml
argfile=config_d_hydro.xml





    #
    # Set the directory containing delftflow.exe here
    #
export ARCH=intel
export D3D_HOME=/home/critechuser/DELFT3D/5740/bin/lnx64
exedir=$D3D_HOME/flow2d3d/bin
libdir=$D3D_HOME/flow2d3d/bin
LD_LIBRARY_PATH=/usr/local/lib
 
    #
    # No adaptions needed below
    #

    # Set some (environment) parameters
export LD_LIBRARY_PATH=$exedir:$libdir:$LD_LIBRARY_PATH 

    # Run
mpiexec -np 4 $exedir/d_hydro.exe $argfile
