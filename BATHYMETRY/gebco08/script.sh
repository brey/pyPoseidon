for idx in `ls *.nc`; do # Bourne Shell
  ncpdq -a lon,lat,depth $idx foo_${idx}  # Make x record dimension
done
ncrcat foo_*.nc out.nc       # Concatenate along x
#ncpdq -a time,x out.nc out.nc # Revert to time as record dimension

