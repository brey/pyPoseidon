import numpy as np
import datetime
import os

def meteo2delft3d(p,u,v,lat0,lon0,dlat,dlon,time,nt,path='./',curvi=False):

 nodata=-9999.000

 if not os.path.exists(path):
    os.makedirs(path)

# open files
 pfid = open(path+'p.amp','w')
 ufid = open(path+'u.amu','w')
 vfid = open(path+'v.amv','w')

 fi=[pfid,ufid,vfid]
 wi=[ufid,vfid]

# write file headers
 for f in fi:
   f.write('FileVersion      = 1.03\n')
 if curvi :
  for f in fi:
    f.write('Filetype         = meteo_on_curvilinear_grid\n')
    f.write('grid_file        = wind.grd\n')
    f.write('first_data_value = grid_ulcorner\n')
    f.write('data_row         = grid_row\n')
 else:
  for f in fi:
    f.write('Filetype         = meteo_on_equidistant_grid\n')
    f.write('n_cols           = {}\n'.format(u.shape[2]))
    f.write('n_rows           = {}\n'.format(u.shape[1]))
    f.write('grid_unit        = degree\n')
    # code currently assumes lon and lat are increasing
    f.write('x_llcenter       = {:g}\n'.format(lon0))
    f.write('dx               = {:g}\n'.format(dlon))
    f.write('y_llcenter       = {:g}\n'.format(lat0))
    f.write('dy               = {:g}\n'.format(dlat))

 for f in fi:
  f.write('NODATA_value     = {:.3f}\n'.format(nodata))
  f.write('n_quantity       = 1\n')

 ufid.write('quantity1        = x_wind\n')
 vfid.write('quantity1        = y_wind\n')
 pfid.write('quantity1        = air_pressure\n')

 for f in wi:
  f.write('unit1            = m s-1\n')

 pfid.write('unit1            = Pa\n')

 time0=datetime.datetime.strptime('2016-01-01 00:00:00','%Y-%m-%d %H:%M:%S')

# write time blocks
 for it in range(nt+1): # nt + 0 hour
 
   ntime=time+datetime.timedelta(hours=it)
   dt=(ntime-time0).total_seconds()/3600.


   for f in fi:
     f.write('TIME = {} hours since 2016-01-01 00:00:00 +00:00\n'.format(dt))

   np.savetxt(pfid,np.flipud(p[it,:,:]/0.01),fmt='%.3f')
   np.savetxt(ufid,np.flipud(u[it,:,:]),fmt='%.3f')
   np.savetxt(vfid,np.flipud(v[it,:,:]),fmt='%.3f')

# write the same values for the end time
#dt=dt+nt*60.
#for f in fi:
# f.write('TIME = {} hours since 1900-01-01 00:00:00 +00:00\n'.format(dt))

#np.savetxt(pfid,np.flipud(p/0.01),fmt='%.3f')
#np.savetxt(ufid,np.flipud(u),fmt='%.3f')
#np.savetxt(vfid,np.flipud(v),fmt='%.3f')

# close files
 for f in fi:
   f.close()
