import numpy as np
from netCDF4 import Dataset
from bilinear import bilinear_interpolation


def bl4(tideset,plon,plat):


  #read from med.nc
  dmed=Dataset(tideset)

  lat=dmed['lat'][:]
  lon=dmed['lon'][:]

  tidal_c=dmed['tidal_constituents'][:]

  tidal_c=[''.join(k).upper().strip() for k in tidal_c]

  amp=dmed['tidal_amplitude_h']
  ph=dmed['tidal_phase_h']


  i=np.abs(lon-np.float(plon)).argmin()
  j=np.abs(lat-np.float(plat)).argmin()

  x0,y0 = lon[i], lat[j]

  # retrieve the 4 nearest diagonal points of the i,j
  lon4=lon[i-1:i+2:2]
  lon4=np.vstack([lon4,lon4])
  lat4=lat[j-1:j+2:2]
  lat4=np.vstack([lat4,lat4]).T
# print '==================='
  # define the quandrant
  A=plon-x0
  B=plat-y0
  for xx,yy in zip(lon4.flatten(),lat4.flatten()):
#        print xx,yy
         C=np.sign([xx-plon,A])
         D=np.sign([yy-plat,B])
         if C[0]==C[-1] and D[0] == D[-1] : 
           corner=[xx,yy]
           indx=[np.abs(lon-xx).argmin(),np.abs(lat-yy).argmin()]

  phv=np.zeros(ph.shape[-1])
  amv=np.zeros(amp.shape[-1])
  for k in range(amp.shape[-1]):

      p1=[x0,y0,amp[i,j,k]]
      p2=[x0,corner[1],amp[i,indx[1],k]]
      p3=[corner[0],corner[1],amp[indx[0],indx[1],k]]
      p4=[corner[0],y0,amp[indx[0],j,k]]
  
      points=[p1,p2,p3,p4]

      amv[k]= bilinear_interpolation(plon,plat,points)

      p1=[x0,y0,ph[i,j,k]]
      p2=[x0,corner[1],ph[i,indx[1],k]]
      p3=[corner[0],corner[1],ph[indx[0],indx[1],k]]
      p4=[corner[0],y0,ph[indx[0],j,k]]
  
      points=[p1,p2,p3,p4]

      phv[k]= bilinear_interpolation(plon,plat,points)

  return tidal_c,amv,phv

if __name__ == "__main__":
    plon=11.3243
    plat=45.7255
    tideset='../TIDES/med.nc'
    tidal_c,amv,phv=bl4(tideset,plon,plat)


