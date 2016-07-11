import numpy as np
from osgeo import gdal,gdal_array
import osr
from bunch import bunchify

dataType={1:'Boolean',2:'Nominal',3:'Ordinal',4:'Scalar',5:'Directional',6:'Ldd'}
dataTypeformat={1:np.byte,2:np.int32,3:np.int32,4:np.float32,5:np.float32,6:np.byte}
nodataType={'Ldd':255,'Boolean':255,'Nominal':-2147483647,'Ordinal':-2147483647,'Scalar':np.nan,'Directional':np.nan}
VSType={1:'VS_BOOLEAN',2:'VS_NOMINAL',3:'VS_ORDINAL',4:'VS_SCALAR',5:'VS_DIRECTION',6:'VS_LDD'}


def viewmap(filename):
     figure()

     m=gdal.Open(filename)
     data=m.ReadAsArray()
     mp=np.ma.where(data>-1.e-30) 
     dmin, dmax = data[mp].min(), data[mp].max()
     im=imshow(data, interpolation='nearest', vmin=dmin, vmax=dmax)
     im.cmap.set_under('grey')
     title(np.str(filename))
     colorbar()
     draw()

     show(block=False)



def getmap(filename):
     dic={}
     dataset=gdal.Open(filename)
     dic['DRIVER']= dataset.GetDriver().LongName
     dic['NCOLS']=dataset.RasterXSize
     dic['NROWS']=dataset.RasterYSize 
     dic['GeoTr']=dataset.GetGeoTransform()
     dic['Proj']=dataset.GetProjection()
     geotransform = dataset.GetGeoTransform()
     if not geotransform is None:
        dic['XSTART']=geotransform[0]
        dic['YSTART']=geotransform[3]
        dic['CELLSIZE']=geotransform[1]

     v = dataset.GetRasterBand(1)
     dic['data']=v.ReadAsArray().astype(np.float32)
     dic['nan']=v.GetNoDataValue()

     dataset=None

     return bunchify(dic)


def putmap(filename,var,geo,TYPE):
     driver=gdal.GetDriverByName('PCRaster')
     varw=var.astype(dataTypeformat[TYPE])
     gtype=gdal_array.NumericTypeCodeToGDALTypeCode(varw.dtype)
     NROWS,NCOLS = var.shape
     VS='PCRASTER_VALUESCALE={}'.format(VSType[TYPE])
     dst_ds=driver.Create(filename,NCOLS,NROWS,1,gtype,[VS])
     proj=osr.SpatialReference()
     proj.SetWellKnownGeogCS('EPSG:4326')
#    dst_ds.SetProjection(proj.ExportToWkt())
     dst_ds.SetProjection('')
#    dst_ds.SetGeoTransform((57.0, 0.10000000149011612, 0.0, -6.0, 0.0, -0.10000000149011612))
     dst_ds.SetGeoTransform(geo)
     dst_ds.GetRasterBand(1).WriteArray(varw)
#    dst_ds.FlushCache()
     dst_ds=None
     return


def puttif(filename,var,geot,sp,nband,nodatav):

  DRIVER=gdal.GetDriverByName('GTiff')
  NROWS,NCOLS=var.shape
  gtype=gdal_array.NumericTypeCodeToGDALTypeCode(var.dtype)
 #dataset = DRIVER.Create(filename, NCOLS, NROWS, nband, gtype)
  dataset = DRIVER.Create(filename, NCOLS, NROWS, nband, 6)
  dataset.SetGeoTransform(geot)
  dataset.SetProjection(sp)
  dataset.GetRasterBand(1).WriteArray(var)
  dataset.GetRasterBand(1).SetNoDataValue(nodatav)
# dataset.FlushCache()
  dataset=None



#def writemap(filename,var,wformat):
#     os.system('mapattr -R {} -C {} {}'.format(var.shape,filename)
