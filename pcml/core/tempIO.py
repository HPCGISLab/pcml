import pcml.core.PCMLConfig as PCMLConfig
#from linecache import getline
import fileinput
import numpy as np
try:
   PCMLConfig.osgeoenabled=1
   from osgeo import gdal
   from osgeo import ogr
   from osgeo import osr
except ImportError as e:
   PCMLConfig.osgeoenabled=0
   if e.message != 'No module named osgeo':
      raise
#This method should be in OperationIO
#This method takes in a filename and layer as input parameter and creates a blank tiff as output
def WriteGeoTIFFNew(filename, layer):
    if PCMLConfig.osgeoenabled==0:
       return None
    offset=1      #rows to write . This can be increased to write more rows at a time
    driver = gdal.GetDriverByName('GTiff')
    out = driver.Create(filename, layer.ncols, layer.nrows, 1, gdal.GDT_Float64)
    #out = driver.Create(filename, layer.ncols, layer.nrows, 1, gdal.GDT_Int32)
    out.SetGeoTransform((layer.x, layer.cellsize, 0, layer.y, 0, layer.cellsize))
    outband = out.GetRasterBand(1)
    #set no data value
    outband.SetNoDataValue(layer.nodata_value)
    #adding rows one by one,creating a big array in a single step creates Memory Error
    #adding temp array to improve performance
    newarr=np.zeros((offset,layer.ncols))
    for i in xrange(layer.nrows):
        outband.WriteArray(newarr,0,i*offset)
    outSRS = osr.SpatialReference()
    outSRS.ImportFromEPSG(4326)
    out.SetProjection(outSRS.ExportToWkt())
    outband.FlushCache()
