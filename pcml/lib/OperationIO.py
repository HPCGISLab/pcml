"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from ..core.Layer import *
import pcml.core.PCMLConfig as PCMLConfig

import numpy as np
#from linecache import getline
import fileinput


try: 
   PCMLConfig.osgeoenabled=1
   from osgeo import gdal
   from osgeo import ogr
   from osgeo import osr
except ImportError as e:
   PCMLConfig.osgeoenabled=0
   if e.message != 'No module named osgeo':
      raise

# Assumes the following format
'''
ncols         4
nrows         6
xllcorner     0.0
yllcorner     0.0
cellsize      50.0
NODATA_value  -9999
'''


def ReadASCIIGrid(filename):

    nrows=ncols=None
    x=y=None
    cellsize=None
    nodata_value=None
    for line in fileinput.input([filename]):
        arg, val = str.split(line)
        if arg == "nrows":
            nrows = int(val)
        if arg == "ncols":
            ncols = int(val)
        if arg == "xllcorner":
            x = float(val)
        if arg == "yllcorner":
            y = float(val)
        if arg == "cellsize":
            cellsize=float(val)
        if arg == "NODATA_value":
            nodata_value=float(val)
        if fileinput.filelineno()>=6:
            break
    fileinput.close()

    assert(nrows!=None)
    assert(ncols!=None)
    assert(y!=None)
    assert(x!=None)
    assert(cellsize!=None)
    assert(nodata_value!=None)

    arg=None
    val=None

    # TODO: Here we should check to see if all 6 values are set
    h=float(nrows)*cellsize
    w=float(ncols)*cellsize
    layer=Layer(y,x,h,w,filename)
    nparr=np.loadtxt(filename, skiprows=6)
    layer.set_nparray(nparr,cellsize,nodata_value)

    del nparr

    return layer


# Ugly but functional
def WriteASCIIGrid(filename, layer):
    assert(layer.data_structure==Datastructure.array)
    string = "ncols        %i\n" % layer.ncols
    string += "nrows        %i\n" % layer.nrows
    string += "xllcorner    %f\n" % layer.x
    string += "yllcorner    %f\n" % layer.y
    string += "cellsize     %.15f\n" % layer.cellsize
    string += "NODATA_value %f\n" % layer.nodata_value
    arr = layer.get_nparray()
    assert(layer.data_structure==Datastructure.array)
    for i in xrange(layer.nrows):
        for j in xrange(layer.ncols):
            string += (PCMLConfig.value_precision + ' ') % arr[i][j]
        string += "\n"
    asciigridfile = open(filename, "w")
    asciigridfile.write(string)
    asciigridfile.close()

def ReadGeoTIFF(filename,bandnumber=1):
    if PCMLConfig.osgeoenabled==0:
       PCMLUserInformation("ReadGeoTIFF is disabled, because PCML could not find osgeo or gdal library")
       return None 

    ds = gdal.Open(filename) # Open gdal dataset
    if ds is None:
        raise PCMLException("Cannot open "+filename+" in ReadGeoTIFF")

    # By default get the first band
    band = ds.GetRasterBand(bandnumber)
    ncols = ds.RasterXSize
    nrows = ds.RasterYSize

    if band is None:
        raise PCMLException("Cannot read selected band in "+filename+" in ReadGeoTIFF")

    nodata_value = band.GetNoDataValue()

    transform = ds.GetGeoTransform()
    cellsize = transform[1]
    origin_x = transform[0]
    origin_y = transform[3]
    x=origin_x
    y=origin_y-nrows*cellsize
    if(abs(transform[1])!=abs(transform[5])): # pixelwidth=1, pixelheight=5
        PCMLUserInformation("Cells of different height and width selecting width not height")

    h=float(nrows)*cellsize
    w=float(ncols)*cellsize
    layer=Layer(y,x,h,w,filename)
    nparr=band.ReadAsArray(0,0,ncols,nrows) 
    layer.set_nparray(nparr,cellsize,nodata_value)

    del transform
    del nparr
    del band
    del ds
    nparr=None
    ds=None # Close gdal dataset
    band=None

    return layer

def WriteGeoTIFF(filename, layer):
    if PCMLConfig.osgeoenabled==0:
       PCMLUserInformation("WriteGeoTIFF is disabled, because PCML could not find osgeo or gdal library")
       return None 

    assert(layer.data_structure==Datastructure.array)

    driver = gdal.GetDriverByName('GTiff')
    out = driver.Create(filename, layer.ncols, layer.nrows, 1, gdal.GDT_CFloat64)
    if out is None:
        raise PCMLException("Cannot open '"+filename+"' to write")
    out.SetGeoTransform((layer.x, layer.cellsize, 0, layer.y, 0, layer.cellsize))
    outband = out.GetRasterBand(1)
    outband.WriteArray(layer.get_nparray())
    outSRS = osr.SpatialReference()
    outSRS.ImportFromEPSG(4326)
    out.SetProjection(outSRS.ExportToWkt())
    outband.FlushCache()

