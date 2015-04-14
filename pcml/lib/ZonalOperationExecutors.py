"""Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Jayakrishnan Ajayakumar (jajayaku@kent.edu)  
"""
from ..core.Operation import *
from ..core.Scheduler import *
from ..util.OperationBuilder import *
import numpy as np
import types
import math
from scipy import stats
@executor
@zonaloperation
#Calculate the zonal sum based on two input subdomains with raster data and zonal data
def ZonalSum_exec(self, subdomains):
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()
    zoneslice=subdomains[2].slice_nparray(subdomains[0].r,subdomains[0].c,outsubdomain.nrows,outsubdomain.ncols)
    zonalarray=subdomains[2].get_nparray()
    rasterarray=subdomains[1].get_nparray()
    zones=np.unique(zoneslice)
    zoneindicesdict={}
    #code changes to reuse cache.Looping and applying numpy operation to enhance performance
    for zone in zones:
        totalsum=0
        for j in xrange(subdomains[2].nrows):
            totalsum+=np.sum(rasterarray[j,][np.where(zonalarray[j,]==zone)])
        zoneindicesdict[zone]=totalsum
    vals=[]
    for zonedata in zoneslice.flat:
        vals.append(zoneindicesdict[zonedata])
    outarr[:,:]=np.asarray(vals).reshape(outarr.shape)

@executor
@zonaloperation
#Calculate the zonal mean based on two input subdomains with raster data and zonal data
def ZonalMean_exec(self, subdomains):
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()
    zoneslice=subdomains[2].slice_nparray(subdomains[0].r,subdomains[0].c,outsubdomain.nrows,outsubdomain.ncols)
    zonalarray=subdomains[2].get_nparray()
    rasterarray=subdomains[1].get_nparray()
    zones=np.unique(zoneslice)
    zoneindicesdict={}
    #code changes to reuse cache.Looping and applying numpy operation to enhance performance
    for zone in zones:
        totalsum=0
        for j in xrange(subdomains[2].nrows):
            totalsum+=np.sum(rasterarray[j,][np.where(zonalarray[j,]==zone)])
        zoneindicesdict[zone]=totalsum/zonalarray.size
    vals=[]
    for zonedata in zoneslice.flat:
        vals.append(zoneindicesdict[zonedata])
    outarr[:,:]=np.asarray(vals).reshape(outarr.shape)

@executor
@zonaloperation
#Calculate the zonal maximum based on two input subdomains with raster data and zonal data
def ZonalMaximum_exec(self, subdomains):
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()
    zoneslice=subdomains[2].slice_nparray(subdomains[0].r,subdomains[0].c,outsubdomain.nrows,outsubdomain.ncols)
    zonalarray=subdomains[2].get_nparray()
    rasterarray=subdomains[1].get_nparray()
    zones=np.unique(zoneslice)
    zoneindicesdict={}
    #code changes to reuse cache.Looping and applying numpy operation to enhance performance
    for zone in zones:
        maxval=np.NINF
        for j in xrange(subdomains[2].nrows):
            maxval=max(np.sum(rasterarray[j,][np.where(zonalarray[j,]==zone)]),maxval)
        zoneindicesdict[zone]=maxval
    vals=[]
    for zonedata in zoneslice.flat:
        vals.append(zoneindicesdict[zonedata])

    outarr[:,:]=np.asarray(vals).reshape(outarr.shape)


@executor
@zonaloperation
#Calculate the zonal minimum based on two input subdomains with raster data and zonal data
def ZonalMinimum_exec(self, subdomains):
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()
    zoneslice=subdomains[2].slice_nparray(subdomains[0].r,subdomains[0].c,outsubdomain.nrows,outsubdomain.ncols)
    zonalarray=subdomains[2].get_nparray()
    rasterarray=subdomains[1].get_nparray()
    zones=np.unique(zoneslice)
    zoneindicesdict={}
    #code changes to reuse cache.Looping and applying numpy operation to enhance performance
    for zone in zones:
        minval=np.inf
        for j in xrange(subdomains[2].nrows):
            minval=min(np.sum(rasterarray[j,][np.where(zonalarray[j,]==zone)]),minval)
        zoneindicesdict[zone]=minval
    vals=[]
    for zonedata in zoneslice.flat:
        vals.append(zoneindicesdict[zonedata])

    outarr[:,:]=np.asarray(vals).reshape(outarr.shape)


@executor
@zonaloperation
#Calculate the zonal majority based on two input subdomains with raster data and zonal data
def ZonalMajority_exec(self, subdomains):
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()
    zoneslice=subdomains[2].slice_nparray(subdomains[0].r,subdomains[0].c,outsubdomain.nrows,outsubdomain.ncols)
    zonalarray=subdomains[2].get_nparray()
    rasterarray=subdomains[1].get_nparray()
    zones=np.unique(zoneslice)
    zoneindicesdict={}
    for zone in zones:
        zoneindicesdict[zone]=stats.mode(rasterarray[np.where(zonalarray==zone)])[0][0]
    vals=[]
    for zonedata in zoneslice.flat:
        vals.append(zoneindicesdict[zonedata])
    outarr[:,:]=np.asarray(vals).reshape(outarr.shape)

@executor
@zonaloperation
#Calculate the zonal minority based on two input subdomains with raster data and zonal data
def ZonalMinority_exec(self, subdomains):
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()
    zoneslice=subdomains[2].slice_nparray(subdomains[0].r,subdomains[0].c,outsubdomain.nrows,outsubdomain.ncols)
    zonalarray=subdomains[2].get_nparray()
    rasterarray=subdomains[1].get_nparray()
    zones=np.unique(zoneslice)
    zoneindicesdict={}
    for zone in zones:
        frequency=stats.itemfreq(rasterarray[np.where(zonalarray==zone)])
        zoneindicesdict[zone]=frequency[np.argmin(frequency[:,1])][0]
    vals=[]
    for zonedata in zoneslice.flat:
        vals.append(zoneindicesdict[zonedata])
    outarr[:,:]=np.asarray(vals).reshape(outarr.shape)

