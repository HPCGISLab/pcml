"""Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
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
def ZonalSum_exec(self, subdomains):
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()
    zoneslice=subdomains[2].slice_nparray(subdomains[0].r,subdomains[0].c,outsubdomain.nrows,outsubdomain.ncols)
    zonalarray=subdomains[2].get_nparray()
    rasterarray=subdomains[1].get_nparray()
    zones=np.unique(zoneslice)
    zoneindicesdict={}
    for zone in zones:
        zoneindicesdict[zone]=np.sum(rasterarray[np.where(zonalarray==zone)])
    vals=[]
    for zonedata in zoneslice.flat:
        vals.append(zoneindicesdict[zonedata])
    outarr[:,:]=np.asarray(vals).reshape(outarr.shape)

@executor
@zonaloperation
def ZonalMean_exec(self, subdomains):
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()
    zoneslice=subdomains[2].slice_nparray(subdomains[0].r,subdomains[0].c,outsubdomain.nrows,outsubdomain.ncols)
    zonalarray=subdomains[2].get_nparray()
    rasterarray=subdomains[1].get_nparray()
    zones=np.unique(zoneslice)
    zoneindicesdict={}
    for zone in zones:
        zoneindicesdict[zone]=np.mean(rasterarray[np.where(zonalarray==zone)])
    vals=[]
    for zonedata in zoneslice.flat:
        vals.append(zoneindicesdict[zonedata])
    outarr[:,:]=np.asarray(vals).reshape(outarr.shape)

@executor
@zonaloperation
def ZonalMaximum_exec(self, subdomains):
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()
    zoneslice=subdomains[2].slice_nparray(subdomains[0].r,subdomains[0].c,outsubdomain.nrows,outsubdomain.ncols)
    zonalarray=subdomains[2].get_nparray()
    rasterarray=subdomains[1].get_nparray()
    zones=np.unique(zoneslice)
    zoneindicesdict={}
    for zone in zones:
        zoneindicesdict[zone]=np.max(rasterarray[np.where(zonalarray==zone)])
    vals=[]
    for zonedata in zoneslice.flat:
        vals.append(zoneindicesdict[zonedata])
    outarr[:,:]=np.asarray(vals).reshape(outarr.shape)


@executor
@zonaloperation
def ZonalMinimum_exec(self, subdomains):
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()
    zoneslice=subdomains[2].slice_nparray(subdomains[0].r,subdomains[0].c,outsubdomain.nrows,outsubdomain.ncols)
    zonalarray=subdomains[2].get_nparray()
    rasterarray=subdomains[1].get_nparray()
    zones=np.unique(zoneslice)
    zoneindicesdict={}
    for zone in zones:
        zoneindicesdict[zone]=np.min(rasterarray[np.where(zonalarray==zone)])
    vals=[]
    for zonedata in zoneslice.flat:
        vals.append(zoneindicesdict[zonedata])
    outarr[:,:]=np.asarray(vals).reshape(outarr.shape)


@executor
@zonaloperation
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

