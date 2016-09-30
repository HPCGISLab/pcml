"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Jayakrishnan Ajayakumar (jajayaku@kent.edu) 
"""
from ..core.Operation import *
from ..core.Scheduler import *
from ..util.OperationBuilder import *
#import numpy as np
#import types
#import math

if PCMLConfig.scipyenabled:
    from scipy import stats
	


#Calculate the zonal mean at a location based on two input subdomains with raster data and zonal data
@zonaloperation
def zonalmean(self,locations,subdomains):
    rasterarray=subdomains[0].bufferedlocgetarr(locations[0],self.buffersize)
    zonalarray=subdomains[1].bufferedlocgetarr(locations[1],self.buffersize)
    return np.mean(rasterarray[np.where(zonalarray==locations[1]['v'])])

#Calculate the zonal sum at a location based on two input subdomains with raster data and zonal data
@zonaloperation
def zonalsum(self,locations,subdomains):
    rasterarray=subdomains[0].bufferedlocgetarr(locations[0],self.buffersize)
    zonalarray=subdomains[1].bufferedlocgetarr(locations[1],self.buffersize)
    return np.sum(rasterarray[np.where(zonalarray==locations[1]['v'])])

#Calculate the zonal maximum at a location based on two input subdomains with raster data and zonal data
@zonaloperation
def zonalmaximum(self,locations,subdomains):
    rasterarray=subdomains[0].bufferedlocgetarr(locations[0],self.buffersize)
    zonalarray=subdomains[1].bufferedlocgetarr(locations[1],self.buffersize)
    return np.max(rasterarray[np.where(zonalarray==locations[1]['v'])])

#Calculate the zonal minumum at a location based on two input subdomains with raster data and zonal data
@zonaloperation
def zonalminimum(self,locations,subdomains):
    rasterarray=subdomains[0].bufferedlocgetarr(locations[0],self.buffersize)
    zonalarray=subdomains[1].bufferedlocgetarr(locations[1],self.buffersize)
    return np.min(rasterarray[np.where(zonalarray==locations[1]['v'])])

if PCMLConfig.scipyenabled:
    #Calculate the zonal majority at a location based on two input subdomains with raster data and zonal data
    @zonaloperation
    def zonalmajority(self,locations,subdomains):
        rasterarray=subdomains[0].bufferedlocgetarr(locations[0],self.buffersize)
        zonalarray=subdomains[1].bufferedlocgetarr(locations[1],self.buffersize)
        print("outzm=",stats.mode(rasterarray[np.where(zonalarray==locations[1]['v'])])[0][0])
        return stats.mode(rasterarray[np.where(zonalarray==locations[1]['v'])])[0][0]

if PCMLConfig.scipyenabled:
    #Calculate the zonal minority at a location based on two input subdomains with raster data and zonal data
    @zonaloperation
    def zonalminority(self,locations,subdomains):
        rasterarray=subdomains[0].bufferedlocgetarr(locations[0],self.buffersize)
        zonalarray=subdomains[1].bufferedlocgetarr(locations[1],self.buffersize)
        frequency=stats.itemfreq(rasterarray[np.where(zonalarray==locations[1]['v'])])
        return frequency[np.argmin(frequency[:,1])][0]
