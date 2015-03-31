"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
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

@zonaloperation
def zonalmean(self,locations,subdomains):
    rasterarray=subdomains[0].bufferedlocgetarr(locations[0],self.buffersize)
    zonalarray=subdomains[1].bufferedlocgetarr(locations[1],self.buffersize)
    return np.mean(rasterarray[np.where(zonalarray==locations[1]['v'])])
@zonaloperation
def zonalsum(self,locations,subdomains):
    rasterarray=subdomains[0].bufferedlocgetarr(locations[0],self.buffersize)
    zonalarray=subdomains[1].bufferedlocgetarr(locations[1],self.buffersize)
    return np.sum(rasterarray[np.where(zonalarray==locations[1]['v'])])

@zonaloperation
def zonalmaximum(self,locations,subdomains):
    rasterarray=subdomains[0].bufferedlocgetarr(locations[0],self.buffersize)
    zonalarray=subdomains[1].bufferedlocgetarr(locations[1],self.buffersize)
    return np.max(rasterarray[np.where(zonalarray==locations[1]['v'])])

@zonaloperation
def zonalminimum(self,locations,subdomains):
    rasterarray=subdomains[0].bufferedlocgetarr(locations[0],self.buffersize)
    zonalarray=subdomains[1].bufferedlocgetarr(locations[1],self.buffersize)
    return np.min(rasterarray[np.where(zonalarray==locations[1]['v'])])

@zonaloperation
def zonalmajority(self,locations,subdomains):
    rasterarray=subdomains[0].bufferedlocgetarr(locations[0],self.buffersize)
    zonalarray=subdomains[1].bufferedlocgetarr(locations[1],self.buffersize)
    return stats.mode(rasterarray[np.where(zonalarray==locations[1]['v'])])[0][0]

@zonaloperation
def zonalminority(self,locations,subdomains):
    rasterarray=subdomains[0].bufferedlocgetarr(locations[0],self.buffersize)
    zonalarray=subdomains[1].bufferedlocgetarr(locations[1],self.buffersize)
    frequency=stats.itemfreq(rasterarray[np.where(zonalarray==locations[1]['v'])])
    return frequency[np.argmin(frequency[:,1])][0]
