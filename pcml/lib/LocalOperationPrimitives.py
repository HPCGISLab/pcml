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

@localoperation
def LocalSum(self, locations, subdomains):
    # LocalSum will return the sum of values at each location for two or more layers.
    val = 0
    for loc in locations:
        val += loc['v']
    return val

@localoperation
def LocalMult(self, locations, subdomains):
    val = 1
    for loc in locations:
        val *= loc['v']
    return val

@localoperation
def LocalDivision(self, locations, subdomains):
    # NOTE: Assumes 2 subdomains for division.
    location0=locations[0]
    location1=locations[1]
    return location0['v']/location1['v']


@localoperation
def LocalSubtraction(self, locations, subdomains):
    # NOTE: Assumes 2 subdomains for subtraction.
    location0=locations[0]
    location1=locations[1]
    return location0['v']-location1['v']


@localoperation
def LocalSin(self, locations, subdomains):
    # NOTE: Assumes 1 subdomain
    return math.sin(locations[0]['v'])



