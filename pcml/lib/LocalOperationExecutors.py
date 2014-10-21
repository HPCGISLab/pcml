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

@executor
@localoperation
def LocalSum_np(self, subdomains):
    # NOTE: Assumes 3 subdomains, first is output, second and third should be added
    # Get the array from the output subdomain
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()

    # Apply numpy operation to arrays from second and third subdomains
    arr=np.add(subdomains[1].get_nparray(),subdomains[2].get_nparray())
    # Copy values to outarr (outsubdomain)
    outarr[:,:]=arr
    # Notice we don't need to return anything, because the resulting array (arr) is copied to outsubdomain through outarr 



@executor
@localoperation
def LocalMult_np(self, subdomains):
    # NOTE: Assumes 3 subdomains, first is output, second and third should be added
    # Get the array from the output subdomain
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()

    # Apply numpy operation to arrays from second and third subdomains
    arr=np.multiply(subdomains[1].get_nparray(),subdomains[2].get_nparray())
    # Copy values to outarr (outsubdomain)
    outarr[:,:]=arr


