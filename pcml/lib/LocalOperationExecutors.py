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

@executor
@localoperation
def LocalMaximum_np(self, subdomains):
#finding Local Maximum among the given locations
	outsubdomain = subdomains[0]
	outarr = outsubdomain.get_nparray()	
	arr= np.maximum(subdomains[1].get_nparray(),subdomains[2].get_nparray())
	outarr[:,:]=arr

@executor
@localoperation
#finding Local Minimum among the given locations
def LocalMinimum_np(self, subdomains):
	outsubdomain = subdomains[0]
	outarr = outsubdomain.get_nparray()
	arr = np.minimum(subdomains[1].get_nparray(),subdomains[2].get_nparray())
	outarr[:,:]=arr
	
@executor
@localoperation
#finding Local Mean among the given locations
def LocalMean_np(self, subdomains):
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()
    arr=np.add(subdomains[1].get_nparray(),subdomains[2].get_nparray())
    denom=len(subdomains)-1
    outarr[:,:]=arr/denom


@executor
@localoperation
#finding Local Difference among the given locations
def LocalDifference_np(self, subdomains):
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()
    arr1=np.array(subdomains[1].get_nparray())
    arr2=np.array(subdomains[2].get_nparray())
    arr=arr1-arr2
    outarr[:,:]=arr
    
@executor
@localoperation
#finding Local Product among the given locations
def LocalProduct_np(self, subdomains):
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()
    arr=np.multiply(subdomains[1].get_nparray(),subdomains[2].get_nparray())
    outarr[:,:]=arr

