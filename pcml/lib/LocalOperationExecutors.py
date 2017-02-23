"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu); Gowtham Kukkadapu(gkukkada@kent.edu)
"""

from ..core.Operation import *
from ..core.Scheduler import *
from ..util.OperationBuilder import *
import numpy as np
import types
import math

# Try to load numba, which is a just-in-time (jit) compiler for Python,
# which can be used to dramatically speedup data processing even though the code is written in native python
numbaenabled=True
try:
   from numba import jit 
except ImportError as e:
   numbaenabled=False


if numbaenabled:
    # FIXME: Look at http://numba.pydata.org/numba-doc/0.3/doc/examples.html#id1
    #        numba example to make it optional for compilation
    #        I think I can say if numba_exists: ... else just call same python function without the jit recompilation
    #        to do this I need to remove the @jit decorator below and change it to an independent function call
    @jit(nopython=True)
    def numba_sum(outarr,leftarr,rightarr,nrows,ncols):
        for i in xrange(nrows):
            for j in xrange(ncols):
                outarr[i,j]=leftarr[i,j]+rightarr[i,j]
else:
    def numba_sum(outarr,leftarr,rightarr,nrows,ncols):
        for i in xrange(nrows):
            for j in xrange(ncols):
                outarr[i,j]=leftarr[i,j]+rightarr[i,j]

@executor
@localoperation
def LocalSum_numba(self, subdomains):
    # NOTE: Assumes 3 subdomains, first is output, second and third should be added
    # Get the array from the output subdomain
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray() 
    leftarr=subdomains[1].get_nparray()
    rightarr=subdomains[2].get_nparray()
    numba_sum(outarr,leftarr,rightarr,subdomains[1].nrows,subdomains[1].ncols)



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
	
    #get the arrays from subdomains . Assumes 2 subdomains    outsubdomain = subdomains[0]
	outsubdomain = subdomains[0]
	outarr = outsubdomain.get_nparray()
	arr = np.minimum(subdomains[1].get_nparray(),subdomains[2].get_nparray())
	outarr[:,:]=arr
	
    # Notice we don't need to return anything, because the resulting array (arr) is copied to outsubdomain through outarr

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

@executor
@localoperation
# Local Divide Executor
def LocalDivision_np(self, subdomains):
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()
    arr=np.divide(subdomains[1].get_nparray(),subdomains[2].get_nparray())
    outarr[:,:]=arr

@executor
@localoperation
# Local Exponent Executor
def LocalExponent_np(self, subdomains):
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()
    arr=np.exp(subdomains[1].get_nparray())
    outarr[:,:]=arr

@executor
@localoperation
# Local Natural Log Executor
def LocalLog_np(self, subdomains):
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()
    arr=np.log(subdomains[1].get_nparray())
    outarr[:,:]=arr

@executor
@localoperation
# Local Log base-2 Executor
def LocalLog_np(self, subdomains):
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()
    arr=np.log2(subdomains[1].get_nparray())
    outarr[:,:]=arr

# More utility operations
# Generally used to simplify output for visualization or final analysis

@executor
@localoperation
# Local Floor Executor
def LocalFloor_np(self, subdomains):
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()
    arr=np.floor(subdomains[1].get_nparray())
    outarr[:,:]=arr


@executor
@localoperation
# Local Ceiling Executor
def LocalCeil_np(self, subdomains):
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()
    arr=np.ceil(subdomains[1].get_nparray())
    outarr[:,:]=arr

@executor
@localoperation
# Local Truncatation Executor
def LocalTrunc_np(self, subdomains):
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()
    arr=np.trunc(subdomains[1].get_nparray())
    outarr[:,:]=arr



@executor
@localoperation
#creating Local Zones with parameter value of 20
def LocalClassify20_np(self, subdomains):
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()
    arr=np.floor(subdomains[1].get_nparray()/20)
    outarr[:,:]=arr

@executor
@localoperation
#finding Local Product among the given locations
def LocalClassify(self, subdomains):
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()

    classifyval = self.kwargs.get('classtype',0)

    arr=np.where(subdomains[1].get_nparray() == classifyval, 1, 0)
    outarr[:,:]=arr




