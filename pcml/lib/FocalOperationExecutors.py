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
#from scipy import stats

@executor
@focaloperation
def FocalMean_np_exec(self, subdomains):
    # Get the array from the output subdomain as well as the output subdomain
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray() 
    # Get the input subdomain
    insubdomain = subdomains[1]
    inarr = insubdomain.get_nparray() 

    # Get the buffersize for the focal operation
    buffersize=self.buffersize

    locind={}

    for roworig in xrange(outsubdomain.nrows):
        for colorig in xrange(outsubdomain.ncols):
            # Calculate the row and column positions (global for the study area)
            # which is why we add outsubdomain.{r,c}
            locind['r']=roworig+outsubdomain.r
            locind['c']=colorig+outsubdomain.c

            # The following 2 lines of code are the easy-to-read version of the code below
            # Essentially we get the location for the input array and then average it
            # Get the array from the location
            #arr=insubdomain.bufferedlocgetarr(locind,buffersize)
            # Calculate the mean and set it
            #outarr[roworig][colorig]=np.average(arr)

            # Add the offsets for the subdomains (difference etween out/input subdomains) 
            row=locind['r']-insubdomain.r
            col=locind['c']-insubdomain.c
            # Calculate the starting row and column accounting for buffersize and the edge of layer
            r=max(0,row-buffersize)
            c=max(0,col-buffersize)

            # Calculate the h accounting for buffersize and edge of layer
            h=buffersize+(row-r)+1
            if (r+h > insubdomain.nrows):
                h=insubdomain.nrows-r

            # Same for width 
            w=buffersize+(col-c)+1
            if (c+w > insubdomain.ncols):
                w=insubdomain.ncols-c

            # Take a slice of the array (inarr) based on those calculations, average the results
            outarr[roworig][colorig]=np.average(inarr[r:r + h, c:c + w])

