"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from .Layer import *
from .Subdomain import *
import pcml.core.PCMLConfig as PCMLConfig

import math

# Defines the number of subdomains (e.g., rows) to decompose a single layer into. (DEPRECATED)
#numchunksgoal=2
#numchunksgoal=16

def globalpointlistdecomposition(layer, buffersize):
    # Row decomposition supports pointlist only for globalclass operations
    if layer.data_structure==Datastructure.pointlist:

        if buffersize>=0: # Then it is not globalclass operation
            raise PCMLNotSupported("Currently globalpointlistdecomposition only supports globalclass+pointlist")

        # If this layer is a pointlist, then it is assumed to be global operation
        # so just copy layer information and duplicate the subdomain
        subdomain = Subdomain(layer.y, layer.x, layer.h, layer.w, layer.title+" subdomain pointlist")
        subdomain.set_pointlist(layer.get_pointlist())
        subdomainlist=[]
        for sdind in xrange(numsubdomains):
            subdomainlist.append(subdomain)
        return subdomainlist
    else:
        raise PCMLNotSupported("globalpointlistdecomposition only supports pointlist datastructures")

# Take a layer and return a list of subdomains 
def rowdecomposition(layer, buffersize):

    #print("Row decomposition")

    # Row decomposition supports pointlist only for globalclass operations
    if layer.data_structure==Datastructure.pointlist:
        globalpointlistdecomposition(layer,buffersize)

    assert(layer.data_structure==Datastructure.array),"Data structure is not an array"

    # If global then buffer size is infinite as all subdomains will have all data
    if buffersize<0: # This indicates the buffer should be infinite sized (global/zonal operation)
        buffersize=9999999999999
        # FIXME: I should do the same global subdomain as pointlist here

    # Sanity check nrows and ncols
    # FIXME: In the future this check will happen in operation._decompositioninit
    assert(layer.nrows!=None),"Layer number of rows (nrows) is None"
    assert(layer.ncols!=None),"Layer number of columns (ncols) is None"

    subdomainlist = []

    # Numer of rows per subdomain given suggested decomposition granularity (think number of chunks)
    #rowspersubdomain = int(math.ceil(float(layer.nrows)/float(PCMLConfig.decomposition_granularity)))

    # Number of subdomains to create when given rowspersubdomain
    numsubdomains = int(math.ceil(float(layer.nrows)/float(PCMLConfig.decomposition_granularity)))

    # For each subdomain indexed by sdind, calculate the size
    for sdind in xrange(numsubdomains):
        # First row in the subdomain
        r = PCMLConfig.decomposition_granularity*sdind

        # Default number of rows for this subdomain
        nrows = PCMLConfig.decomposition_granularity # Number of rows for this sudomain

        if buffersize>0: # If we have a buffer (e.g., focal operation), then add the buffer 
           # A buffer will generally reduce r by buffersize and increase nrows by buffersize*2
           # However, r and r+nrows must be contained within the range 0-layer.nrows
           new_r=max(0,r-buffersize) # Calculate new r value making sure it is not negative
           new_h=min(layer.nrows,r+nrows+buffersize) # calculate new height making sure it is <= layer.nrows

           # Replace original r and nrows with new values
           nrows=new_h-new_r
           r=new_r
           #print("new_r",new_r,"new_h",new_h)
        else: # Ensure that we don't allocate more rows past the number of layer rows
            nrows=min(layer.nrows-r,nrows)

        # Sanity check
        #print("r",r,"nrows",nrows,"layer.nrows",layer.nrows)
        assert(r+nrows<=layer.nrows),"Number of rows for layer is less than total for subdomains"

        # In row decomposition, column index is always 0 and ncols never changes
        c = 0 
        ncols = layer.ncols  

        # Now derive y, x, h, w
        y = layer.y + r * layer.cellsize
        h = nrows * layer.cellsize 
        # In row decomposition: x and w always remain the same
        x = layer.x
        w = layer.w

        # Create a subdomain and populate it with the correct attribute values
        subdomain = Subdomain(y, x, h, w, layer.title+" subdomain "+str(sdind))
        subdomain.cellsize=layer.cellsize
        subdomain.nodata_value=layer.nodata_value
        subdomain.r=r
        subdomain.c=c
        subdomain.nrows=nrows
        subdomain.ncols=ncols

        # Extract an array slice (reference to data in a layer for lower memory overhead)
        # from the layer and set the data reference for the subdomain to use
        arrslice=layer.slice_nparray(r,0,nrows,ncols)
        subdomain.set_data_ref(arrslice)

        # Add the subdomain to the list
        subdomainlist.append(subdomain)

    return subdomainlist 


# Take a layer and return a list of subdomains 
def columndecomposition(layer, buffersize):

    #print("Column decomposition")

    # Col decomposition supports pointlist only for globalclass operations
    if layer.data_structure==Datastructure.pointlist:
        globalpointlistdecomposition(layer,buffersize)

    assert(layer.data_structure==Datastructure.array),"Data structure is not an array"

    # If global then buffer size is infinite as all subdomains will have all data
    if buffersize<0: # This indicates the buffer should be infinite sized (global/zonal operation)
        buffersize=9999999999999
        # FIXME: I should do the same global subdomain as pointlist here

    # Sanity check nrows and ncols
    # FIXME: In the future this check will happen in operation._decompositioninit
    assert(layer.nrows!=None),"Layer number of rows (nrows) is None"
    assert(layer.ncols!=None),"Layer number of columns (ncols) is None"

    subdomainlist = []

    # Numer of columns per subdomain given suggested decomposition granularity (think number of chunks) 
    #colspersubdomain = int(math.ceil(float(layer.ncols)/float(PCMLConfig.decomposition_granularity)))

    # Number of subdomains to create when given colspersubdomain
    numsubdomains = int(math.ceil(float(layer.ncols)/float(PCMLConfig.decomposition_granularity)))

    # For each subdomain indexed by sdind, calculate the size
    for sdind in xrange(numsubdomains):
        # First col in the subdomain
        c = PCMLConfig.decomposition_granularity*sdind

        # Default number of columns for this subdomain
        ncols = PCMLConfig.decomposition_granularity # Number of columns for this sudomain

        if buffersize>0: # If we have a buffer (e.g., focal operation), then add the buffer 
           # A buffer will generally reduce c by buffersize and increase ncols by buffersize*2
           # However, c and c+ncols must be contained within the range 0-layer.ncols
           new_c=max(0,c-buffersize) # Calculate new c value making sure it is not negative
           new_w=min(layer.ncols,c+ncols+buffersize) # calculate new width making sure it is <= layer.ncols

           # Replace original c and ncols with new values
           ncols=new_w-new_c
           c=new_c
        else: # Ensure that we don't allocate more cols than the cols in a layer
            ncols=min(layer.ncols-c,ncols)

        # Sanity check
        assert(c+ncols<=layer.ncols),"Number of columns in layer is less than total for subdomains"

        # In column decomposition, row index is always 0 and nrows never changes
        r = 0 
        nrows = layer.nrows 

        # Now derive y, x, h, w
        x = layer.x + c * layer.cellsize
        w = ncols * layer.cellsize
        # In column decomposition: y and h always remain the same
        y = layer.y
        h = layer.h 

        # Create a subdomain and populate it with the correct attribute values
        subdomain = Subdomain(y, x, h, w, layer.title+" subdomain "+str(sdind))
        subdomain.cellsize=layer.cellsize
        subdomain.nodata_value=layer.nodata_value
        subdomain.r=r
        subdomain.c=c
        subdomain.nrows=nrows
        subdomain.ncols=ncols

        # Extract an array slice (reference to data in a layer for lower memory overhead)
        # from the layer and set the data reference for the subdomain to use
        arrslice=layer.slice_nparray(0,c,nrows,ncols)
        subdomain.set_data_ref(arrslice)

        # Add the subdomain to the list
        subdomainlist.append(subdomain)

    return subdomainlist 



