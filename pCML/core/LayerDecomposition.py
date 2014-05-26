"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from .Layer import *
from .Subdomain import *

import math

# Defines the number of subdomains (e.g., rows) to decompose a single layer into.
numchunksgoal=2
numchunksgoal=16

# Different decomposition methods
class DecompositionMethod():
    """ Enumeration class. It defines the methods supported for layer decomposition.
    """
    row = 1
    column = 2
    quadtree = 3
    grid = 4

# Take a layer and return a list of subdomains 
def rowdecomposition(layer, buffersize):

    if layer.data_structure!=Datastructure.array:
        if layer.data_structure==Datastructure.pointlist and buffersize < 0: # globalclass+pointlist 
            pass
        else:
            PCMLNotSupported("Currently rowdecomposition assumes array data_structure or globalclass+pointlist")

    if buffersize<0: # This indicates the buffer should be infinite sized (global/zonal operation)
        buffersize=99999999

    assert(layer.nrows!=None)
    assert(layer.ncols!=None)

    subdomainlist = []

    # Number of rows in array
    numrows=layer.nrows

    # Numer of rows per subdomain 
    rowspersubdomain = int(math.ceil(float(numrows)/float(numchunksgoal)))

    # Number of subdomains to create
    numsubdomains = int(math.ceil(float(numrows)/float(rowspersubdomain)))

    if layer.data_structure==Datastructure.pointlist:
        # If this layer is a pointlist, then it is assumed to be global operation
        # so just copy layer information and duplicate the subdomain
        subdomain = Subdomain(layer.y, layer.x, layer.h, layer.w, layer.title+" subdomain pointlist")
        subdomain.set_pointlist(layer.get_pointlist())
        subdomainlist=[]
        for sdind in xrange(numsubdomains):
            subdomainlist.append(subdomain)
        return subdomainlist

    for sdind in xrange(numsubdomains):
        r = rowspersubdomain*sdind
        nrows = rowspersubdomain # Number of rows for this sudomain

        if r+nrows >= layer.nrows: # Problem, because the subdomain is too big
           nrows = layer.nrows-r # Set nrows to be the difference between layer.nrows and r 

        if buffersize>0: # Focal operation, need to adjust r,nrows
           # Need to calculate how many rows we added to r based on buffersize
           # When accounting for 0 boundary, this is stored in rdiff
           oldr=r                # Record old r value
           r=max(0,r-buffersize) # Calculate new r value
           oldh=nrows            # Record old height
           newh=min(layer.nrows,oldr+oldh+buffersize) # Calculate new height (bounded by layer.nrows)
           nrows=newh-r # New nrows is based on newh and new r

        # In row decomposition, column index is always 0 and ncols never changes
        c = 0 
        ncols = layer.ncols  

        # Now derive y, x, h, w
        y = layer.y + r * layer.cellsize
        h = nrows * layer.cellsize 
        # In row decomposition: x and w always remain the same
        x = layer.x
        w = layer.w

        subdomain = Subdomain(y, x, h, w, layer.title+" subdomain "+str(sdind))

        subdomain.cellsize=layer.cellsize
        subdomain.nodata_value=layer.nodata_value
        subdomain.r=r
        subdomain.c=c
        subdomain.nrows=nrows
        subdomain.ncols=ncols

        arrslice=layer.slice_nparray(r,0,nrows,ncols)
        subdomain.set_data_ref(arrslice)

        subdomainlist.append(subdomain)

    return subdomainlist 

