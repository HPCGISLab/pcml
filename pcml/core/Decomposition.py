"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from .Layer import *
from .Subdomain import *
import pcml.core.PCMLConfig as PCMLConfig
from .PCMLPrims import *
import math

# Take a list of layers and a buffer size, 
# create a list of (buffered) bounding boxes following decomposition strategy
def rowdecomposition(layers, buffersize):
    boundingboxlist = []
    # If output is raster, then use raster row decomposition
    if layers[0].data_structure==Datastructure.array:
        boundingboxlist = rowdecomposition_raster(layers[0],buffersize)
    # If output is a pointlist, then use point row decomposition
    elif layers[0].data_structure==Datastructure.pointlist:
        boundingboxlist = rowdecomposition_pointlist(layers[0],buffersize)
    else: 
        PCMLNotSupported("Currently rowdecomposition does not support this datatype")

    return boundingboxlist

def rowdecomposition_raster(layer, buffersize):

    # Sanity check nrows and ncols
    # FIXME: In the future this check will happen in operation._decompositioninit
    assert(layer.nrows!=None),"Layer number of rows (nrows) is None"
    assert(layer.ncols!=None),"Layer number of columns (ncols) is None"

    boundingboxlist = []

    # Number of subdomains to create when given rowspersubdomain
    numberof_boundingboxes = int(math.ceil(float(layer.nrows)/float(PCMLConfig.decomposition_granularity)))

    # Loop over each subdomain, and create a bounding box
    for bb_index in xrange(numberof_boundingboxes):
        # Calculate first row in subdomain using the formula below
        r = bb_index * PCMLConfig.decomposition_granularity

        # Default number of rows for this subdomain
        nrows = PCMLConfig.decomposition_granularity # Number of rows for this sudomain

        # If we exceed the number of rows in the last subdomain
        if r + nrows >= layer.nrows:
            nrows = layer.nrows - r # Resize the number of rows to "fit"
        
        # Now that we have r and nrows calculated we can calculate y,x h,w
        # In row decomposition: x and w always remain the same
        y = layer.y + r * layer.cellsize
        x = layer.x
        h = nrows * layer.cellsize 
        w = layer.w

        # Create the bounding box
        bb = BoundingBox(y, x, h, w)

        # Fill in the remaining data
        bb.nrows = nrows
        bb.ncols = layer.ncols
        bb.r = r
        bb.c = 0 

        # Add this bounding box to the list
        boundingboxlist.append(bb)

    cell_buffersize = int(math.ceil((buffersize) / layer.cellsize) ) # To convert the buffer to a cellsize, must do this division 

    # Loop over the newly created bounding box list and add buffers
    # This loop is separate to improve readability
    for bb_index in xrange(numberof_boundingboxes):
        bb = boundingboxlist[bb_index] # We want an index so we can modify the bounding boxes

        # In row decomposition, x, w, c, ncols are not buffered because they are full width
        bb.x_buf = bb.x
        bb.w_buf = bb.w
        bb.c_buf = bb.c
        bb.ncols_buf = bb.ncols
        
        # For most cases this will be the y, h, r, nrows
        bb.y_buf = bb.y - cell_buffersize * layer.cellsize 
        bb.r_buf = bb.r - cell_buffersize 
        bb.nrows_buf = bb.nrows + cell_buffersize * 2 # Remember we buffer both sides
        bb.h_buf = bb.nrows_buf * layer.cellsize

        # Test if we are over the bottom edge of the layer
        if bb.r_buf < 0:
            bb.nrows_buf += bb.r_buf # Remove those additional rows (remember r_buf is negative) 
            bb.r_buf = 0 # Set r_buf to edge of layer (zero)

        # Test if we are over the top edge of the layer
        if bb.r_buf + bb.nrows_buf >= layer.nrows:
            # Make it shorter by the difference of the two
            bb.nrows_buf -= (bb.r_buf + bb.nrows_buf) - layer.nrows

        # Recalculate y and h
        bb.y_buf = layer.y + bb.r_buf * layer.cellsize 
        bb.h_buf = bb.nrows_buf * layer.cellsize



        # FIXME: The _buf variables should be removed for now
        # FIXME: This may change when we implement MPI with ghost zones, but it doesn't make sense now.
        bb.x=bb.x_buf
        bb.y=bb.y_buf
        bb.h=bb.h_buf
        bb.w=bb.w_buf
        bb.nrows=bb.nrows_buf
        bb.ncols=bb.ncols_buf
        bb.r=bb.r_buf
        bb.c=bb.c_buf

        #print bb

    return boundingboxlist

def rowdecomposition_pointlist(layer, buffersize):
    PCMLNotSupported("rowdecomposition_pointlist is not yet supported")

# Take a list of layers and a buffer size, 
# create a list of (buffered) bounding boxes following decomposition strategy
def columndecomposition(layers, buffersize):
    boundingboxlist = []
    # If output is raster, then use raster row decomposition
    if layers[0].data_structure==Datastructure.array:
        boundingboxlist = coldecomposition_raster(layers[0],buffersize)
    else: 
        PCMLNotSupported("Currently coldecomposition does not support this datatype")

    return boundingboxlist

def coldecomposition_raster(layer, buffersize):

    # Sanity check nrows and ncols
    # FIXME: In the future this check will happen in operation._decompositioninit
    assert(layer.nrows!=None),"Layer number of rows (nrows) is None"
    assert(layer.ncols!=None),"Layer number of columns (ncols) is None"

    boundingboxlist = []

    # Number of subdomains to create when given decomposition granularity 
    numberof_boundingboxes = int(math.ceil(float(layer.ncols)/float(PCMLConfig.decomposition_granularity)))

    # Loop over each subdomain, and create a bounding box
    for bb_index in xrange(numberof_boundingboxes):
        # Calculate first column in subdomain using the formula below
        c = bb_index * PCMLConfig.decomposition_granularity

        # Default number of columns for this subdomain
        ncols = PCMLConfig.decomposition_granularity # Number of columns for this sudomain

        # If we exceed the number of columns in the last subdomain
        if c + ncols >= layer.ncols:
            ncols = layer.ncols - c # Resize the number of columns to "fit"
        
        # Now that we have c and ncols calculated we can calculate y,x h,w
        # In column decomposition: y and h always remain the same
        y = layer.y
        x = layer.x + c * layer.cellsize
        h = layer.h
        w = ncols * layer.cellsize 

        # Create the bounding box
        bb = BoundingBox(y, x, h, w)

        # Fill in the remaining data
        bb.nrows = layer.nrows
        bb.ncols = ncols
        bb.r = 0
        bb.c = c 

        # Add this bounding box to the list
        boundingboxlist.append(bb)

    cell_buffersize = int(math.ceil((buffersize) / layer.cellsize) ) # To convert the buffer to a cellsize, must do this division 

    # Loop over the newly created bounding box list and add buffers
    # This loop is separate to improve readability
    for bb_index in xrange(numberof_boundingboxes):
        bb = boundingboxlist[bb_index] # We want an index so we can modify the bounding boxes

        # In row decomposition, x, w, c, ncols are not buffered because they are full width
        bb.y_buf = bb.y
        bb.h_buf = bb.h
        bb.r_buf = bb.r
        bb.nrows_buf = bb.nrows
        
        # For most cases this will be the x, w, c, ncols
        bb.x_buf = bb.x - cell_buffersize * layer.cellsize 
        bb.c_buf = bb.c - cell_buffersize 
        bb.ncols_buf = bb.ncols + cell_buffersize * 2 # Remember we buffer both sides
        bb.w_buf = bb.ncols_buf * layer.cellsize

        # Test if we are over the left edge of the layer
        if bb.c_buf < 0:
            bb.ncols_buf += bb.c_buf # Remove those additional columns (remember r_buf is negative) 
            bb.c_buf = 0 # Set c_buf to edge of layer (zero)

        # Test if we are over the right edge of the layer
        if bb.c_buf + bb.ncols_buf >= layer.ncols:
            # Make it shorter by the difference of the two
            bb.ncols_buf -= layer.ncols - (bb.c_buf + bb.ncols_buf)

        # Recalculate x and w
        bb.x_buf = layer.x + bb.c_buf * layer.cellsize 
        bb.w_buf = bb.ncols_buf * layer.cellsize

        # FIXME: Same as row decomposition 
        bb.x=bb.x_buf
        bb.y=bb.y_buf
        bb.h=bb.h_buf
        bb.w=bb.w_buf
        bb.nrows=bb.nrows_buf
        bb.ncols=bb.ncols_buf
        bb.r=bb.r_buf
        bb.c=bb.c_buf
        '''
        '''

        print bb

    return boundingboxlist

'''
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
def rowdecompositionold(layer, buffersize):

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
def columndecompositionold(layer, buffersize):

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

#point decomposition using row strategy
def pointrowdecomposition(layer,buffersize):
    subdomainlist = []
    totalsubdomains=PCMLConfig.numsubdomains
    currenty=layer.y
    currentwithbuffy=layer.y
    layerblockheight=layer.h/totalsubdomains
    for subdindex in xrange(totalsubdomains):
        buffh=buffy=0
        if buffersize>0:
            buffh=min(layer.h,currenty+layerblockheight+buffersize)
            buffy=currentwithbuffy
        else:
            buffh=layerblockheight
            buffy=currenty
        subdomain=Subdomain(currenty,layer.x,layerblockheight,layer.w,layer.title+" subdomain "+str(subdindex))
        subdomain.buffx=layer.x
        subdomain.buffw=layer.w
        subdomain.buffh=buffh
        subdomain.buffy=buffy
        pointlist=[]
        for point in layer.get_pointlist():
            if subdomain.isinsidebounds(point,usehalo=True):
                pointlist.append(point.copy())
        #if serial execution then subdomains will need only ordinary list or else a multiprocessing list implementation
        if PCMLConfig.exectype==ExecutorType.serialpython:
            subdomain.set_pointlist(pointlist)
        else:
            subdomain.set_pointlist(pointlist,ref=True)
        subdomainlist.append(subdomain)
        currenty=currenty+layerblockheight
        currentwithbuffy=max(currenty-buffersize,layer.y)
    return subdomainlist

#Create a point subdomain by using raster layer as model from point layer
def pointsubdomainsfromrastersubdomains(pointlayer,rasterlayer,buffersize):
    subdomainlist = []
    rowspersubdomain = float(PCMLConfig.decomposition_granularity)
    numsubdomains = int(math.ceil(float(rasterlayer.nrows)/float(rowspersubdomain)))
    for sdind in xrange(numsubdomains):
        r = rowspersubdomain*sdind
        nrows = rowspersubdomain
        hwithoutbuff=min(rasterlayer.nrows-r,nrows)*rasterlayer.cellsize
        ywithoutbuff=rasterlayer.y + r * rasterlayer.cellsize
        if buffersize>0:
           new_r=max(0,r-buffersize)
           new_h=min(rasterlayer.nrows,r+nrows+buffersize)
           nrows=new_h-new_r
           r=new_r
        else:
            nrows=min(rasterlayer.nrows-r,nrows)
        y = rasterlayer.y + r * rasterlayer.cellsize
        h = nrows * rasterlayer.cellsize
        x = rasterlayer.x
        w = rasterlayer.w
        subdomain = Subdomain(ywithoutbuff, x, hwithoutbuff, w, pointlayer.title+" subdomain "+str(sdind))
        subdomain.buffx=x
        subdomain.buffw=w
        subdomain.buffh=h
        subdomain.buffy=y
        pointlist=[]
        for point in pointlayer.get_pointlist():
            if subdomain.isinsidebounds(point,usehalo=True):
                pointlist.append(point.copy())
        subdomain.set_pointlist(pointlist)
        subdomainlist.append(subdomain)
    return subdomainlist

def pointrasterrowdecomposition(layer,buffersize,layerlist=None):
    if layer.data_structure==Datastructure.array:
        return rowdecomposition(layer,buffersize)
    elif layer.data_structure==Datastructure.pointlist and layerlist is not None:
        return pointsubdomainsfromrastersubdomains(layer,layerlist[1],buffersize)

'''

