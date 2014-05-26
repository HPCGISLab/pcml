"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from ..util.SharedMemory import *
from ..util.Messaging import *
import copy

class Datastructure():
    """ Enumeration class. Defines the type of data structure.
    """
    array = 1
    pointlist = 2

class BoundingBox(object):
    """
    BoundingBox defines a rectangular location (y,x) + (h,w) and may contain data describing something within its boundaries.
    This class is also used as the parent class of Layer and Subdomain.
    """

    def __init__(self, y, x, h, w):
        """Create a new BoundingBox.
            :param y (double): The y location (lower left corner typically) of the :class:`BoundingBox`.
            :param x (double): The x location (lower left corner typically) of the :class:`BoundingBox`.
            :param h (double): The height of the :class:`BoundingBox`.
            :param w (double): The width of the :class:`BoundingBox`.
        """
        self.y = y
        self.x = x
        self.h = h
        self.w = w

        # Data is held within the internal structure _data
        # and the data structure (e.g., array, list) and type (e.g., location, float, int) must also be described
        self._data = None
        self.data_structure = Datastructure.array # FIXME: For now we assume the data_structure is an array
        self.data_type = None

        # TODO: This should be looked at with dependency on data_structure/type. Perhaps a new class should be created to encapsulate data.
        # TODO: If data_structure is an array, then must describe the cellsize and set a nodata_value
        # By default set to none
        self.nodata_value=None
        self.cellsize=None
        self.nrows=None
        self.ncols=None

    def __repr__(self):
        return "<BoundingBox: (%f,%f) [%f,%f]>" % (self.y, self.x, self.h, self.w)

    def _reset_dim(self):
        nparr=self.get_nparray()
        # Set nrows and ncols based on dimensions of array
        self.nrows=len(nparr)
        self.ncols=len(nparr[0])

        h=self.nrows*self.cellsize
        w=self.ncols*self.cellsize
        if h != self.h:
            PCMLUserInformation("Updating height from "+str(self.h)+" to "+str(h))
            self.h=h
        if w != self.w:
            PCMLUserInformation("Updating width from "+str(self.w)+" to "+str(w))
            self.w=w

    def set_nparray(self,nparr,cellsize,nodata_value):
        self.data_structure=Datastructure.array
        PCMLNotImplemented("self.data_type is not set")
        # Unfortunately you cannot hide shared memory in another module as rawarray will return the same reference
        # It may have to do with saving the __shmem_data variable, could explore this later
        self.__shmem_data=mp.RawArray(ctypes.c_double,nparr.size)
        self._data=shmem_as_ndarray(self.__shmem_data).reshape(nparr.shape)
        self._data[:,:]=nparr
        self.cellsize=cellsize
        self.nodata_value=nodata_value
        self._reset_dim()

    def set_pointlist(self,pointlist):
        self.data_structure=Datastructure.pointlist
        # FIXME: Should check if pointlist is a list datastructure
        self._data=pointlist

    def get_pointlist(self):
        assert(self.data_structure==Datastructure.pointlist)
        return self._data

    def get_nparray(self):
        return self._data

    def set_data_ref(self,ref):
        """
        Set the _data variable to a reference (used for shared memory accesses - particularly subdomains)
        """
        self._data=ref
        self._reset_dim()

    def get_locval(self,loc):
        newloc = copy.copy(loc)
        newloc['v']=self._data[loc['r']-self.r][loc['c']-self.c]
        return newloc
    def get_ind_from_loc(self,loc):
        return { \
                        'r':loc['r']-self.r, \
                        'c':loc['c']-self.c}
    def get_yxloc(self,locind):
        return { \
                'y':self.y+self.cellsize*locind['r'], \
                'x':self.x+self.cellsize*locind['c'], \
                'z':self._data[locind['r']][locind['c']]}

    def bufferedlocgetarr(self,loc,buffersize):
        ind=self.get_ind_from_loc(loc)
        # FIXME: NEED TO DOUBLE CHECK THIS LOGIC, NEED TO RECORD RDIFF,CDIFF FOR H/W
        r=max(0,ind['r']-buffersize)
        c=max(0,ind['c']-buffersize)
        h=buffersize+(ind['r']-r)+1
        if (r+h > self.nrows):
            h=self.nrows-r
        w=buffersize+(ind['c']-c)+1
        if (c+w > self.ncols):
            w=self.ncols-c
        return self.slice_nparray(r,c,h,w)

    def slice_nparray(self,r,c,h,w):
        # IMPORTANT: Notice slice uses r,c NOT x,y!
        return self._data[r:r + h, c:c + w]

    def print_data(self):
        """Print out all data."""
        if(self.data_structure==Datastructure.array):
            print self._data
        else:
            raise PCMLNotImplemented("print_data other data_structured are not supported")
