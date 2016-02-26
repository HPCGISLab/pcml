# Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory.
# All rights reserved. Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file. Contributors and authors
# can be found in the AUTHORS file and we always welcome new contributors.

"""
A space-time bounding box (box) is a three-dimensional rectangular box that
bounds both spatial (X and Y) and temporal (T) dimensions.
A box may contain spatial, temporal, or spatio-temporal data in many formats.
A box can be reduced by one dimension to create a layer (spatial or temporal)
or by two dimensions to create a series.
A single series (space-time point) is a special case of a series.
"""

from ..util.SharedMemory import *
from ..util.Messaging import *
from .PCMLPrims import *
import copy
import multiprocessing as mp
from scipy.spatial import cKDTree


class BoundingBox(object):
    """
    BoundingBox defines a three-dimensional rectangular box (y,x,t) + (h,w,d)
    and may contain data describing something within its boundaries.
    This class is also used as the parent class of Layer and Subdomain.
    Currently, the temporal dimension time (t) and duration (d) are set to 0
    by default, but in future versions these parameters may be required.
    """

    def __init__(self, y, x, h, w, t=0, d=0):
        """
        Create a new :class:`BoundingBox`.
        :param y (double): The y location (lower left corner typically).
        :param x (double): The x location (lower left corner typically).
        :param h (double): The height of the :class:`BoundingBox`.
        :param w (double): The width of the :class:`BoundingBox`.
        :param t (double): The start time of the :class:`BoundingBox`.
        :param d (double): The duration of the :class:`BoundingBox`.
        """

        self.y = y
        self.x = x
        self.h = h
        self.w = w
        self.t = t
        self.d = d

        # Set optional buffer data to none by default 
        self.x_buf = None
        self.y_buf = None
        self.h_buf = None
        self.w_buf = None

        self.r_buf = None
        self.c_buf = None
        self.nrows_buf = None
        self.ncols_buf = None

        # Sanity check
        if(h <= 0):
            raise PCMLInvalidInput(
                    "BoundingBox does not support a negative or zero height",
                    h)
        if(w <= 0):
            raise PCMLInvalidInput(
                    "BoundingBox does not support a negative or zero width",
                    w)

        # Data is held within the internal structure _data
        # and the data structure (e.g., array, list)
        # and type (e.g., location, float, int) must also be described
        self._data = None
        # FIXME: For now we assume the data_structure is an array
        self.data_structure = Datastructure.array
        self.data_type = None
        self.tree = None
        # TODO: This should be looked at with dependency on
        # data_structure/type. Perhaps a new class should be created
        # to encapsulate data.
        # TODO: If data_structure is an array, then must describe
        # the cellsize and set a nodata_value

        # By default set to none
        self.nodata_value = None
        self.cellsize = None
        self.nrows = None
        self.ncols = None
        self.r = None
        self.c = None

    def __repr__(self):
        return "<BoundingBox: (%f,%f) [%f,%f]>" % \
               (self.y, self.x, self.h, self.w)

    def set_nparray(self, nparr, cellsize, nodata_value):
        """
        Takes a Numpy array, size of a raster cell, and a nodata value
        as input, which is then saved as part of the BoundingBox
        """
        if nparr is None:
            raise PCMLInvalidInput("BoundingBox.set_nparray \
                                    does not support a nparr of None",
                                   nparr)
        if cellsize is None:
            raise PCMLInvalidInput("BoundingBox.set_nparray \
                                    does not support a cellsize of None",
                                   cellsize)
        if cellsize <= 0:
            raise PCMLInvalidInput("BoundingBox.set_nparray \
                                    does not support cellsize<=0",
                                   cellsize)

        self.data_structure = Datastructure.array
        PCMLNotImplemented("self.data_type is not set")
        # Unfortunately you cannot hide shared memory in another module
        # as rawarray will return the same reference
        # It may have to do with saving the __shmem_data variable,
        # could explore this later
        self.__shmem_data = mp.RawArray(ctypes.c_double, nparr.size)
        self._data = shmem_as_ndarray(self.__shmem_data).reshape(nparr.shape)
        self._data[:, :] = nparr
        self.cellsize = cellsize
        self.nodata_value = nodata_value
        self._reset_dim()

    def get_nparray(self):
        """
        Returns a nparray of the data.
        """
        return self._data

    def set_pointlist(self, pointlist, ref=False):
        """
        Establish a list of points based on pointlist variable.
        If ref is true, then use shared memory MP.manager for
        a shared memory pointlist (slower performance).
        """
        self.data_structure = Datastructure.pointlist
        # FIXME: Should check if pointlist is a list datastructure
        if not ref:
            self._data = pointlist
        else:
            self._data = mp.Manager().list(pointlist)

    def get_pointlist(self):
        """
        Return list of points that are within the BoundingBox
        including inside the ghost zone.
        """
        assert(self.data_structure == Datastructure.pointlist), "Cannot get \
               point list if datastructure is not a point list"
        return self._data

    # Bounding box check
    def isinsidebounds(self, point, usehalo=False):
        """
        Return boolean whether point is within the bounds of a BoundingBox,
        which may optionally could use a halo or buffer distance.
        """
        x, y, w, h = self.x, self.y, self.w, self.h
        if usehalo:
            x, y, w, h = self.x_buf, self.y_buf, self.w_buf, self.h_buf
        # FIXME: These should be removed.
        #a = y
        #b = y+h
        #c = point['y']
        #if point['x'] <= x+w and point['y'] <= y+h: 
        #if point['x'] >= x and point['y'] >= y:
        #    return True
        #print "py",point['y'],"y",y
        #return False
        if point['x'] < x+w and point['x'] >= x and \
           point['y'] < y+h and point['y'] >= y:
            return True
        else:
            return False

    # FIXME: This should be included in the isinsidebounds
    #        by using an optional parameter
    def get_pointlistwithouthalozone(self):
        """
        Return a list of points that are within the boundaries
        of the BoundingBox that does not include the ghost zone.
        """
        pointswithouthalozone = []
        for point in self._data:
            if self.isinsidebounds(point):
                pointswithouthalozone.append(point)
        return pointswithouthalozone

    def _reset_dim(self):
        """
        Reset the dimensions of BoundingBox based on
        the number of rows and columns in the nparray
        """
        nparr = self.get_nparray()
        self.nrows = len(nparr)
        self.ncols = len(nparr[0])

        h = self.nrows * self.cellsize
        w = self.ncols * self.cellsize

        # Warn user if the dimensions are changing
        if h != self.h:
            PCMLUserInformation("Updating height from " +
                                str(self.h) + " to " + str(h))
            self.h = h
        if w != self.w:
            PCMLUserInformation("Updating width from " +
                                str(self.w) + " to " + str(w))
            self.w = w

    def set_data_ref(self, ref):
        """
        Set the _data variable to a reference
        (used for shared memory accesses - particularly subdomains)
        """
        self._data = ref
        self._reset_dim()

    def get_locval(self, loc):
        """
        Return a new location with 'v' set using r,c location
        """
        newloc = copy.copy(loc)
        newloc['v'] = self._data[loc['r'] - self.r][loc['c'] - self.c]
        return newloc

    def get_ind_from_loc(self, loc):
        """
        Return an index (r,c) from a location (r,c),
        which is calculated by the offset of the subdomain/layer/box
        """
        return {
                        'r': loc['r'] - self.r,
                        'c': loc['c'] - self.c}

    def get_yxloc(self, locind):
        """
        Return a x,y location with z value when given
        a locind in r,c format.
        """
        return {
                'y': self.y + self.cellsize * locind['r'],
                'x': self.x + self.cellsize * locind['c'],
                'z': self._data[locind['r']][locind['c']]}

    def bufferedlocgetarr(self, loc, buffersize):
        """
        Get a slice of a nparray surrounding a location (loc)
        with number of cells on each size defined by buffersize
        """
        ind = self.get_ind_from_loc(loc)
        # FIXME: DOUBLE CHECK THIS LOGIC, NEED TO RECORD RDIFF,CDIFF FOR H/W
        r = max(0, ind['r'] - buffersize)
        c = max(0, ind['c'] - buffersize)
        h = buffersize + (ind['r'] - r) + 1
        if (r + h > self.nrows):
            h = self.nrows - r
        w = buffersize + (ind['c'] - c) + 1
        if (c+w > self.ncols):
            w = self.ncols - c
        return self.slice_nparray(r, c, h, w)

    def copy_boundingbox_data(self, bbox):
        """
        Will accept a bounding box (bbox) and will copy data from the bbox
        to this bounding box. At least what fits within the boundaries of
        the boundingbox and self.
        """

        # Currently only array structures are supported,
        # in the future points may be supported also
        assert(self.data_structure == Datastructure.array)

        # Find the largest lower left coordinate (r,c) that are shared
        # This coordinate is relative to self
        llr_index = max(0, bbox.r - self.r) # Make sure it stays non-negative
        llc_index = max(0, bbox.c - self.c) # Make sure it stays non-negative

        # Find the largest lower left coordinate (r,c) that are shared
        # This coordinate is relative to bbox 
        llr_index_bbox = max(0, self.r - bbox.r) # Make sure it stays non-negative
        llc_index_bbox = max(0, self.c - bbox.c) # Make sure it stays non-negative

        # Find the smallest upper right coordinate (r,c)
        ur_coord_r = min(bbox.r + bbox.nrows, self.r + self.nrows)
        ur_coord_c = min(bbox.c + bbox.ncols, self.c + self.ncols)

        # Calculate the height/width of the shared bounding box
        nrows = (ur_coord_r - self.r) - llr_index
        ncols = (ur_coord_c - self.c) - llc_index

        bbox_slice = bbox.slice_nparray(llr_index_bbox, llc_index_bbox, nrows, ncols)
        self._data[llr_index:llr_index + nrows, llc_index:llc_index + ncols] = bbox_slice

    def slice_nparray(self, r, c, h, w):
        """
        Return a slice (or rectangular portion) of the nparray of data.
        It is important to note that slice uses r,c and not x,y
        """
        return self._data[r:r + h, c:c + w]

    def print_data(self):
        """
        Print out all data.
        """
        print(self._data)

    def getneighbors(self, location, count=1, radius=np.inf,
                     excludesearchlocation=False, distreq=False):
        """
        Get neighbors for points using cKDTree
        within radius (if radius is not provided then return all points
        """
        if excludesearchlocation:
            count += 1
        if self.tree is None:
            pointlist = self.get_pointlist()
            pointdata = []
            if len(pointlist) == 0:
                if not distreq:
                    return []
                else:
                    return [[], []]
            for point in pointlist:
                pointdata.append([point['x'], point['y']])
            self.tree = cKDTree(pointdata)
        points = [location['x'], location['y']]
        dist, neighbors = self.tree.query(points,
                                          k=count,
                                          distance_upper_bound=radius)
        if count == 1:
            if neighbors == self.tree.n:
                if not distreq:
                    return []
                else:
                    return [[], []]
            else:
                if not distreq:
                    return [neighbors]
                else:
                    return [[neighbors], [dist]]
        if not distreq:
            return neighbors[neighbors != self.tree.n]
        else:
            return [neighbors[neighbors != self.tree.n], dist[dist != np.inf]]
