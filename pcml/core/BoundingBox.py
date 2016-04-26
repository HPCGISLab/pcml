"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from ..util.SharedMemory import *
from ..util.Messaging import *
from .PCMLPrims import *
import PCMLConfig as PCMLConfig
import copy
import multiprocessing as mp
try:
   PCMLConfig.scipyenabled = 1
   from scipy.spatial import cKDTree
except ImportError as e:
   PCMLConfig.scipyenabled = 0



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
        # Adding buffer data for points
        self.buffx = None
        self.buffy = None
        self.buffh = None
        self.buffw = None

        # Sanity check
        if(h <= 0):
            raise PCMLInvalidInput("BoundingBox does not support a negative or zero height", h)
        if(w <= 0):
            raise PCMLInvalidInput("BoundingBox does not support a negative or zero width", w)

        # Data is held within the internal structure _data
        # and the data structure (e.g., array, list) and type (e.g., location, float, int) must also be described
        self._data = None
        self.data_structure = Datastructure.array  # FIXME: For now we assume the data_structure is an array
        self.data_type = None
        self.tree = None
        # TODO: This should be looked at with dependency on data_structure/type. Perhaps a new class should be created to encapsulate data.
        # TODO: If data_structure is an array, then must describe the cellsize and set a nodata_value
        # By default set to none
        self.nodata_value = None
        self.cellsize = None
        self.nrows = None
        self.ncols = None

    def __repr__(self):
        return "<BoundingBox: (%f,%f) [%f,%f]>" % (self.y, self.x, self.h, self.w)

    def set_nparray(self, nparr, cellsize, nodata_value):
        if nparr is None:
            raise PCMLInvalidInput("BoundingBox.set_nparray does not support a nparr of None", nparr)
        if cellsize is None:
            raise PCMLInvalidInput("BoundingBox.set_nparray does not support a cellsize of None", cellsize)
        if cellsize <= 0:
            raise PCMLInvalidInput("BoundingBox.set_nparray does not support cellsize<=0", cellsize)

        self.data_structure = Datastructure.array
        PCMLNotImplemented("self.data_type is not set")
        # Unfortunately you cannot hide shared memory in another module as rawarray will return the same reference
        # It may have to do with saving the __shmem_data variable, could explore this later
        self.__shmem_data = mp.RawArray(ctypes.c_double, nparr.size)
        self._data = shmem_as_ndarray(self.__shmem_data).reshape(nparr.shape)
        self._data[:, :] = nparr
        self.cellsize = cellsize
        self.nodata_value = nodata_value
        self._reset_dim()

    def get_nparray(self):
        return self._data

    def set_pointlist(self, pointlist, ref=False):
        self.data_structure = Datastructure.pointlist
        # FIXME: Should check if pointlist is a list datastructure
        if not ref:
            self._data = pointlist
        else:
            self._data = mp.Manager().list(pointlist)

    def get_pointlist(self):
        assert(self.data_structure == Datastructure.pointlist), "Cannot get point list if datastructure is not a point list"
        return self._data

    # Bounding box check
    def isinsidebounds(self, point, usehalo=False):
        x, y, w, h = self.x, self.y, self.w, self.h
        if usehalo:
            x, y, w, h = self.buffx, self.buffy, self.buffw, self.buffh
        a = y
        b = y + h
        c = point['y']
        if point['x'] < x+w and point['x'] >= x and point['y'] < y+h and point['y'] >= y:
            return True
        else:
            return False

    # Get points with out including halozone
    def get_pointlistwithouthalozone(self):
        pointswithouthalozone = []
        for point in self._data:
            if self.isinsidebounds(point):
                pointswithouthalozone.append(point)
        return pointswithouthalozone

    def _reset_dim(self):
        nparr = self.get_nparray()
        # Set nrows and ncols based on dimensions of array
        self.nrows = len(nparr)
        self.ncols = len(nparr[0])

        h = self.nrows * self.cellsize
        w = self.ncols * self.cellsize
        if h != self.h:
            PCMLUserInformation("Updating height from " + str(self.h) + " to " + str(h))
            self.h = h
        if w != self.w:
            PCMLUserInformation("Updating width from " + str(self.w) + " to " + str(w))
            self.w = w

    def set_data_ref(self, ref):
        """
        Set the _data variable to a reference (used for shared memory accesses - particularly subdomains)
        """
        self._data = ref
        self._reset_dim()

    def get_locval(self, loc):
        newloc = copy.copy(loc)
        newloc['v'] = self._data[loc['r'] - self.r][loc['c'] - self.c]
        return newloc

    def get_ind_from_loc(self, loc):
        return {
                        'r': loc['r'] - self.r,
                        'c': loc['c'] - self.c}

    def get_yxloc(self, locind):
        return {
                'y': self.y + self.cellsize * locind['r'],
                'x': self.x + self.cellsize * locind['c'],
                'z': self._data[locind['r']][locind['c']]}

    def bufferedlocgetarr(self, loc, buffersize):
        ind = self.get_ind_from_loc(loc)
        # FIXME: NEED TO DOUBLE CHECK THIS LOGIC, NEED TO RECORD RDIFF,CDIFF FOR H/W
        r = max(0, ind['r'] - buffersize)
        c = max(0, ind['c'] - buffersize)
        h = buffersize + (ind['r'] - r) + 1
        if (r + h > self.nrows):
            h = self.nrows - r
        w = buffersize + (ind['c'] - c) + 1
        if (c + w > self.ncols):
            w = self.ncols - c
        return self.slice_nparray(r, c, h, w)

    def slice_nparray(self, r, c, h, w):
        # IMPORTANT: Notice slice uses r,c NOT x,y!
        return self._data[r:r + h, c:c + w]

    def print_data(self):
        """Print out all data."""
        print(self._data)

    # Get neighbors for points using cKDTree
    def getneighbors(self, location, count=1, radius=np.inf, excludesearchlocation=False, distreq=False):
        if PCMLConfig.scipyenabled == 0:
            PCMLNotSupported("SciPy module required for getneighbors()")
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
        dist, neighbors = self.tree.query(points, k=count, distance_upper_bound=radius)
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
