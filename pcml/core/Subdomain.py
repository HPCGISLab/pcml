# Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory.
# All rights reserved. Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file. Contributors and authors
# can be found in the AUTHORS file and we always welcome new contributors.

from .BoundingBox import *
import itertools
import PCMLConfig
from .PCMLPrims import *
import multiprocessing as mp


class Subdomain(BoundingBox):
    """
    Subdomain class represents a subdomain (portion) of a layer.
    """

    def __init__(self, y, x, h, w, title):
        """
        Create a new subdomain object.
        :class:`Subdomain` objects are generators,
               iterating over a subdomain object will give you all locations
               in the subdomain.
        Args:
            :param y (double):  The y location (lower left corner typically)
                                of the :class:`BoundingBox`.
            :param x (double):  The x location (lower left corner typically)
                                of the :class:`BoundingBox`.
            :param h (double):  The height of the :class:`BoundingBox`.
            :param w (double):  The width of the :class:`BoundingBox`.
            :param title (str): title(name) of the portion
        """
        super(Subdomain, self).__init__(y, x, h, w)
        self.title = title
        self.layer = None # By default do not set a layer relationship

        # By default do not set iteration_number # FIXME: This should be removed.
        #self.iteration_number = None

        # Only valid when data_structure==Datastructure.array
        #self.r = None
        #self.c = None
        # iteration_number is used to keep track of re-processing subdomains
        # this is only necessary for a small number of operations 
        if PCMLConfig.exectype == ExecutorType.serialpython:
            self._iteration_number = 0
        elif PCMLConfig.exectype == ExecutorType.parallelpythonqueue:
            self._iteration_number = mp.Value('i', 0)
            print "iteration number is",self.iteration_number
        else:
            PCMLNotSupported("This type of executor is not supported")

    def get_iteration_number(self):
        """
        Get the number of iterations for this subdomain.
        """
        if PCMLConfig.exectype == ExecutorType.serialpython:
            return self._iteration_number
        else:
            return self._iteration_number.value

    def set_iteration_number(self, val):
        """
        Set the number of iterations for this subdomain.
        """
        if PCMLConfig.exectype == ExecutorType.serialpython:
            self._iteration_number = val
        else:
            self._iteration_number.value = val

    def __repr__(self):
        """
        Print Subdomain information
        """
        return "<Subdomain: (%f,%f) [%f,%f] : %s>" % (self.y, self.x, self.h,
                                                      self.w, self.title)

    iteration_number = property(get_iteration_number,set_iteration_number) 


