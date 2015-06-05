"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from .BoundingBox import *
import itertools
import PCMLConfig
from .PCMLPrims import *
import multiprocessing as mp
class Subdomain(BoundingBox):
    """Subdomain class represents a subdomain (portion) of a layer."""

    def __init__(self, y, x, h, w, title):
        """Create a new subdomain object.
        :class:`Subdomain` objects are generators, iterating over a subdomain object will give you all locations in the subdomain.
        Args:
            :param y (double): The y location (lower left corner typically) of the :class:`BoundingBox`.
            :param x (double): The x location (lower left corner typically) of the :class:`BoundingBox`. 
            :param h (double): The height of the :class:`BoundingBox`. 
            :param w (double): The width of the :class:`BoundingBox`. 
            :param title (str): title(name) of the portion
        """
        super(Subdomain, self).__init__(y, x, h, w)
        self.title=title

        # Only valid when data_structure==Datastructure.array
        self.r=None
        self.c=None
        #iter count is used to check whether the processing for subdomain is done at scheduler level. If it is 0 then the subdomain won't be processed any more
        if PCMLConfig.exectype==ExecutorType.serialpython:
            self.itercount=0
        elif PCMLConfig.exectype==ExecutorType.parallelpythonqueue:
            self.itercount=mp.Value('i',0)
    #getter for iter count
    def get_itercount(self):
        if PCMLConfig.exectype==ExecutorType.serialpython:
            return self.itercount
        else:
            return self.itercount.value
    def set_itercount(self,val):
        if PCMLConfig.exectype==ExecutorType.serialpython:
            self.itercount=val
        else:
            self.itercount.value=val

    def __repr__(self):
        return "<Subdomain: (%f,%f) [%f,%f] : %s>" % (self.y,self.x,self.h,self.w,self.title)

