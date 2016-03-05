"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from .BoundingBox import *
import itertools


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
        self.title = title

        # Only valid when data_structure==Datastructure.array
        self.r = None
        self.c = None

    def __repr__(self):
        return "<Subdomain: (%f,%f) [%f,%f] : %s>" % (self.y, self.x, self.h, self.w, self.title)
