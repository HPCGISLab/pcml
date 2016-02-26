# Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory.
# All rights reserved. Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file. Contributors and authors
# can be found in the AUTHORS file and we always welcome new contributors.

from .BoundingBox import *
from ..lib.LocalOperationPrimitives import *
from .Decomposition import *
from ..util.Messaging import *
from .PCMLPrims import *


class Layer(BoundingBox):
    def __init__(self, y, x, h, w, title):
        """
        Layer class a map layer.
        Args:
            :param y (double):  The y location (lower left corner typically)
                                of the :class:`Layer`.
            :param x (double):  The x location (lower left corner typically)
                                of the :class:`Layer`.
            :param h (double):  The height of the :class:`Layer`.
            :param w (double):  The width of the :class:`Layer`.
            :param title (str): The title of the :class:`Layer`
        """
        super(Layer, self).__init__(y, x, h, w)
        self.title = title

    def __repr__(self):
        """
        Print layer infromation.
        """
        return "<Layer: (%f,%f) [%f,%f] : %s>" % \
               (self.y, self.x, self.h, self.w, self.title)

    def duplicate(self):
        """
        Create a new layer, based on this layer
        :returns: a new layer
        """
        # FIXME: A function such as setfromlayer() should be defined
        #        to do this automatically
        # WARNING: This does not get all values
        newlayer = Layer(self.y, self.x, self.h, self.w,
                         self.title + " (duplicate)")
        if self.data_structure == Datastructure.array:
            newlayer.set_nparray(np.zeros((self.nrows, self.ncols)),
                                 self.cellsize, self.nodata_value)
        elif self.data_structure == Datastructure.pointlist:
            newlayer.set_pointlist(self.get_pointlist())
        # TODO: PCMLTODO("Double check that all of the values are copied over")
        return newlayer

    def tosubdomains(self,boundingboxlist):
        '''
        Given a list of bounding boxes, return a list of subdomains from this layer
        '''

        subdomainlist = []

        # Loop over bounding boxes
        for bb_index in xrange(len(boundingboxlist)):
            bb = boundingboxlist[bb_index]

            # FIXME: This is done in other locations so this code should be removed
            # If the buffer is not disabled, then replace coordinates with buffer coordinates
            if False:
                bb.x=bb.x_buf
                bb.y=bb.y_buf
                bb.h=bb.h_buf
                bb.w=bb.w_buf
                bb.nrows=bb.nrows_buf
                bb.ncols=bb.ncols_buf
                bb.r=bb.r_buf
                bb.c=bb.c_buf

            #subdomain = Subdomain(bb.y,bb.x,bb.h,bb.w,"Subdomain of " + self.title)
            subdomain = Subdomain(bb.y, bb.x, bb.h, bb.w, self.title+" subdomain "+str(bb_index))

            subdomain.cellsize=self.cellsize
            subdomain.nodata_value=self.nodata_value
            subdomain.r=bb.r
            subdomain.c=bb.c
            subdomain.nrows=bb.nrows
            subdomain.ncols=bb.ncols

            # FIXME: Eventually this will be removed
            subdomain.r_buf=bb.r_buf
            subdomain.c_buf=bb.c_buf
            subdomain.h_buf=bb.h_buf
            subdomain.w_buf=bb.w_buf
            subdomain.nrows_buf=bb.nrows_buf
            subdomain.ncols_buf=bb.ncols_buf
            subdomain.y_buf=bb.y_buf
            subdomain.x_buf=bb.x_buf

            if self.data_structure == Datastructure.array: # Copy array data to subdomain
                # Extract an array slice (reference to data in a layer for lower memory overhead)
                # from the layer and set the data reference for the subdomain to use
                # FIXME: Return perhaps arrslice=self.slice_nparray(bb.r_buf,0,bb.nrows_buf,bb.ncols_buf)
                arrslice=self.slice_nparray(bb.r,0,bb.nrows,bb.ncols)
                subdomain.set_data_ref(arrslice)
            elif self.data_structure == Datastructure.pointlist: # Copy points to subdomain
                # FIXME: This should be abstracted out in the future
                pointlist=[]
                for point in self.get_pointlist():
                    #if subdomain.isinsidebounds(point,usehalo=True):
                    if subdomain.isinsidebounds(point):
                        pointlist.append(point.copy())
                subdomain.set_pointlist(pointlist)

            subdomainlist.append(subdomain)

        return subdomainlist

    # Syntax Sugar

    def __add__(self, right):
        """
        Enables layer1 + layer2 functionality using LocalSum
        """
        print('__add__ %s + %s' % (self.title, right.title))
        return LocalSum(self, right)

    def __mul__(self, right):
        """
        Enables layer1 * layer2 functionality using LocalMult
        """
        print('__mul__ %s * %s' % (self.title, right.title))
        return LocalMult(self, right)

    def __div__(self, right):
        """
        Enables layer1 / layer2 functionality using LocalDivision
        """
        print('__div__ %s / %s' % (self.title, right.title))
        return LocalDivision(self, right)

    def __truediv__(self, right):
        """
        Enables layer1 / layer2 functionality using LocalDivision
        """
        print('__div__ %s / %s' % (self.title, right.title))
        return LocalDivision(self, right)

    def __sub__(self, right):
        """
        Enables layer1 - layer2 functionality using LocalSubtraction
        """
        print('__sub__ %s - %s' % (self.title, right.title))
        return LocalSubtraction(self, right)
