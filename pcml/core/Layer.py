"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from .BoundingBox import *
from ..lib.LocalOperationPrimitives import *
from .Decomposition import *
from ..util.Messaging import *
from .PCMLPrims import *
class Layer(BoundingBox):
    def __init__(self,y,x,h,w,title):
        """Layer class a map layer.
        Args:
            :param y (double): The y location (lower left corner typically) of the :class:`Layer`.
            :param x (double): The x location (lower left corner typically) of the :class:`Layer`. 
            :param h (double): The height of the :class:`Layer`. 
            :param w (double): The width of the :class:`Layer`. 
            :param title (str): The title of the :class:`Layer`
        """
        super(Layer, self).__init__(y, x, h, w)
        self.title = title

    def __repr__(self):
        return "<Layer: (%f,%f) [%f,%f] : %s>" % (self.y,self.x,self.h,self.w,self.title)

    def duplicate(self):
        """
        Create a new layer, based on this layer
        :returns: a new layer
        """
        # FIXME: A function such as setfromlayer() should be defined that will do this automatically
        newlayer=Layer(self.y,self.x,self.h,self.w,self.title+" (duplicate)")
        if self.data_structure == Datastructure.array:
            newlayer.set_nparray(np.zeros((self.nrows,self.ncols)),self.cellsize,self.nodata_value)
        elif self.data_structure == Datastructure.pointlist:
            newlayer.set_pointlist(self.get_pointlist())
        #TODO: PCMLTODO("Double check that all of the values are copied over")
        return newlayer

    '''
    def decomposition(self, method, buffersize):
        """Return a list of subdomains based on decomposition method.
        """
        if method == DecompositionMethod.row:
            return rowdecomposition(self, buffersize)
        else:
            raise PCMLNotImplemented("Other decomposition methods are not implemented yet")
    '''

    ### Syntax Sugar

    def __add__(self,right):
        print('__add__ %s + %s' %(self.title,right.title))
        return LocalSum(self,right)

    def __mul__(self,right):
        print('__mul__ %s * %s' %(self.title,right.title))
        return LocalMult(self,right)
    
    def __div__(self,right):
        print('__div__ %s / %s' %(self.title,right.title))
        return LocalDivision(self,right)

    def __truediv__(self,right):
        print('__div__ %s / %s' %(self.title,right.title))
        return LocalDivision(self,right)

    def __sub__(self,right):
        print('__sub__ %s - %s' %(self.title,right.title))
        return LocalSubtraction(self,right)



