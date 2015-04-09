"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from ..util.Messaging import *
from .Decomposition import *
from .BoundingBox import *
from abc import ABCMeta, abstractmethod

class OpClass():
    """ Enumeration class. Classifies operations as local, focal, zonal, or global.
    """
    localclass = 1
    focalclass = 2
    zonalclass = 3
    globalclass = 4


class Operation(object):
    __metaclass__ = ABCMeta
    def __init__(self,name,*args,**kwargs):
        """Operations are applied to layers.
        Args:
            :param name (string): String representation of Operation name
            :param layers (tuple): Tuple of layers to apply operation on
            :param opclass (OpClass): Operation classification (local, focal, zonal, global)
        """

        # Derive class name from operation name
        self.name = name

        PCMLTODO("Only row decomposition method supported, so hardcoding for now")
        #self.decomposition_method=DecompositionMethod.row

        _layerstuple = kwargs.get('layers', None)
        if _layerstuple!=None:
            self._layers = list(_layerstuple)
        self.opclass = kwargs.get('opclass', OpClass.localclass)
        self.buffersize = kwargs.get('buffersize', 0)
        self.decomposition = kwargs.get('decomposition',rowdecomposition) # By default use row decomposition
        self.iterator = kwargs.get('iterator', 'rowiterator')
        if self.opclass==OpClass.localclass and self.buffersize != 0:
            raise PCMLOperationError("Buffersize should be 0 for localclass currently %s" % self.buffersize)
	#If zonal operation we want the entire layer data
	if self.opclass==OpClass.zonalclass:
            self.buffersize=999999999999

    def __repr__(self):
        return "<Operation: %s : %i layers>" % (self.name,len(self._layers))

    def getOutputLayers(self):
        PCMLTODO("Need to support more than one output layer")
        return self._layers[0]

    def _decompositioninit(self):
        # Duplicate a layer to create an output layer with the correct dimensions
        # Get the first layer
        firstlayer=self._layers[0]

        # The output layer is a duplicate of the first layer
        outputlayer=firstlayer.duplicate()
        outputlayer.title="Output for operation %s"%self.name

        self._layers.insert(0, outputlayer)  # Add the output layer to the front of the layers list

    # By default we use rowdecomposition as our decomposition method
    # Users may override decomposition with any other method they would like
    #def decomposition(self,layer,buffersize):
    #    return rowdecomposition(layer,buffersize)

    def _decompositionrun(self):
        """ Divides the :member:_layers into subdomains for further processing.
        The decomposition method is defined by :member:`decompositionmethod`.
        You can also define you own decomposition algorithm by overriding this method.
        """

        PCMLTODO("Need to support multiple output layers, this can be done by overriding decomposition and inserting multiple output layers")
        listofsubdomains = []

        self._decompositioninit()

        # The output layer is the first layer in the layers list (self.layers[0])
        # Decompose it with a 0 buffer
        #listofsubdomains.append(self._layers[0].decomposition(self.decomposition_method, 0))
        listofsubdomains.append(self.decomposition(self._layers[0], 0))

        for layer in self._layers:
            if layer != self._layers[0]: # Skip the output layer, because it was already decomposed and added
               #listofsubdomains.append(layer.decomposition(self.decomposition_method, self.buffersize)) # buffer size is set based on classification (L,F,Z,G)
               listofsubdomains.append(self.decomposition(layer, self.buffersize)) # buffer size is set based on classification (L,F,Z,G)

        # The listofsubdomains is inverted using zip and map to create a list of lists
        # so that each subdomain is grouped with the corresponding subdomain from each layer (see example below)
        subdomainlists = map(list,zip(*listofsubdomains))

        # listofsubdomains = ( (layer1subdomain1 , layer1subdomain2) , (layer2subdomain1 , layer2subdomain2) )
        # subdomainlists =   ( (layer1subdomain1 , layer2subdomain1) , (layer1subdomain2 , layer2subdomain2) )

        return subdomainlists

    def executor(self,subdomains):
        """ Executor handles processing of the function by iterating over locations in a subdomain
        :return: #TODO: Undefined return value.
        """
        PCMLTODO("executor assumes single subdomain as output, which is not universal for all operations")

        outsubdomain = subdomains.pop(0)
        outarr = outsubdomain.get_nparray()
        outsubdomain.set_iterator(self.iterator)
        if outsubdomain.data_structure!=Datastructure.array:
           print "datatype",outsubdomain.data_type,"arraydt",Datastructure.array
           PCMLNotSupported("Executor currently assumes an array data structure")

        PCMLTODO("Sanity check subdomains are all the same dimensions")

        # Iterate over locations in the outsubdomain and apply function to each location
        #for locind in outsubdomain:
        for loc in outsubdomain:
            l = [] # Create an empty list to store locations
            for sd in subdomains:
                if sd.data_structure!=Datastructure.array: # Skip non array subdomains
                    continue
                # Get a location in this subdomain with same coordinates as locind
                locv=sd.get_locval(loc)
                l.append(locv) # append to list of locations
            val = self.function(l,subdomains) # Apply function to all locations
            outarr[loc['r']-outsubdomain.r][loc['c']-outsubdomain.c]=val # Set val to outarr at locind

    def function(self,locations,subdomains):
        raise PCMLOperationError("Operation function is not implemented")

