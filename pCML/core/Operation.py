"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from ..util.Messaging import *
from .LayerDecomposition import *
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
        self.decomposition_method=DecompositionMethod.row

        _layerstuple = kwargs.get('layers', None)
        if _layerstuple!=None:
            self._layers = list(_layerstuple)
        self.opclass = kwargs.get('opclass', OpClass.localclass)
        self.buffersize = kwargs.get('buffersize', 0)

    def __repr__(self):
        return "<Operation: %s : %i layers>" % (self.name,len(self._layers))

    def getOutputLayers(self):
        PCMLTODO("Need to support more than one output layer")
        return self._layers[0]



    def decomposition(self):
        """ Divides the :member:_layers into subdomains for further processing.
        The decomposition method is defined by :member:`decompositionmethod`.
        You can also define you own decomposition algorithm by overriding this method.
        """

        listofsubdomains = []

        PCMLTODO("Need to support multiple output layers")

        # Duplicate a layer to create an output layer with the correct dimensions
        # Get the first layer
        firstlayer=self._layers[0]

        # The output layer is a duplicate of the first layer
        outputlayer=firstlayer.duplicate()
        outputlayer.title="Output for operation %s"%self.name

        # IMPORTANT : first entry in the listoflayerportions is the output portions list!
        if self.opclass == OpClass.localclass:  # If a local operation
            self._layers.insert(0, outputlayer)  # Add the output layer to the front of the layers list
            for layer in self._layers:
                listofsubdomains.append(layer.decomposition(self.decomposition_method, 0)) # 0 is the buffer size, which is 0 for local operations
        elif self.opclass == OpClass.focalclass:  # If a focal operation
            # For non-local operations, we don't want a buffer for our output layer otherwise locations may overlap
            PCMLTODO("Need to ensure that buffer versus non-buffer decompositions are identical for this to work")
            listofsubdomains.append(outputlayer.decomposition(self.decomposition_method, 0))
            for layer in self._layers:
                listofsubdomains.append(layer.decomposition(self.decomposition_method, self.buffersize))
            self._layers.insert(0, outputlayer)  # Add the output layer to the front of the layers list after loop
                                                 # so it is not divided with a buffersize (in previous for loop)
        elif self.opclass == OpClass.globalclass:  # If a focal operation
            # For non-local operations, we don't want a buffer for our output layer otherwise locations may overlap
            PCMLTODO("Need to ensure that buffer versus non-buffer decompositions are identical for this to work")
            listofsubdomains.append(outputlayer.decomposition(self.decomposition_method, 0))
            for layer in self._layers:
                listofsubdomains.append(layer.decomposition(self.decomposition_method, -1))  # -1 means the buffer is infinite sized
            self._layers.insert(0, outputlayer)  # Add the output layer to the front of the layers list after loop
                                                 # so it is not divided with a buffersize (in previous for loop)
        else:
            raise PCMLClassificationNotSupported(self.opclass)

        # listofsubdomains = ( (layer1subdomain1 , layer1subdomain2) , (layer2subdomain1 , layer2subdomain2) )
        # subdomainlists =   ( (layer1subdomain1 , layer2subdomain1) , (layer1subdomain2 , layer2subdomain2) )

        # The listofsubdomains is inverted using zip and map to create a list of lists
        # so that each subdomain is grouped with the corresponding subdomain from each layer (see above for example)
        subdomainlists = map(list,zip(*listofsubdomains))

        return subdomainlists

    def executor(self,subdomains):
        """ Executor handles processing of the function by iterating over locations in a subdomain
        :return: #TODO: Undefined return value.
        """
        PCMLTODO("executor assumes single subdomain as output, which is not universal for all operations")

        outsubdomain = subdomains.pop(0)
        outarr = outsubdomain.get_nparray()

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
        PCMLOperationError("Operation function is not implemented")

