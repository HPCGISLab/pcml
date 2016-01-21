"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from ..util.Messaging import *
from .Decomposition import *
from .BoundingBox import *
from .Iteration import *
from .PCMLPrims import *
from .PCMLConfig import *
from abc import ABCMeta, abstractmethod
import multiprocessing as mp
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
        self.iteration = kwargs.get('iteration',rowmajoriteration) # By default use product-based iteration 
        self.duration = kwargs.get('duration',None) # By default use no duration 
        self.curtime = kwargs.get('curtime',None) # By default use no curtime (t) 
        #adding this to get the operation specified parameter
        self.kwargs=kwargs
        #Boolean to check whether the operation gives back a layer as output or not. Defaults to true
        self.isnonlayeroutput=kwargs.get('nonlayeroutput',False)
        #Attribute to store single output
        self.nonlayeroutput=None
        if self.opclass==OpClass.localclass and self.buffersize != 0:
            raise PCMLOperationError("Buffersize should be 0 for localclass currently %s" % self.buffersize)
        #If zonal operation we want the entire layer data
        if self.opclass==OpClass.zonalclass:
            self.buffersize=999999999999
        #attribute to store partial results created by sudomains which can be used between processes
        # FIXME: This will be used for operations that do not return layers, but instead global results.
        # This feature is currently not supported
        #if PCMLConfig.exectype==ExecutorType.serialpython:
        #    self.globalresults=[]
        #elif PCMLConfig.exectype==ExecutorType.parallelpythonqueue:
        #    self.globalresults=mp.Manager().list()

    def __repr__(self):
        return "<Operation: %s : %i layers>" % (self.name,len(self._layers))

    def getOutputLayers(self):
        PCMLTODO("Need to support more than one output layer")
        #check whether output is a nonlayer output or not
        if self.isnonlayeroutput:
            return self.nonlayeroutput
        return self._layers[0]

    def _decompositioninit(self):
        # Duplicate a layer to create an output layer with the correct dimensions
        # Get the first layer
        firstlayer=self._layers[0]

        # The output layer is a duplicate of the first layer
        outputlayer=firstlayer.duplicate()
        outputlayer.title="Output for operation %s"%self.name

        self._layers.insert(0, outputlayer)  # Add the output layer to the front of the layers list

    def _decompositionrun(self):
        """ Divides the :member:_layers into subdomains for further processing.
        The decomposition method is defined by :member:`decompositionmethod`.
        You can also define you own decomposition algorithm by overriding this method.
        """

        PCMLTODO("Need to support multiple output layers, this can be done by overriding decomposition and inserting multiple output layers")

        # Run decomposition to generate a list of (buffered) bounding boxes to decompose each layer into
        bufboundingboxlist = self.decomposition(self._layers,self.buffersize)

        # Initialize decomposition by creating the output layer and inserting it into the layer's list of subdomains
        self._decompositioninit()

        # Pass in a list of layers
        unbuf_boundingboxlist = self.decomposition([self._layers[0]],0)

        # Sanity check: We expect at least one bounding box from decomposition
        assert(len(bufboundingboxlist)>0)

        #print "Buffered bounding boxes",bufboundingboxlist
        #print "Unbuffered bb boxes    ",unbuf_boundingboxlist
       
        listofsubdomains = [] 
        for layer in self._layers:
            # Transform the layer to a list of subdomains
            # Add the list of 'listofsubdomains' to create a list of lists
            if layer == self._layers[0]: # If it is the output subdomain, disable buffer
                subdomains = layer.tosubdomains(unbuf_boundingboxlist)
            else: # Otherwise don't
                subdomains = layer.tosubdomains(bufboundingboxlist)

            # Set the layer for each subdomain
            for sd in subdomains:
                sd.layer = layer

            listofsubdomains.append(subdomains)


        # FIXME: This is deprecated code that needs to be removed.
        ''' REMOVE THEM 
        # The output layer is the first layer in the layers list (self.layers[0])
        # Decompose it with a 0 buffer
        #listofsubdomains.append(self._layers[0].decomposition(self.decomposition_method, 0))
        listofsubdomains.append(self.decomposition(self._layers[0], 0))
        if self._layers[0].data_structure==Datastructure.pointlist:
            self._layers[0].set_pointlist([])
        for layer in self._layers:
            if layer != self._layers[0]: # Skip the output layer, because it was already decomposed and added
                #listofsubdomains.append(layer.decomposition(self.decomposition_method, self.buffersize)) # buffer size is set based on classification (L,F,Z,G)
                if self.decomposition.__name__=='pointrasterrowdecomposition':
                    listofsubdomains.append(self.decomposition(layer,self.buffersize,layerlist=self._layers))
                else:
                    # Create a subdomain and populate it with the correct attribute values
                    listofsubdomains.append(self.decomposition(layer, self.buffersize)) # buffer size is set based on classification (L,F,Z,G)
        '''

        # The listofsubdomains is inverted using zip and map to create a list of lists
        # so that each subdomain is grouped with the corresponding subdomain from each layer (see example below)
        subdomainlists = map(list,zip(*listofsubdomains))

        # listofsubdomains = ( (layer1subdomain1 , layer1subdomain2) , (layer2subdomain1 , layer2subdomain2) )
        # subdomainlists =   ( (layer1subdomain1 , layer2subdomain1) , (layer1subdomain2 , layer2subdomain2) )

        return subdomainlists

    # By default we use rowdecomposition as our decomposition method
    # Users may override decomposition with any other method they would like using kwargs (see __init__)
    #def decomposition(layer,buffersize):

    # By default we use rowmajoriteration as our iteration method
    # Users may override iteration with any other method they would like using kwargs (see __init__)
    #def iteration(subdomain):

    def executor(self,subdomains):
        """ Executor handles processing of the function by iterating over locations in a subdomain
        :return: #TODO: Undefined return value.
        """
        PCMLTODO("executor assumes single subdomain as output, which is not universal for all operations")
        outsubdomain = subdomains.pop(0)

        #print "subdomains",subdomains
        #for sd in subdomains:
        #    print "sdlen(",sd,")=",len(sd._data)
        print "start r,c",outsubdomain.r,outsubdomain.c

        if outsubdomain.data_structure==Datastructure.pointlist:
            pointlist=outsubdomain.get_pointlist()
            for i in xrange(len(pointlist)):
                val=self.function([pointlist[i]],subdomains)
                newdict=pointlist[i].copy()
                newdict['v']=val
                pointlist[i]=newdict
            if PCMLConfig.exectype==ExecutorType.serialpython:
                self._layers[0].get_pointlist().extend(outsubdomain.get_pointlist())

            print "pointlist datastructure"
        elif outsubdomain.data_structure==Datastructure.array:
            if outsubdomain.r == 0:
                print "in elif : r",outsubdomain.r
            outarr = outsubdomain.get_nparray()
            # Iterate over locations in the outsubdomain using iteration method and apply function to each location
            for loc in self.iteration(outsubdomain):
                l = [] # Create an empty list to store locations
                for sd in subdomains:
                    if sd.data_structure!=Datastructure.array: # Skip non array subdomains
                        continue
                    # Get a location in this subdomain with same coordinates as locind
                    locv=sd.get_locval(loc)
                    l.append(locv) # append to list of locations
                val = self.function(l,subdomains) # Apply function to all locations
                outarr[loc['r']-outsubdomain.r][loc['c']-outsubdomain.c]=val # Set val to outarr at locind

            print "r,c",outsubdomain.r,outsubdomain.c
        print "Executor finished",outsubdomain

    def writepointdatatooutputlayer(self,subdomainlists):
        for subdomains in subdomainlists:
            self._layers[0].get_pointlist().extend(subdomains[0].get_pointlist())

    def function(self,locations,subdomains):
        raise PCMLOperationError("Operation function is not implemented")

