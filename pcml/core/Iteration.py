"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); 
"""
from .Layer import *
from .Subdomain import *
import pcml.core.PCMLConfig as PCMLConfig

def rowmajoriteration(subdomain):
        if subdomain.data_structure!=Datastructure.array:
           PCMLNotSupported("subdomain.__iter__ currently assumes an array data structure")

        # Passes a locind with r,c set to *absolute* cell location in the layer
        # This makes it possible to reference buffered boxes (halo/ghost zones)
        for locind in itertools.product(xrange(subdomain.nrows), xrange(subdomain.ncols)):
            yield {'r':locind[0]+subdomain.r,'c':locind[1]+subdomain.c}

    
def columnmajoriteration(subdomain):
        if subdomain.data_structure!=Datastructure.array:
           PCMLNotSupported("subdomain.__iter__ currently assumes an array data structure")

        # Passes a locind with r,c set to *absolute* cell location in the layer
        # This makes it possible to reference buffered boxes (halo/ghost zones)
        for locind in itertools.product(xrange(subdomain.ncols), xrange(subdomain.nrows)):
            yield {'c':locind[0]+subdomain.c,'r':locind[1]+subdomain.r}

