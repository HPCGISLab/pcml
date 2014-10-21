"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from ..core.Operation import *
from ..core.Scheduler import *
from ..util.OperationBuilder import *
import numpy as np
import types
import math

@globaloperation
def GlobalMinMHDistance(self, locations, subdomains):
    pointlist=subdomains[1].get_pointlist()

    locind=locations[0]
    loc=subdomains[0].get_yxloc(locind) # Convert from array coordinates (r,c) to (y,x) coordinates

    mindst=999999999.0
    mindstindex=-1
    for index in xrange(len(pointlist)):
        point=pointlist[index]
        dst=abs(loc['y'] - point['y']) + abs(loc['x'] - point['x'])
        if dst<mindst:
            mindst=dst
            mindstindex=index
    return mindst

@globaloperation
def GlobalMinMHDistanceIndex(self, locations, subdomains):
    pointlist=subdomains[1].get_pointlist()

    locind=locations[0]
    loc=subdomains[0].get_yxloc(locind) # Convert from array coordinates (r,c) to (y,x) coordinates

    mindst=999999999.0
    mindstindex=-1
    for index in xrange(len(pointlist)):
        point=pointlist[index]
        dst=abs(loc['y'] - point['y']) + abs(loc['x'] - point['x'])
        if dst<mindst:
            mindst=dst
            mindstindex=index
    return mindstindex

@globaloperation
def GlobalMinDistanceIndex(self, locations, subdomains):
    pointlist=subdomains[1].get_pointlist()

    locind=locations[0]
    loc=subdomains[0].get_yxloc(locind) # Convert from array coordinates (r,c) to (y,x) coordinates
    #print "loc",loc

    mindst=999999999.0
    mindstindex=-1
    for index in xrange(len(pointlist)):
        point=pointlist[index]
        #print "point",point
        dst=math.sqrt((loc['y'] - point['y']) ** 2 + (loc['x'] - point['x']) ** 2)
        if dst<mindst:
            mindst=dst
            mindstindex=index
    #print "mindst",mindst
    return mindstindex

@globaloperation
def GlobalMinDistance(self, locations, subdomains):
    pointlist=subdomains[1].get_pointlist()

    locind=locations[0]
    loc=subdomains[0].get_yxloc(locind) # Convert from array coordinates (r,c) to (y,x) coordinates

    mindst=999999999.0
    for point in pointlist:
        #print "point",point
        dst=math.sqrt((loc['y'] - point['y']) ** 2 + (loc['x'] - point['x']) ** 2)
        if dst<mindst:
            mindst=dst
    #print "mindst",mindst
    return mindst

