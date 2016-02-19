"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu);
Sandeep Vutla(svutla@kent.edu, sandeepvutla@yahoo.in); Gowtham Kukkadapu(gkukkada@kent.edu)
"""
from ..core.Operation import *
from ..core.Scheduler import *
from ..util.OperationBuilder import *
import numpy as np
import types
import math

@localoperation
def LocalSum(self, locations, subdomains):
    # LocalSum will return the sum of values at each location for two or more layers.
    val = 0
    for loc in locations:
        val += loc['v']
    return val

@localoperation
def LocalMult(self, locations, subdomains):
    val = 1
    for loc in locations:
        val *= loc['v']
    return val

@localoperation
def LocalDivision(self, locations, subdomains):
    # NOTE: Assumes 2 subdomains for division.
    location0=locations[0]
    location1=locations[1]
    return location0['v']/location1['v']


@localoperation
def LocalSubtraction(self, locations, subdomains):
    # NOTE: Assumes 2 subdomains for subtraction.
    location0=locations[0]
    location1=locations[1]
    return location0['v']-location1['v']


@localoperation
def LocalSin(self, locations, subdomains):
    # NOTE: Assumes 1 subdomain
    return math.sin(locations[0]['v'])

@localoperation
def LocalMaximum(self,locations,subdomains):
    #Takes number of layers and returns the maximum value among the layers in particular location.
    mx=locations[0]['v']
    for i in locations:
        if mx< i['v']:
            mx= i['v']
    return mx
    
@localoperation
def LocalMinimum(self,locations,subdomains):
    #Takes number of layers and returns the minimum value among the layers in particular location.
    mi=locations[0]['v']
    for i in locations:
        if mi> i['v']:
            mi= i['v']
    return mi

@localoperation
def LocalMean(self,locations,subdomains):
    # Takes two or more layers and returns mean
    val = 0
    count=0
    for loc in locations:
        count+=1
        val += loc['v']
    return val/count

@localoperation
def LocalDifference(self,locations,subdomains):
    # LocalDifference will return the Difference of values at each location for two or more layers.
    val = locations[1]['v']
    for loc in locations:
        val -= loc['v']
        #print loc['v']
    return val

@localoperation
def LocalVariety(self,locations,subdomains):
    # Takes two or more layers and returns number of dissimilar values associated with respective location
    temp=[]
    count=0
    for loc in locations:
        count+=1
        if count<=len(locations):
            temp.append(loc['v'])
        else:
            count=0
    return len(variety(temp))    
def variety(i):
    return [] if i==[] else [i[0]] + variety(filter(lambda x: x!= i[0], i[1:]))

@localoperation
def LocalSine(self, locations, subdomains):
    # NOTE: Assumes 1 subdomain
     #takes single layer as input and returns the computed Sine of each location's value
    return math.sin(locations[0]['v'])
    
@localoperation
def LocalCosine(self,locations,subdomains):
    # NOTE: Assumes 1 subdomain 
    #takes single layer as input and returns the computed Cosine of each location's value
    return math.cos(locations[0]['v'])
    
@localoperation
def LocalTangent(self, locations, subdomains):
    # NOTE: Assumes 1 subdomain
     #takes single layer as input and returns the computed Tangent of each location's value
    return math.tan(locations[0]['v'])

@localoperation
def LocalArcCosine(self,locations,subdomains):
    # NOTE: Assumes 1 subdomain, takes layer as input and returns the computed ArcCosine of each location's value
     #takes single layer as input and returns the computed ArcCosine of each location's value
    return math.acos(locations[0]['v'])

@localoperation
def LocalArcSine(self,locations,subdomains):
    # NOTE: Assumes 1 subdomain, takes layer as input and returns the computed ArcSine of each location's value
    #takes single layer as input and returns the computed ArcSine of each location's value
    return math.asin(locations[0]['v'])

@localoperation
def localArcTangent(self,locations,subdomains):
    # NOTE: Assumes 1 subdomain, takes layer as input and returns the computed ArcTangent of each location's value
    #takes single layer as input and returns the computed ArcTangent of each location's value
    return math.atan(locations[0]['v'])

@localoperation
def LocalNotEqual(self, locations, subdomains):
    # NOTE: Assumes 2 subdomains for not equal.
    # Cheng Zhang 02/19/2016
    return location[0]['v']!=location[1]['v']

@localoperation
def LocalLessThan(self, locations, subdomains):
    # NOTE: Assumes 2 subdomains for less than.
    # Cheng Zhang 02/19/2016
    return location[0]['v']<location[1]['v']

@localoperation
def LocalGreaterThan(self, locations, subdomains):
    # NOTE: Assumes 2 subdomains for greater than.
    # Cheng Zhang 02/19/2016
    return location[0]['v']>location[1]['v']
