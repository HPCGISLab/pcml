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
    for i in locations:
        a=i['v']
        for j in locations:
            b=j['v']
            if a>b:
                c = a
    return c
    
@localoperation
def LocalMinimum(self,locations,subdomains):
    for i in locations:
        a=i['v']
        for j in locations:
            b=j['v']
            if a<b:
                c = a
    return c

@localoperation
def LocalMean(self,locations,subdomains):
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
def LocalTangent(self, locations, subdomains):
    # NOTE: Assumes 1 subdomain
    return math.tan(locations[0]['v'])

@localoperation
def LocalVariety(self,locations,subdomains):
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
def LocalMajority(self,locations,subdomains):
    #returns most frequent of each location
    temp=[]
    count=0
    #print locations
    for loc in locations:
        count+=1
        if count<=len(locations):
            temp.append(loc['v'])
        else:
            count=0
    return most_common(temp)
def most_common(lst):
    c_len = 0
    max_len = 0
    c_i = 0
    max_i = 0
    c_item = None
    max_item = None
    for i, item in sorted(enumerate(lst), key=lambda x: x[1]):
        if c_item is None or c_item != item:
            if c_len > max_len or (c_len == max_len and c_i < max_i):
                max_len = c_len
                max_i = c_i
                max_item = c_item
            c_len = 1
            c_i = i
            c_item = item
        else:
            c_len += 1
    if c_len > max_len or (c_len == max_len and c_i < max_i):
        return c_item
    return max_item
    

