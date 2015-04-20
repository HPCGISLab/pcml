"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu);
Sandeep Vutla (svutla@kent.edu, sandeepvutla@yahoo.in); Gowtham Kukkadapu (gkukkada@kent.edu)
"""
from ..core.Operation import *
from ..core.Scheduler import *
from ..util.OperationBuilder import *
import numpy as np
import types
import math

@focaloperation
def FocalMean(self, locations, subdomains):
    # Begin with a zero sum and count to calculate the mean
    # This function assumes one location and one subdomain are passed in
    sumct=count=0.0
    arr=subdomains[0].bufferedlocgetarr(locations[0],self.buffersize)
    for i in arr.flat:
        sumct+=i
        count+=1.0
    return sumct/count

@focaloperation
def FocalMean_np(self, locations, subdomains):
    # Begin with a zero sum and count to calculate the mean
    # This function assumes one location and one subdomain are passed in
    arr=subdomains[0].bufferedlocgetarr(locations[0],self.buffersize)

    return np.average(arr)

@focaloperation
def FocalContourLines(self, locations, subdomains): # Experimental
    # This function assumes one location and one subdomain are passed in
    cntrline=5.0
    sumct=count=0.0
    arr=subdomains[0].bufferedlocgetarr(locations[0],self.buffersize)
    if(arr.size<9):
        return subdomains[0].nodata_value
    (q,r)=divmod(arr[1][1],cntrline)

    cutoff=(q+1)*cntrline # If any neighboring value is above this, then this point is part of a contour line
    if(arr[2][1]>=cutoff or arr[1][2]>=cutoff or arr[0][1]>=cutoff or arr[1][0]>=cutoff):
        return 1
    return 0
    for i in arr.flat:
        sumct+=i
        count+=1.0
    avg=sumct/count
    if(int(avg)%10):
        return 1
    return 0

@focaloperation
def FocalAspect(self, locations, subdomains): 
    # Begin with a zero sum and count to calculate the mean
    # This function assumes one location and one input subdomain are passed in
    sumct=count=0
    arr=subdomains[0].bufferedlocgetarr(locations[0],1) # Set buffersize to 1
    if(arr.size<9):
        return subdomains[0].nodata_value
    dzdx=((arr[2][0]+2*arr[1][0]+arr[0][0]) - (arr[2][2]+2*arr[1][2]+arr[0][2])) / 8
    dzdy=((arr[0][0]+2*arr[0][1]+arr[0][2]) - (arr[2][0]+2*arr[2][1]+arr[2][2])) / 8
    # Flipped
    #dzdx=((arr[2][2]+2*arr[1][2]+arr[0][2]) - (arr[2][0]+2*arr[1][0]+arr[0][0])) / 8
    #dzdy=((arr[2][0]+2*arr[2][1]+arr[2][2]) - (arr[0][0]+2*arr[0][1]+arr[0][2])) / 8
    # 57.2... = 360/(2*pi)
    #aspect=270 - 57.29577951308232 * math.atan2(dzdy,-dzdx)
    #if(aspect>360):
        #aspect=aspect-360.0
    #return aspect
    aspect=57.29578 * math.atan2(dzdy,-dzdx)
    if(aspect<0):
        return 90.0-aspect
    elif(aspect>90.0):
        return 360-aspect+90
    else:
        return 90.0-aspect


@focaloperation
def HillShade(self, locations, subdomains): # Experimental
    altitude=45
    #altitude=60 # Override standard
    azimuth=315
    rtod=0.017453292519943295 # ( pi / 180.0 )

    arr=subdomains[0].bufferedlocgetarr(locations[0],1) # Set buffersize to 1
    if(arr.size<9):
        return subdomains[0].nodata_value
    zenith_rad=(90-altitude)*rtod
    azimuth_rad=(360.0-azimuth+90)*rtod

    dzdx=((arr[2][2]+2*arr[1][2]+arr[0][2]) - (arr[2][0]+2*arr[1][0]+arr[0][0])) / (8*subdomains[0].cellsize)
    dzdy=((arr[2][0]+2*arr[2][1]+arr[2][2]) - (arr[0][0]+2*arr[0][1]+arr[0][2])) / (8*subdomains[0].cellsize)

    aspect=math.atan2(dzdy,-dzdx)
    slope=math.atan( 1 * math.sqrt(dzdx*dzdx + dzdy*dzdy) )

    shade = 255.0 * ( (math.cos(zenith_rad) * math.cos(slope) ) +
                    (  math.sin(zenith_rad) * math.sin(slope) * math.cos(azimuth_rad-aspect) ) )
    return shade

@focaloperation
def FocalMajority(self, locations, subdomains):
    #Note: Takes two or more layers as input and returns majority element respective to the location
    chunk=[]
    for i in range(0,len(locations)):
        arr=subdomains[i].bufferedlocgetarr(locations[i],self.buffersize)
        chunk+=arr.tolist()
    return majority( np.asarray( chunk).flatten('C'))
def majority(arr):
    c_len = 0
    max_len = 0
    c_i = 0
    max_i = 0
    c_item = None
    maj_item = None
    for i, item in sorted(enumerate(arr), key=lambda x: x[1]):
        if c_item is None or c_item != item:
            if c_len > max_len or (c_len == max_len and c_i < max_i):
                max_len = c_len
                max_i = c_i
                maj_item = c_item
            c_len = 1
            c_i = i
            c_item = item
        else:
            c_len += 1
    if c_len > max_len or (c_len == max_len and c_i < max_i):
        return c_item
    return maj_item

@focaloperation
def FocalMinority(self, locations, subdomains):
    #takes two or more layers as input and returns minority element in respective locaitons
    chunk=[]
    for i in range(0,len(locations)):
        arr=subdomains[i].bufferedlocgetarr(locations[i],self.buffersize)
        chunk+=arr.tolist()
    print np.asarray( chunk).flatten('C')
    return minority(np.asarray( chunk).flatten('C'))
                    
def minority(temp):
    ctr = Counter(temp)
    keys = ctr.keys()
    vals = ctr.values()
    least = []
    m = min(vals)
    for i in range(0,len(vals)):
        if vals[i] == m:
            least.append(keys[i])
    if len(least)==1:
        return least
    else:
        return least[:1] 
        
@focaloperation
def FocalMaximum(self,locations,subdomains):
    #can take multiple layers and returns maximum number 
    chunk=[]
    for i in range(0,len(locations)):
        arr=subdomains[i].bufferedlocgetarr(locations[i],self.buffersize)
        chunk+=arr.tolist()
    return maximum( np.asarray( chunk).flatten('C'))
def maximum(arr):
    iterator = iter(arr)
    m = next(iterator)
    for item in iterator:
        if item > m:
            m = item
    return m
    
@focaloperation
def FocalMinimum(self,locations,subdomains):
    #can take multiple layers and returns minimum number 
    chunk=[]
    for i in range(0,len(locations)):
        arr=subdomains[i].bufferedlocgetarr(locations[i],self.buffersize)
        chunk+=arr.tolist()
    return minimum( np.asarray( chunk).flatten('C'))
def minimum(arr):
    iterator = iter(arr)
    m = next(iterator)
    for item in iterator:
        if item < m:
            m = item
    return m
    
@focaloperation
def FocalMaximum_np(self, locations, subdomains):
    #returns maximum number in focal proximity of multiple layers
    chunk=[]
    for i in range(0,len(locations)):
        arr=subdomains[i].bufferedlocgetarr(locations[i],self.buffersize)
        chunk+=arr.tolist()
    return np.max( np.asarray( chunk).flatten('C'))

@focaloperation 
def FocalMinimum_np(self, locations, subdomains):
    #returns minimum number in focal proximity of multiple layers
    chunk=[]
    for i in range(0,len(locations)):
        arr=subdomains[i].bufferedlocgetarr(locations[i],self.buffersize)
        chunk+=arr.tolist()
    return np.min( np.asarray( chunk).flatten('C'))

@focaloperation
def FocalSum(self, locations, subdomains):
    arr=subdomains[0].bufferedlocgetarr(locations[0],self.buffersize) 
    sum=0.0
    for i in arr.flatten('F'):
        sum+=i
    sum=sum
    return sum

@focaloperation
def FocalPercentage(self, locations, subdomains):
    #Begin with zero numerator and denominator
    #Increment deominator by one for each loop and numerator by one if the neighborhood value matches the location value
    num=denom=0.0
    arr=subdomains[0].bufferedlocgetarr(locations[0],self.buffersize)
    for i in arr.flat:
        if locations[0]['v']==i:
                num+=1.0
        denom+=1.0
    return (num/denom)*100
