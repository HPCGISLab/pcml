#!/usr/bin/python
"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from pcml import *

import os.path as path
import sys,getopt
import timeit
import numpy as np
import itertools
import timeit
if __name__ == '__main__':
    # Default data directory
    #modified the queue sleep time in scheduler as it was affecting speed.Link for reference http://programming.nullanswer.com/forum/3177367
    datadir='/home/jajayaku/pcmlprodfortest/pcml/data/'
    corearray=[1,2,4,8,16]
    decompositionstrategy=[columndecomposition]
    granularityvalues=[2532]
    alldata=[corearray,decompositionstrategy,granularityvalues]   
    layerA=ReadGeoTIFF(path.join(datadir,"datahuge.img"))
    fileout=open(path.join(datadir,"process5GBCol1subpercore.csv"),'w')
    configstring='Cores,Decomposition,Granularity,Operation,Time\n'
    fileout.write(configstring)
    #layerA=ReadASCIIGrid(path.join(datadir,"datab.asc"))
    layerB=layerA.duplicate()  
    layerO=layerA.duplicate()
    layerBarray=layerB.get_nparray()
    layerBarray+=1
    layerC=LocalClassify20_np(layerA)
    allcomb=list(itertools.product(*alldata)) 
    for config in allcomb:
	valdict={}
        core=config[0]
        decomp=config[1]
        granularity=config[2]
        PCMLConfig.num_procs=core
        PCMLConfig.decomposition_granularity=granularity
        #test local sum
        start=timeit.default_timer()
        layer_o=ZonalMean_exec(layerA,layerC,decomposition=decomp,outputlayer=layerO)
        end=timeit.default_timer()
        valdict['ZonalMeanExec']=end-start
        start=timeit.default_timer()
        layer_o=LocalSum_np(layerA,layerC,decomposition=decomp,outputlayer=layerO)
        end=timeit.default_timer()
        valdict['LocalSumExec']=end-start
        start=timeit.default_timer()
        layer_o=FocalMean_np_exec_filter(layerA,decomposition=decomp,buffersize=2,outputlayer=layerO)
        end=timeit.default_timer()
        valdict['FocalMeanExec']=end-start
        for key in valdict:
            resultstring=str(core)+','+decomp.__name__+','+str(granularity)+','+key+','+str(valdict[key])+'\n'
            fileout.write(resultstring)
            fileout.flush()
    fileout.close()

