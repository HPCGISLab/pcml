#!/usr/bin/python
"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu)
"""
from pcml import *

import os.path as path
import sys,getopt

import numpy as np

if __name__ == '__main__':

    # See if number of cores was provided as a command-line parameter
    # This script can be re-executed with 1, 2, 4, 8, 16 cores 
    try:
        opts,args=getopt.getopt(sys.argv[1:],"n:",["ncores"])
    except getopt.GetoptError:
        print("bed-and-breakfast-procedure.py -n <cores>")
        sys.exit(1)

    for opt, arg in opts:
        if opt in ("-n", "--num"):
            PCMLConfig.num_procs=int(arg)
            # Swith to serial evalution if n is 1
            if PCMLConfig.num_procs == 1:
                PCMLConfig.exectype = ExecutorType.serialpython

    # The decomposition granularity can be re-configured here
    PCMLConfig.decomposition_granularity = 500
    #PCMLConfig.decomposition_granularity = 800
    
    print("Number of cores used for spatial data processing: ",PCMLConfig.num_procs)

    thebuffersize=600

    # Default data directory
    datadir="./ijgis-data/"

    # Read 2 layers
    landcover = ReadGeoTIFF(path.join(datadir,"ohio_nlcd.tif"))
    transportation = ReadGeoTIFF(path.join(datadir,"ohio_transportation.tif"))

    # Classification type 11 is open water in NLCD
    #near_water=BufferedClassify(landcover,buffersize=thebuffersize,classtype=11)

    # Buffer+Classify
    water_buf=FocalBuffer(landcover,buffersize=thebuffersize,classtype=11)
    near_water=LocalClassify(water_buf,classtype=11)

    in_forest=LocalClassify(landcover,classtype=41)

    # Classification type 1400 in Transportation is non-major roadway (local or rural road)
    # Buffered Classify
    #near_road=BufferedClassify(transportation,buffersize=thebuffersize,classtype=1400)

    # Buffer+Classify
    road_buf=FocalBuffer(transportation,buffersize=thebuffersize,classtype=1400)
    near_road=LocalClassify(road_buf,classtype=1400)
    

    suitable = LocalAnd(near_water, in_forest, near_road)

    WriteGeoTIFF(path.join(datadir,"suitable_map.tif"),suitable)

    exit()

