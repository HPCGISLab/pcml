#!/usr/bin/python
"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from pcml import *

import os.path as path
import sys,getopt

import numpy as np

if __name__ == '__main__':

    # See if number of cores was provided as a command-line parameter
    try:
        opts,args=getopt.getopt(sys.argv[1:],"n:",["ncores"])
    except getopt.GetoptError:
        print "test.py -n <cores>"
        sys.exit(1)

    for opt, arg in opts:
        if opt in ("-n", "--num"):
            PCMLConfig.num_procs=int(arg)
            # Swith to serial evalution if n is 1
            if PCMLConfig.num_procs == 1:
                PCMLConfig.exectype = ExecutorType.serialpython
    
    print "Number of cores used for spatial data processing: ",PCMLConfig.num_procs

    bb=BoundingBox(0,0,10,10)
 
    # Default data directory
    datadir="./data/"

    # Read 2 test layers
    layer1=ReadASCIIGrid(path.join(datadir,"datab.asc"))
    layer2=ReadASCIIGrid(path.join(datadir,"dataa.asc"))
   
    # Print out the layers 
    print "layer1",layer1
    layer1.print_data()
    print "layer2",layer2
    layer2.print_data()

    print "\nTesting layer division" 
    layero=layer1/layer2
    layero.print_data()

    print "\nTesting layer division" 
    layero=layer1/layer1
    layero.print_data()

    print "\nTesting layer addition and multiplication"
    layero=(layer1+layer1)*layer2
    print "layero",layero
    layero.print_data()

    print "\nTesting layer subtraction layer1-layer1 "
    layero=layer1-layer1
    print "layero",layero
    layero.print_data()


    layero=LocalSum_np(layer1,layer2)
    print "layero", layero 
    layero.print_data()

    print "\nTesting FocalMean"
    layer1=FocalMean(layer1) 
    layer1=LocalSin(layer1) 
    layer3=FocalMean(layer1, buffersize=3)

    print "layer1 (FocalMean with buffer = 1)", layer1
    layer1.print_data()
    print "layer3 (FocalMean with buffer = 3)", layer3
    layer3.print_data()

    print "\nWriting layer1 and layer3 out as an ASCII grid file and GeoTIFF"
    print "Writing",datadir,"layer1.asc"
    WriteASCIIGrid(path.join(datadir,"layer1.asc"), layer1)

    # Print please note that GeoTIFF writing will only work with the GDAL package installed.
    print "Writing",datadir,"layer3.tif"
    WriteGeoTIFF(path.join(datadir,"layer3.tif"), layer3)

