"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
import os
from .PCMLPrims import *
from ..util.Messaging import *


if os.name == 'nt':
    print(PCMLUserInformation("Windows cannot support parallel computation in PCML, enabling serial implementation"))
    # Number of processes to run in Windows is 1
    num_procs = 1

    # Use the serial implementation in Windows
    exectype = ExecutorType.serialpython

    # Defines the default granularity for spatial domain decomposition
    decomposition_granularity = 1
    # Defines number of subdomains for point decomposition
    numsubdomains = 16
else:
    # Default number of processes to run in Linux
    num_procs = 4

    # Use the parallel queue implementation by default
    exectype = ExecutorType.parallelpythonqueue

    # Defines the default granularity for spatial domain decomposition
    decomposition_granularity = 16
    # Defines number of subdomains for point decomposition
    numsubdomains = 16


# Number of processes to run
# num_procs=4

# exectype=ExecutorType.serialpython
# exectype=ExecutorType.parallelpythonqueue

# The precision used in formatting floating values into strings
value_precision = "%f"

# By default osgeo including gdal, ogr, and osr are not available
# In OperationIO we try to import them and if successful osgeoenabled=1
osgeoenabled = 0
