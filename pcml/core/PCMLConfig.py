"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from Scheduler import *

# Number of processes to run
num_procs=4

exectype=ExecutorType.serialpython
exectype=ExecutorType.parallelpythonqueue

# The precision used in formatting floating values into strings
value_precision="%f"
