"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from pcml import *
import unittest

class PCMLSerialTestCase(unittest.TestCase):
    def setUp(self):
        PCMLConfig.num_procs = 1
        PCMLConfig.exectype = ExecutorType.serialpython

class PCMLParallelTestCase(unittest.TestCase):
    def setUp(self):
        PCMLConfig.num_procs = 4
        PCMLConfig.exectype = ExecutorType.parallelpythonqueue
