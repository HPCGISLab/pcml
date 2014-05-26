"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from pCML import *
from pCML.util.Messaging import PCMLNotImplemented
from pCML.util.LayerBuilder import lst_to_layer
from numpy.ma import allequal
import numpy as np
import unittest

class TestLayer(unittest.TestCase):
    def test_layer_repr(self):
        l1 = Layer(0,0,0,0,'test_layer1')
        r1 = repr(l1)
        r1_res = "<Layer: (%f,%f) [%f,%f] : %s>" % (0,0,0,0,'test_layer1')
        self.assertEqual(r1, r1_res)

        not_supported_decomposition_method = -1
        self.assertRaises(PCMLNotImplemented, l1.decomposition, not_supported_decomposition_method, 2)
