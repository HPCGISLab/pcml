"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from pCML import *
from pCML.util.LayerBuilder import lst_to_layer
from numpy.ma import allequal
import numpy as np
import unittest

class TestBoundingBox(unittest.TestCase):
    def test_boundingbox_repr(self):
        b1 = BoundingBox(0,0,0,0)
        r1 = repr(b1)
        r1_res = "<BoundingBox: (%f,%f) [%f,%f]>" % (0,0,0,0)
        self.assertEqual(r1,
                r1_res,
                "BoundingBox geometory does not match")

        b2 = BoundingBox(3, 4, 7, 1)
        r2 = repr(b2)
        self.assertEqual(r2,
                "<BoundingBox: (%f,%f) [%f,%f]>" % (3,4,7,1),
                "Boundingbox geometory does not match")

        self.assertEqual(b2._data, None)
        self.assertEqual(b2.print_data(), None)

        b2.set_pointlist([{'x': 1, 'y':1, 'z':1}])
        self.assertRaises(PCMLNotImplemented, b2.print_data)

    def test_set_nparray(self):
        b1 = BoundingBox(0,0,100,2)
        arr = np.asarray([[5] * 4] * 10)
        b1.set_nparray(arr, 1, 0)
        self.assertTrue(allequal(b1._data, arr))

