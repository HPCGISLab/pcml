"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from pCML import *
from pCML.util.LayerBuilder import *
from numpy.ma import allequal
import numpy as np
import unittest

class TestOperationBuilder(unittest.TestCase):
    def setUp(self):
        pass

    def test_operation(self):
        @operation(opclass=OpClass.localclass, override='function')
        def LocalOperation(self, locations, subdomains):
            pass
        self.assertTrue(LocalOperation._PCML_exported)

    def test_localoperation(self):
        @localoperation
        def LocalOperation(self, locations, subdomains):
            pass
        self.assertTrue(LocalOperation._PCML_exported)

        @localoperation
        @operation()
        def LocalOperation2(self, locations, subdomains):
            pass
        self.assertTrue(LocalOperation2._PCML_exported)
        self.assertEqual(LocalOperation2.opclass, OpClass.localclass)

        test_arr = [[1,2,3], [4,5,6], [7,8,9]]
        test_layer = lst_to_layer(test_arr)
        def test_equal(a, b):
            self.assertEqual(a, b)

        @localoperation
        @operation(opclass=OpClass.localclass, override='function', test_mark=True)
        def LocalAddOne(self, locations, subdoamins):
            test_equal(self.test_mark, True)
            for i in xrange(len(locations)):
                loc = locations[i]
                test_equal(loc['v'], test_arr[loc['r']][loc['c']])
            return locations[0]['v']+1

        test_res = [[2,3,4], [5,6,7], [8,9,10]]
        res = LocalAddOne(test_layer)
        self.assertTrue(allequal(test_res, res._data))


    def test_focaloperation(self):
        @focaloperation
        def FocalOperation(self, locations, subdomains):
            pass
        self.assertTrue(FocalOperation._PCML_exported)

        @focaloperation
        @operation()
        def FocalOperation2(self, locations, subdomains):
            pass

        self.assertTrue(FocalOperation2._PCML_exported)
        self.assertEqual(FocalOperation2.opclass, OpClass.focalclass)

    def test_globaloperation(self):
        @globaloperation
        def GlobalOperation(self, locations, subdomains):
            pass
        self.assertTrue(GlobalOperation._PCML_exported)

        @globaloperation
        @operation()
        def GlobalOperation2(self, locations, subdomains):
            pass

        self.assertTrue(GlobalOperation2._PCML_exported)
        self.assertEqual(GlobalOperation2.opclass, OpClass.globalclass)

    def test_zonaloperation(self):
        @zonaloperation
        def ZonalOperation(self, locations, subdomains):
            pass
        self.assertTrue(ZonalOperation._PCML_exported)

        @zonaloperation
        @operation()
        def ZonalOperation2(self, locations, subdomains):
            pass
        self.assertTrue(ZonalOperation2._PCML_exported)
        self.assertEqual(ZonalOperation2.opclass, OpClass.zonalclass)

    def test_decorator_errors(self):
        def FocalOperation(self, locations, subdomains):
            pass
        def test():
            return executor(FocalOperation)
        self.assertRaises(PCMLOperationError, test)

        ZonalOperation = Operation(name='Test', opclass=OpClass.zonalclass)
        ZonalOperation._PCML_exported = True
        def test2():
            return executor(ZonalOperation)
        self.assertRaises(PCMLOperationError, test2)
