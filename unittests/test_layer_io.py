"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from pCML import *
from pCML.util.LayerBuilder import *
import numpy as np
from numpy.ma import allequal
import unittest
from os import path, remove

class TestLayerIO(unittest.TestCase):
    def setUp(self):
        self.datadir='./data'

    def test_asciigrid_read_write(self):
        print "Testing asc IO 1..."
        test_file = path.join(self.datadir, 'test_layerio_asc.asc')
        l1 = lst_to_layer([[1] * 4]*4)
        print 'l1='
        print l1._data
        WriteASCIIGrid(test_file, l1)
        l2 = ReadASCIIGrid(test_file)
        print 'l2='
        print l2._data
        self.assertTrue(allequal(l1._data, l2._data))

        print "Testing asc IO 2..."
        value = normolized_value(1.532)
        l3 = lst_to_layer([[value] * 10] * 9)
        print 'l3='
        print l3._data
        WriteASCIIGrid(test_file, l3)
        l4 = ReadASCIIGrid(test_file)
        print 'l4='
        print l4._data
        self.assertTrue(allequal(l3._data, l4._data))
        remove(test_file)

    def test_tiff_read_write(self):
        print "Testing tiff IO 1..."
        test_file = path.join(self.datadir, 'test_layer_tiff.tif')
        l1 = lst_to_layer([[1]*4]*4)
        WriteGeoTIFF(test_file, l1)
        print 'l1='
        print l1._data
        l2 = ReadGeoTIFF(test_file)
        print 'l2='
        print l2._data
        self.assertTrue(allequal(l1._data, l2._data))

        print "Testing tiff IO 2..."
        value = normolized_value(92.123)
        l3 = lst_to_layer([[value] * 5] * 10)
        print 'l3='
        print l3._data
        WriteGeoTIFF(test_file, l3)
        l4 = ReadGeoTIFF(test_file)
        print 'l4='
        print l4._data
        self.assertTrue(allequal(l3._data, l4._data))
        remove(test_file)

    def test_file_not_found(self):
        test_asc = path.join(self.datadir, 'does/not/exist/file.asc')
        test_tiff = path.join(self.datadir, 'does/not/exist/file.tif')
        exceptions = (PCMLException, IOError) # expected exceptions

        with self.assertRaises(exceptions):
            l1 = ReadASCIIGrid(test_asc)
        with self.assertRaises(exceptions):
            # This will raise an exception because the directory is not accessable
            l1 = lst_to_layer([[1] * 8] * 8)
            WriteASCIIGrid(test_asc, l1)

        with self.assertRaises(exceptions):
            l1 = ReadGeoTIFF(test_tiff)
        with self.assertRaises(exceptions):
            l1 = lst_to_layer([[2.1] * 9] * 3)
            WriteGeoTIFF(test_tiff, l1)

    def test_tiff_illegal_band_number(self):
        test_file = path.join(self.datadir, 'test_tiff.tif')
        l1 = lst_to_layer([[10] * 5] * 8)
        WriteGeoTIFF(test_file, l1)
        with self.assertRaises(PCMLException):
            l2 = ReadGeoTIFF(test_file, bandnumber = 10)
        remove(test_file)

