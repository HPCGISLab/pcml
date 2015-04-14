"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from pcml import *
from pcml.util.LayerBuilder import lst_to_layer
from numpy.ma import allequal
import numpy as np
import unittest

class TestBoundingBox(unittest.TestCase):
    def setUp(self):
        # Basic bounding box from origin to h X w
        self.bb1=BoundingBox(0,0,100,100)
        # Bounding box that overlaps the origin including negative coordinates
        self.bb2=BoundingBox(-10,-10,20,20)
        # Bounding box that is entirely negative
        self.bb3=BoundingBox(-100,-100,10,10)
        # Bounding box that is entirely positive 
        self.bb4=BoundingBox(100,100,10,10)

        # Bounding box that has floating point y,x,h,w 
        self.bb5=BoundingBox(-1.2,-3.4,5.6,7.8)

    # Test the return from repr for example bounding boxes
    def test_boundingbox_repr(self):
        bb1_test = "<BoundingBox: (%f,%f) [%f,%f]>" % (0,0,100,100)
        bb2_test = "<BoundingBox: (%f,%f) [%f,%f]>" % (-10,-10,20,20)
        bb3_test = "<BoundingBox: (%f,%f) [%f,%f]>" % (-100,-100,10,10)
        bb4_test = "<BoundingBox: (%f,%f) [%f,%f]>" % (100,100,10,10)
        bb5_test = "<BoundingBox: (%f,%f) [%f,%f]>" % (-1.2,-3.4,5.6,7.8)
        self.assertEqual(bb1_test,repr(self.bb1))
        self.assertEqual(bb2_test,repr(self.bb2))
        self.assertEqual(bb3_test,repr(self.bb3))
        self.assertEqual(bb4_test,repr(self.bb4))
        self.assertEqual(bb5_test,repr(self.bb5))


    # Test the defaults for one bounding box - logic applies to all so only need to test one
    def test_boundingbox_defaults(self):
        self.assertEqual(self.bb5.y, -1.2)
        self.assertEqual(self.bb5.x, -3.4)
        self.assertEqual(self.bb5.h, 5.6)
        self.assertEqual(self.bb5.w, 7.8)
        self.assertEqual(self.bb5._data, None)
        self.assertEqual(self.bb5.data_structure, Datastructure.array)
        self.assertEqual(self.bb5.data_type, None)
        self.assertEqual(self.bb5.nodata_value, None)
        self.assertEqual(self.bb5.cellsize, None)
        self.assertEqual(self.bb5.nrows, None)
        self.assertEqual(self.bb5.ncols, None)
       
    # Ensure no negative or zero heights or widths are accepted 
    def test_boundingbox_initialization(self):
        with self.assertRaises(PCMLInvalidInput):
            failedbb=BoundingBox(0,0,-10,1)
        with self.assertRaises(PCMLInvalidInput):
            failedbb=BoundingBox(0,0,1,-10)
        with self.assertRaises(PCMLInvalidInput):
            failedbb=BoundingBox(0,0,0,1)
        with self.assertRaises(PCMLInvalidInput):
            failedbb=BoundingBox(0,0,1,0)

    def test_boundingbox_set_nparray(self):
        arr = np.asarray([[5] * 10] * 10) # Create a 10x10 array, each element has a value of 5
        self.bb4.set_nparray(arr,1,-10)
        #print arr
        #b1.set_nparray(arr, 1, 0)
        self.assertTrue(allequal(self.bb4._data, arr))
        self.assertEqual(self.bb4.cellsize,1)
        self.assertEqual(self.bb4.nodata_value,-10)
        self.assertEqual(self.bb4.data_structure,Datastructure.array)

        # Negative cellsize
        with self.assertRaises(PCMLInvalidInput):
            self.bb3.set_nparray(arr,-1,-10)
        # Zero cellsize
        with self.assertRaises(PCMLInvalidInput):
            self.bb3.set_nparray(arr,0,-10)
        # None cellsize
        with self.assertRaises(PCMLInvalidInput):
            self.bb3.set_nparray(arr,None,-10)
        # None array 
        nonearray=None
        with self.assertRaises(PCMLInvalidInput):
            self.bb2.set_nparray(nonearray,1,-10)

    def test_boundingbox_get_nparray(self):
        arr = np.asarray([[5] * 10] * 10) # Create a 10x10 array, each element has a value of 5
        self.bb4.set_nparray(arr,1,-10)
        self.assertTrue(allequal(self.bb4.get_nparray(), arr))
        self.assertIsNone(self.bb1.get_nparray())
    
     


    # FIXME: These are tests that we need
    '''
    check data_type when inputting data
    check that nparr is an array (np check)
    pointlist checks
    _reset_dim 
    set_data_ref
    get_locval
    get_ind_from_loc
    get xyloc
    bufferedlocgetarr
    slide_nparray
    print_data
    '''
        

