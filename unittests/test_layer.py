"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Suman Jindam (sjindam@kent.edu)
"""
from pcml import *
from pcml.util.Messaging import PCMLNotImplemented
from pcml.util.LayerBuilder import lst_to_layer
from numpy.ma import allequal
import numpy as np
import unittest

# TODO: Use the data*.asc files as test cases for unit + integration testing.
# They should be tested against gold standard data.

class TestLayer(unittest.TestCase):
    def setUp(self):
        # Basic layer from origin to h X w
        self.layer1=Layer(0,0,100,100,"Layer One")
        # Layer that overlaps the origin including negative coordinates
        self.layer2=Layer(-10,-10,20,20,"Layer Two")
        # Layer that is entirely negative
        self.layer3=Layer(-100,-100,10,10,"Layer Three")
        # Layer that is entirely positive 
        self.layer4=Layer(100,100,10,10,"Layer Four")

        # Layer that has floating point y,x,h,w 
        self.layer5=Layer(-1.2,-3.4,5.6,7.8,"Layer Five")
        # Layers that has predefined array values
        self.layer6=Layer(100,0,10,10,"Layer Six")
        arr111=np.asarray([[1]*3]*3)
        self.layer6.set_nparray(arr111,1,10)
        self.layer7=Layer(100,0,10,10,"Layer seven")
        arr=np.asarray([[5]*3]*3)
        self.layer7.set_nparray(arr,1,10)
        
    # Test the return from repr for example layers 
    def test_layer_repr(self):
        layer1_test = "<Layer: (%f,%f) [%f,%f] : %s>" % (0,0,100,100,"Layer One")
        layer2_test = "<Layer: (%f,%f) [%f,%f] : %s>" % (-10,-10,20,20,"Layer Two")
        layer3_test = "<Layer: (%f,%f) [%f,%f] : %s>" % (-100,-100,10,10,"Layer Three")
        layer4_test = "<Layer: (%f,%f) [%f,%f] : %s>" % (100,100,10,10,"Layer Four")
        layer5_test = "<Layer: (%f,%f) [%f,%f] : %s>" % (-1.2,-3.4,5.6,7.8,"Layer Five")
        self.assertEqual(layer1_test,repr(self.layer1))
        self.assertEqual(layer2_test,repr(self.layer2))
        self.assertEqual(layer3_test,repr(self.layer3))
        self.assertEqual(layer4_test,repr(self.layer4))
        self.assertEqual(layer5_test,repr(self.layer5))


    # Test the defaults for one layer - logic applies to all so only need to test one
    def test_layer_defaults(self):
        self.assertEqual(self.layer5.y, -1.2)
        self.assertEqual(self.layer5.x, -3.4)
        self.assertEqual(self.layer5.h, 5.6)
        self.assertEqual(self.layer5.w, 7.8)
        self.assertEqual(self.layer5.title, "Layer Five")
        self.assertEqual(self.layer5._data, None)
        self.assertEqual(self.layer5.data_structure, Datastructure.array)
        self.assertEqual(self.layer5.data_type, None)
        self.assertEqual(self.layer5.nodata_value, None)
        self.assertEqual(self.layer5.cellsize, None)
        self.assertEqual(self.layer5.nrows, None)
        self.assertEqual(self.layer5.ncols, None)
       
    # Ensure no negative or zero heights or widths are accepted 
    def test_layer_initialization(self):
        with self.assertRaises(PCMLInvalidInput):
            failedlayer=Layer(0,0,-10,1,'f')
        with self.assertRaises(PCMLInvalidInput):
            failedlayer=Layer(0,0,1,-10,'f')
        with self.assertRaises(PCMLInvalidInput):
            failedlayer=Layer(0,0,0,1,'f')
        with self.assertRaises(PCMLInvalidInput):
            failedlayer=Layer(0,0,1,0,'f')

    
    def test_layer_add(self):
    #  Ensure layer addition  gives the correct output
        layer_o=self.layer6+self.layer7
        array1=np.asarray(layer_o.get_nparray())
        res = np.asarray([[6]*3]*3)
        self.assertEqual(np.all(array1==6),True)
        self.assertTrue(allequal(layer_o._data, res))
   # FIXME:Above test should fail for different spatial dimensions

    def test_layer_sub(self):
    # Ensure layer subtraction  gives the correct output
        layer_o=self.layer6-self.layer7
        array1=np.asarray(layer_o.get_nparray())
        res = np.asarray([[-4]*3]*3)
        self.assertEqual(np.all(array1==-4),True)
        self.assertTrue(allequal(layer_o._data, res))
    # FIXME:Above test should fail for different spatial dimensions

    def test_layer_mul(self):
    # Ensure layer multiplication  gives the correct output
        layer_o=self.layer6*self.layer7
        array1=np.asarray(layer_o.get_nparray())
        res = np.asarray([[5]*3]*3)
        self.assertEqual(np.all(array1==5),True)
        self.assertTrue(allequal(layer_o._data, res))
     # FIXME:Above test should fail for different spatial dimensions

    def test_layer_div(self):
   # Ensure layer division  gives the correct output
        layer_o=self.layer6/self.layer7
        array1=np.asarray(layer_o.get_nparray())
        res = np.asarray([[0.2]*3]*3)
        self.assertEqual(np.all(array1==0.2),True)
        self.assertTrue(allequal(layer_o._data, res))
    # FIXME:Above test should fail for different spatial dimensions
    
    
    # To ensure duplicate layer also has same number of rows & columns and same x,y,h,w values
    def test_layer_duplicate(self):
        result=np.asarray([[0]*3]*3)
        newlayer=self.layer7.duplicate()
        array=np.asarray(newlayer.get_nparray())
        self.assertEqual(len(result),len(array))
        self.assertEqual(len(result[0]),len(array[0]))
        self.assertEqual(newlayer.y,self.layer7.y)
        self.assertEqual(newlayer.x,self.layer7.x)
        self.assertEqual(newlayer.h,self.layer7.h)
        self.assertEqual(newlayer.w,self.layer7.w)
        self.assertEqual(newlayer.cellsize,self.layer7.cellsize)
        self.assertEqual(newlayer.nodata_value,self.layer7.nodata_value)
        self.assertEqual(np.all(array==0),True)
        
    # FIXME: Tests to be written
    '''
    decomposition
    rowdecomposition
    '''
