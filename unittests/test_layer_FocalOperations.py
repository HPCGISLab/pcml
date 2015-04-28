"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Suman Jindam (sjindam@kent.edu)
"""
from pcml import *
from pcml.util.LayerBuilder import *
from numpy.ma import allequal
from os import path
import numpy as np
import cml_test
import unittest

# Preliminary tests for spatial operations
class TestLayerOperationsSerial(cml_test.PCMLSerialTestCase):
    def setUp(self):
        super(TestLayerOperationsSerial, self).setUp()
        self.datadir = './data'
        self.l1 = lst_to_layer([[1]*4]*4)
        self.l2 = lst_to_layer([[2]*4]*4)
        self.l3 = lst_to_layer([[5]*4]*4)
        self.l4 = lst_to_layer([[normolized_value(1.53)] * 13] * 9)
        # l5 = l1+(l2+l3)*l4
        self.l5 = lst_to_layer([[normolized_value(11.71)] * 4] * 4)

        self.l6 = lst_to_layer([range(1,5)] * 4)
        self.l7 = lst_to_layer([[2, 2.5, 2.5, 3]] * 4)
        self.l8 = ReadASCIIGrid(path.join(self.datadir, 'data_c.asc'))
        self.l9 = Layer(0,0,100, 100, 'notitle')
        self.l9.set_nparray(np.ones((100,100)),.5,-999)
        self.l10 = Layer(0,0,100,100,'notitle')
        self.l11 = Layer(0,0,100, 100, 'notitle')
        pt_lst = [{'x':-81.4479691,'y':41.0593074,'z':1}
                 ,{'x':-81.5135,'y':41.0293074,'z':1}
                 ,{'x':-81.4779691,'y':41.0503074,'z':1}
                 ,{'x':-81.3779691,'y':41.0303074,'z':1}
                 ,{'x':-81.409691,'y':41.103074,'z':1}
                 ,{'x':-81.51079691,'y':41.08893074,'z':1}
                 ,{'x':-81.4779691,'y':41.0573074,'z':1}]
        self.l10.set_pointlist(pt_lst)
        self.l10.nrows = self.l9.nrows
        self.l10.ncols = self.l9.ncols
        arr11=np.array([[1,2,3,1],[1,2,3,2],[1,3,2,4],[1,3,2,1]])
        self.l11.set_nparray(arr11,5,-999)

        # To ensure FocalMean Operation gives the correct output with different layers
    def test_focalmean(self):
        lo = FocalMean(self.l1, buffersize=1)
        self.assertTrue(allequal(lo._data, self.l1._data), "FocalMean validation failed")
        lo = FocalMean(self.l4, buffersize=1)
        self.assertTrue(np.allclose(lo._data, self.l4._data))
        # To ensure FocalMean Operation gives the correct output with different buffersizes
        lo = FocalMean(self.l6, buffersize=2)
        self.assertTrue(np.allclose(lo._data, self.l7._data))

        # To ensure FocalMean Columndecompostion gives the correct output with different buffer sizes
    def test_focalmean_coldecomp(self):
        lo = FocalMean(self.l1, buffersize=1,decomposition=columndecomposition)
        self.assertTrue(allequal(lo._data, self.l1._data), "FocalMean validation failed")
        lo = FocalMean(self.l4, buffersize=1,decomposition=columndecomposition)
        self.assertTrue(np.allclose(lo._data, self.l4._data))
        lo = FocalMean(self.l6, buffersize=2,decomposition=columndecomposition)
        self.assertTrue(np.allclose(lo._data, self.l7._data))

        # To ensure FocalMean Operation with numpy implementation gives the correct output with different buffer sizes
    def test_focalmean_np(self):
        lo = FocalMean_np(self.l1, buffersize=1)
        self.assertTrue(allequal(lo._data, self.l1._data), "FocalMean_np validation failed")
        lo = FocalMean_np(self.l4, buffersize=1)
        self.assertTrue(np.allclose(lo._data, self.l4._data))
        lo = FocalMean_np(self.l6, buffersize=2)
        self.assertTrue(np.allclose(lo._data, self.l7._data))
        
         # To ensure FocalMean Operation with numpy  gives the correct output with different buffer sizes
    def test_focalmean_np(self):
        lo = FocalMean_np(self.l1, buffersize=1)
        self.assertTrue(allequal(lo._data, self.l1._data), "FocalMean_np validation failed")
        lo = FocalMean_np(self.l4, buffersize=1)
        self.assertTrue(np.allclose(lo._data, self.l4._data))
        lo = FocalMean_np(self.l6, buffersize=2)
        self.assertTrue(np.allclose(lo._data, self.l7._data))

        # To ensure FocalMaximum  Operation gives the correct output with different layers
    def test_focalmaximum(self):
        lo = FocalMaximum(self.l1,self.l2, buffersize=0)
        self.assertTrue(allequal(lo._data, self.l2._data))
        lo = FocalMaximum(self.l1,self.l2, buffersize=2)
        self.assertTrue(allequal(lo._data, self.l2._data))
        lo = FocalMaximum(self.l1,self.l2, buffersize=2,decomposition=columndecomposition)
        self.assertTrue(allequal(lo._data, self.l2._data))

       # To ensure FocalMinimum Operation gives the correct output with different buffer sizes
    def test_focalminimum(self):
        lo = FocalMinimum(self.l1,self.l2, buffersize=0)
        self.assertTrue(allequal(lo._data, self.l1._data))
        lo = FocalMinimum(self.l1,self.l2, buffersize=2)
        self.assertTrue(allequal(lo._data, self.l1._data))
        lo1=FocalMaximum(self.l1,lo, buffersize=1)
        self.assertTrue(allequal(lo1._data,lo._data))
        lo1=FocalMaximum(self.l1,lo, buffersize=1,decomposition=columndecomposition)
        self.assertTrue(allequal(lo1._data,lo._data))

        # To ensure FocalMaximum Operation with numpy gives the correct output with different buffer sizes
    def test_focalmaximum_np(self):
        lo = FocalMaximum_np(self.l1,self.l2, buffersize=0)
        self.assertTrue(allequal(lo._data, self.l2._data))
        lo = FocalMaximum_np(self.l1,self.l2, buffersize=2)
        self.assertTrue(allequal(lo._data, self.l2._data))
        lo1=FocalMaximum_np(self.l1,lo, buffersize=1)
        self.assertTrue(allequal(lo1._data,lo._data))
        lo1=FocalMaximum(self.l1,lo, buffersize=1,decomposition=columndecomposition)
        self.assertTrue(allequal(lo1._data,lo._data))
        
        # To ensure FocalMinimum Operation with numpy gives the correct output with different buffer sizes
    def test_focalminimum_np(self):
        lo = FocalMinimum_np(self.l1,self.l2, buffersize=0)
        self.assertTrue(allequal(lo._data, self.l1._data))
        lo = FocalMinimum_np(self.l1,self.l2, buffersize=2)
        self.assertTrue(allequal(lo._data, self.l1._data))
        lo1=FocalMinimum_np(self.l1,lo, buffersize=1)
        self.assertTrue(allequal(lo1._data,lo._data))
        lo1=FocalMaximum(self.l1,lo, buffersize=1,decomposition=columndecomposition)
        self.assertTrue(allequal(lo1._data,lo._data))

        # To ensure FocalMajority Operation gives the correct output with different buffer sizes
    def test_focalMajority(self):
        lo = FocalMajority(self.l1, buffersize=0)
        self.assertTrue(allequal(lo._data, self.l1._data))
        lo = FocalMajority(self.l1, buffersize=1)
        self.assertTrue(allequal(lo._data, self.l1._data))
        res = np.asarray([[1,1,2,3],[1,1,2,3],[1,1,2,2],[1,1,3,2]])
        lo = FocalMajority(self.l11, buffersize=1)
        self.assertTrue(allequal(res,lo._data))
        lo = FocalMajority(self.l11, buffersize=1,decomposition=columndecomposition)
        self.assertTrue(allequal(res,lo._data))

         # To ensure FocalPercentage Operation gives the correct output with different buffer sizes
    def test_focalpercentage(self):
        lo = FocalPercentage(self.l1, buffersize=1)
        res = np.asarray([[100]*4]*4)
        self.assertTrue(allequal(lo._data, res))
        lo = FocalPercentage(self.l2, buffersize=3,decomposition=columndecomposition)
        self.assertTrue(allequal(lo._data, res))

         # To ensure FocalMean Operation with numpy by executor  gives the correct output with different buffer sizes
    def test_focalmean_np_exec(self):
        lo = FocalMean_np_exec(self.l1, buffersize=1)
        res = np.asarray([[1]*4]*4)
        self.assertTrue(allequal(lo._data, res))
        lo = FocalMean_np_exec(self.l2, buffersize=3,decomposition=columndecomposition)
        res = np.asarray([[2]*4]*4)
        self.assertTrue(allequal(lo._data, res))
        
        # To ensure FocalSum Operation  gives the correct output with different buffer sizes
    def test_focalsum(self):
        lo = FocalSum(self.l1, buffersize=1)
        res = np.asarray([[4,6,6,4],[6,9,9,6],[6,9,9,6],[4,6,6,4]])
        self.assertTrue(allequal(lo._data, res))
        lo = FocalSum(self.l1, buffersize=2)
        res = np.asarray([[9,12,12,9],[12,16,16,12],[12,16,16,12],[9,12,12,9]])
        self.assertTrue(allequal(lo._data, res))

class TestLayerOperationsParallel(TestLayerOperationsSerial):
    def setUp(self):
        super(TestLayerOperationsParallel, self).setUp()
        PCMLConfig.num_procs = 4
        PCMLConfig.exectype = ExecutorType.parallelpythonqueue

if __name__ == '__main__':
    unittest.main()

# FIXME: things to test
'''
    
'''
