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

    def test_LocalSum(self):
        lo = LocalSum(self.l1, self.l2)
        res = np.asarray([[3]*4]*4)
        self.assertTrue(allequal(lo._data, res))

    def test_LocalSum_coldecomp(self):
        lo = LocalSum(self.l1, self.l2,decomposition=columndecomposition)
        res = np.asarray([[3]*4]*4)
        self.assertTrue(allequal(lo._data, res))

    def test_LocalSum_np(self):
        lo = LocalSum_np(self.l1, self.l2)
        res = np.asarray([[3]*4]*4)
        self.assertTrue(allequal(lo._data, res))

    def test_LocalMult(self):
        lo = LocalMult(self.l1, self.l2)
        self.assertTrue(allequal(lo._data, self.l2._data))
        lo = LocalMult(self.l2, self.l3)
        self.assertTrue(allequal(lo._data, [[10]*4]*4))

    def test_LocalMult_np(self):
        lo = LocalMult_np(self.l1, self.l2)
        self.assertTrue(allequal(lo._data, self.l2._data))
        lo = LocalMult_np(self.l2, self.l3)
        self.assertTrue(allequal(lo._data, [[10]*4]*4))

    def test_add(self):
        lo = self.l1 + self.l2
        self.assertTrue(allequal(lo._data, [[3]*4]*4))

    def test_mult(self):
        lo = self.l1 * self.l2
        self.assertTrue(allequal(lo._data, self.l2._data))
        lo = self.l2 * self.l3
        self.assertTrue(allequal(lo._data, [[10]*4]*4))

    def test_mult_add(self):
        lo = self.l1 + (self.l2 + self.l3) * self.l4
        self.assertTrue(np.allclose(lo._data, self.l5._data))

    def test_local_maximum_np(self):
        lo = LocalMaximum_np(self.l1,self.l2)
        res = np.asarray([[2]*4]*4)
        self.assertTrue(allequal(lo._data, res))

    def test_local_minimum_np(self):
        lo = LocalMinimum_np(self.l1,self.l2)
        res = np.asarray([[1]*4]*4)
        self.assertTrue(allequal(lo._data, res))

    def test_localsum_np(self):
        lo = LocalSum_np(self.l1,self.l2)
        res = np.asarray([[3]*4]*4)
        self.assertTrue(allequal(lo._data, res))

    def test_localmult_np(self):
        lo = LocalMult_np(self.l1,self.l2)
        res = np.asarray([[2]*4]*4)
        self.assertTrue(allequal(lo._data, res))

    def test_localmean_np(self):
        lo = LocalMean_np(self.l1,self.l2)
        res = np.asarray([[1.5]*4]*4)
        self.assertTrue(allequal(lo._data, res))

    def test_localproduct_np(self):
        lo = LocalProduct_np(self.l1,self.l2)
        res = np.asarray([[2]*4]*4)
        self.assertTrue(allequal(lo._data, res))

    def test_localdivision(self):
        lo = LocalDivision(self.l1,self.l2)
        res = np.asarray([[0.5]*4]*4)
        self.assertTrue(allequal(lo._data, res))
    

    def test_localsubtraction(self):
        lo = LocalSubtraction(self.l1,self.l2)
        res = np.asarray([[-1]*4]*4)
        self.assertTrue(allequal(lo._data, res))

    def test_local_sin(self):
        lo=LocalSin(self.l1)
        res = np.asarray([[0.84147098]*4]*4)
        boolean =np.allclose(res,lo.get_nparray())
        self.assertEqual(boolean,True)

    def test_local_tangent(self):
        lo=LocalTangent(self.l1)
        res = np.asarray([[1.55740772]*4]*4)
        boolean =np.allclose(res,lo.get_nparray())
        self.assertEqual(boolean,True)


    def test_local_maximum(self):
        l22 = lst_to_layer([[5]*4]*4)
        lo=LocalMaximum(self.l1,self.l2,l22)
        res = np.asarray([[5]*4]*4)
        self.assertTrue(allequal(lo._data,res))

    def test_local_minimum(self):
        l22 = lst_to_layer([[5]*4]*4)
        lo=LocalMinimum(self.l2,self.l1,l22)
        res = np.asarray([[1]*4]*4)
        self.assertTrue(allequal(lo._data,res))

    def test_localmean(self):
        l22 = lst_to_layer([[3]*4]*4)
        lo=LocalMean(self.l1,self.l2,l22)
        res = np.asarray([[2]*4]*4)
        self.assertTrue(allequal(lo._data,res))
  
    def test_localVariety(self):
        l22 = lst_to_layer([[3]*4]*4)
        lo=LocalVariety(self.l1,self.l2,l22)
        res = np.asarray([[3]*4]*4)
        self.assertTrue(allequal(lo._data,res))

    def test_localDifference(self):
        l22 = lst_to_layer([[3]*4]*4)
        lo=LocalDifference(self.l1,self.l2,l22)
        res = np.asarray([[-4]*4]*4)
        self.assertTrue(allequal(lo._data,res))   
        
    def test_localDifference_np(self):
        l22 = lst_to_layer([[3]*4]*4)
        lo=LocalDifference(self.l1,self.l2,l22)
        res = np.asarray([[-4]*4]*4)
        
     def test_local_sine(self):
        lo=LocalSine(self.l1)
        res = np.asarray([[0.84147098]*4]*4)
        boolean =np.allclose(res,lo.get_nparray())
        self.assertEqual(boolean,True)

    def test_local_cosine(self):
        lo=LocalCosine(self.l1)
        res = np.asarray([[0.54030230]*4]*4)
        boolean =np.allclose(res,lo.get_nparray())
        self.assertEqual(boolean,True)

    def test_local_ArcSine(self):
        lo=LocalArcSine(self.l1)
        res = np.asarray([[1.57079633]*4]*4)
        boolean =np.allclose(res,lo.get_nparray())
        self.assertEqual(boolean,True)

    def test_local_ArcCosine(self):
        lo=LocalArcCosine(self.l1)
        res = np.asarray([[0]*4]*4)
        boolean =np.allclose(res,lo.get_nparray())
        self.assertEqual(boolean,True)

    def test_local_ArcTangent(self):
        lo=localArcTangent(self.l1)
        res = np.asarray([[0.78539816]*4]*4)
        boolean =np.allclose(res,lo.get_nparray())
        self.assertEqual(boolean,True)

    def test_LocalClassify20_np(self):
        lo=LocalClassify20_np(self.l1)
        res = np.asarray([[0]*4]*4)
        lo=LocalClassify20_np(self.l2)
        self.assertTrue(allequal(lo._data,res))
        self.assertTrue(allequal(lo._data,res))

    def test_LocalSum_numba(self):
        lo=LocalSum_numba(self.l1,self.l2)
        res = np.asarray([[3]*4]*4)
        self.assertTrue(allequal(lo._data,res))

 
class TestLayerOperationsParallel(TestLayerOperationsSerial):
    def setUp(self):
        super(TestLayerOperationsParallel, self).setUp()
        PCMLConfig.num_procs = 4
        PCMLConfig.exectype = ExecutorType.parallelpythonqueue

if __name__ == '__main__':
    unittest.main()


# FIXME: things to test
'''
    @localoperation has buffersize=0
    @focaloperation has buffersize>0
    @globaloperation has buffersize=-1


'''

