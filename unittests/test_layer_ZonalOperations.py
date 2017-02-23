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

# TODO: Use the dataz.asc file to do unit + integration tests for zonal statistics.

# Preliminary tests for spatial operations
class TestLayerOperationsSerial(cml_test.PCMLSerialTestCase):
    def setUp(self):
        super(TestLayerOperationsSerial, self).setUp()
        self.datadir = './data'
        self.l1 = lst_to_layer([[1]*4]*4)
        self.l2 = lst_to_layer([[2]*4]*4)
        self.l3 = lst_to_layer([[5]*4]*4)
        self.l4 = lst_to_layer([[normolized_value(1.53)] * 13] * 9)
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

        #Layer with all the positive numbers
        self.l12 = lst_to_layer([[5,3,1,6],[4,2,2,1],[3,1,6,7],[8,9,6,1]])
        self.l13 = lst_to_layer([[1,1,1,2],[2,2,3,3],[3,4,4,4],[5,5,4,5]])
        #Layer with all the negative numbers
        self.l14 = lst_to_layer([[-5,-3,-1,-6],[-4,-2,-2,-1],[-3,-1,-6,-7],[-8,-9,-6,-1]])

        #To ensure ZonalSum operation gives the correct output with positive and negative zones
    def test_zonalsum(self):
        lo=zonalsum(self.l12,self.l13)
        res=np.asarray([[9,9,9,12],[12,12,6,6],[6,20,20,20],[18,18,20,18]])
        self.assertTrue(allequal(lo._data,res))
        #To check ZonalSum operation with negative zonal values
        lo=zonalsum(self.l14,self.l13)
        res=np.asarray([[-9,-9,-9,-12],[-12,-12,-6,-6],[-6,-20,-20,-20],[-18,-18,-20,-18]])
        self.assertTrue(allequal(lo._data,res))

        #To ensure ZonalMean operation gives the correct output with positive and negative zones
    def test_zonalmean(self):
        lo=zonalmean(self.l12,self.l13)
        res=np.asarray([[3,3,3,4],[4,4,2,2],[2,5,5,5],[6,6,5,6]])
        self.assertTrue(allequal(lo._data,res))
        #To check ZonalMean operation with negative zonal values
        lo=zonalmean(self.l14,self.l13)
        res=np.asarray([[-3,-3,-3,-4],[-4,-4,-2,-2],[-2,-5,-5,-5],[-6,-6,-5,-6]])
        self.assertTrue(allequal(lo._data,res))

        #To ensure ZonalMaximum operation gives the correct output with positive and negative zones
    def test_zonalmaximum(self):
        lo=zonalmaximum(self.l12,self.l13)
        res=np.asarray([[5,5,5,6],[6,6,3,3],[3,7,7,7],[9,9,7,9]])
        self.assertTrue(allequal(lo._data,res))
        #To check ZonalMaximum operation with negative zonal values
        lo=zonalmaximum(self.l14,self.l13)
        res=np.asarray([[-1,-1,-1,-2],[-2,-2,-1,-1],[-1,-1,-1,-1],[-1,-1,-1,-1]])
        self.assertTrue(allequal(lo._data,res))
        
    def test_zonalminimum(self):
        lo=zonalminimum(self.l12,self.l13)
        res=np.asarray([[1,1,1,2],[2,2,1,1],[1,1,1,1],[1,1,1,1]])
        self.assertTrue(allequal(lo._data,res))
        #To check ZonalMinimum operation with negative zonal values
        lo=zonalminimum(self.l14,self.l13)
        res=np.asarray([[-5,-5,-5,-6],[-6,-6,-3,-3],[-3,-7,-7,-7],[-9,-9,-7,-9]])
        self.assertTrue(allequal(lo._data,res))

        #To ensure ZonalMajority operation gives the correct output with positive and negative zones
    """
    def test_zonalmajority(self):
        lo=zonalmajority(self.l12,self.l13)
        res=np.asarray([[1,1,1,2],[2,2,1,1],[1,6,6,6],[1,1,6,1]])
        self.assertTrue(allequal(lo._data,res))
        l12 = lst_to_layer([[5,5,5,6],[4,2,2,2],[3,6,6,6],[9,9,9,1]])
        l13 = lst_to_layer([[1,1,1,1],[2,2,2,2],[3,3,3,3],[4,4,4,4]])
        lo=zonalmajority(l12,l13)
        res=np.asarray([[5,5,5,5],[2,2,2,2],[6,6,6,6],[9,9,9,9]])
        self.assertTrue(allequal(lo._data,res))
        
        #To ensure ZonalMinority operation gives the correct output with positive and negative zones
    def test_zonalminority(self):
        lo=zonalminority(self.l12,self.l13)
        res=np.asarray([[1,1,1,2],[2,2,1,1],[1,1,1,1],[1,1,1,1]])
        self.assertTrue(allequal(lo._data,res))
        l12 = lst_to_layer([[5,5,5,6],[4,2,2,2],[3,6,6,6],[9,9,9,1]])
        l13 = lst_to_layer([[1,1,1,1],[2,2,2,2],[3,3,3,3],[4,4,4,4]])
        lo=zonalminority(l12,l13)
        res=np.asarray([[6,6,6,6],[4,4,4,4],[3,3,3,3],[1,1,1,1]])
        self.assertTrue(allequal(lo._data,res))
     """
        #To ensure ZonalSum operation by executor gives the correct output with positive and negative values
    def test_ZonalSum_exec(self):
        lo=ZonalSum_exec(self.l12,self.l13)
        res=np.asarray([[9,9,9,12],[12,12,6,6],[6,20,20,20],[18,18,20,18]])
        self.assertTrue(allequal(lo._data,res))
        lo=ZonalSum_exec(self.l14,self.l13)
        #To ensure ZonalSum operation by executor gives the correct output with negative values
        res=np.asarray([[-9,-9,-9,-12],[-12,-12,-6,-6],[-6,-20,-20,-20],[-18,-18,-20,-18]])
        self.assertTrue(allequal(lo._data,res))

        #To ensure ZonalMean operation by executor gives the correct output with positive and negative values
    def test_ZonalMean_exec(self):
        lo=ZonalMean_exec(self.l12,self.l13)
        res=np.asarray([[3,3,3,4],[4,4,2,2],[2,5,5,5],[6,6,5,6]])
        self.assertTrue(allequal(lo._data,res))
        lo=ZonalMean_exec(self.l14,self.l13)
        res=np.asarray([[-3,-3,-3,-4],[-4,-4,-2,-2],[-2,-5,-5,-5],[-6,-6,-5,-6]])
        self.assertTrue(allequal(lo._data,res))

        #To ensure ZonalMaximum operation by executor gives the correct output with positive and negative values
    def test_ZonalMaximum_exec(self):
        lo=ZonalMaximum_exec(self.l12,self.l13)
        res=np.asarray([[5,5,5,6],[6,6,3,3],[3,7,7,7],[9,9,7,9]])
        self.assertTrue(allequal(lo._data,res))
        lo=ZonalMaximum_exec(self.l14,self.l13)
        res=np.asarray([[-1,-1,-1,-2],[-2,-2,-1,-1],[-1,-1,-1,-1],[-1,-1,-1,-1]])
        self.assertTrue(allequal(lo._data,res))

        #To ensure ZonalMinimum operation by executor gives the correct output with positive and negative values
    def test_ZonalMinimum_exec(self):
        lo=ZonalMinimum_exec(self.l12,self.l13)
        res=np.asarray([[1,1,1,2],[2,2,1,1],[1,1,1,1],[1,1,1,1]])
        self.assertTrue(allequal(lo._data,res))
        lo=ZonalMinimum_exec(self.l14,self.l13)
        res=np.asarray([[-5,-5,-5,-6],[-6,-6,-3,-3],[-3,-7,-7,-7],[-9,-9,-7,-9]])
        self.assertTrue(allequal(lo._data,res))
        
        #To ensure ZonalMajority operation by executor gives the correct output with positive and negative values
    def test_ZonalMajority_exec(self):
        lo=ZonalMajority_exec(self.l12,self.l13)
        res=np.asarray([[1,1,1,2],[2,2,1,1],[1,6,6,6],[1,1,6,1]])
        self.assertTrue(allequal(lo._data,res))
        l12 = lst_to_layer([[5,5,5,6],[4,2,2,2],[3,6,6,6],[9,9,9,1]])
        l13 = lst_to_layer([[1,1,1,1],[2,2,2,2],[3,3,3,3],[4,4,4,4]])
        lo=ZonalMajority_exec(l12,l13)
        res=np.asarray([[5,5,5,5],[2,2,2,2],[6,6,6,6],[9,9,9,9]])
        self.assertTrue(allequal(lo._data,res))

        #To ensure ZonalMinority operation by executor gives the correct output with positive values
    def test_ZonalMinority_exec(self):
        lo=ZonalMinority_exec(self.l12,self.l13)
        res=np.asarray([[1,1,1,2],[2,2,1,1],[1,1,1,1],[1,1,1,1]])
        self.assertTrue(allequal(lo._data,res))
        l12 = lst_to_layer([[5,5,5,6],[4,2,2,2],[3,6,6,6],[9,9,9,1]])
        l13 = lst_to_layer([[1,1,1,1],[2,2,2,2],[3,3,3,3],[4,4,4,4]])
        lo=ZonalMinority_exec(l12,l13)
        res=np.asarray([[6,6,6,6],[4,4,4,4],[3,3,3,3],[1,1,1,1]])
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
'''
