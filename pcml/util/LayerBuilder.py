"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from pCML import PCMLConfig
from ..core.Layer import *
from numpy import asarray

def lst_to_layer(lst):
    row = len(lst)
    assert(row > 0)
    col = len(lst[0])
    assert(col > 0)
    l = Layer(0, 0, row, col, "test layer")
    l.set_nparray(asarray(lst), 1, 0)
    return l

def normolized_value(val):
    return float(PCMLConfig.value_precision % val)
