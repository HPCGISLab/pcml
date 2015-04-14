"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
import multiprocessing as mp
import numpy as np
import ctypes


# FROM : https://bitbucket.org/cleemesser/numpy-sharedmem/src/3d6dca4ffd926598c68faa3505df8ff0708989dd/doc/horesh_ctypes_sharedmem.py?at=default
# Another reference : https://chromium.googlesource.com/chromium/deps/python_26/+/f4f5d43a2599abd4ace5dead9d4d0e63e23055b7/Lib/multiprocessing/sharedctypes.py
_ctypes_to_numpy = {
    ctypes.c_char   : np.uint8,
    ctypes.c_wchar  : np.int16,
    ctypes.c_byte   : np.int8,
    ctypes.c_ubyte  : np.uint8,
    ctypes.c_short  : np.int16,
    ctypes.c_ushort : np.uint16,
    ctypes.c_int    : np.int32,
    ctypes.c_uint   : np.uint32,
    ctypes.c_long   : np.int64,
    ctypes.c_ulong  : np.uint64,
    ctypes.c_float  : np.float32,
    ctypes.c_double : np.float64}


def shmem_as_ndarray(raw_array, shape=None ):

    #address = raw_array._obj._wrapper.get_address()
    address = ctypes.addressof(raw_array)
    size = ctypes.sizeof(raw_array)
    size = len(raw_array)
    if (shape is None) or (np.asarray(shape).prod() != size):
        shape = (size,)
    elif type(shape) is int:
        shape = (shape,)
    else:
        shape = tuple(shape)

    #print(ctypes.typeof(raw_array))
    #dtype = _ctypes_to_numpy[raw_array._obj._type_]
    dtype = _ctypes_to_numpy[raw_array._type_]
    class Dummy(object): pass
    d = Dummy()
    d.__array_interface__ = {
        'data' : (address, False),
        'typestr' : np.dtype(dtype).str,
        'descr' : np.dtype(dtype).descr,
        'shape' : shape,
        'strides' : None,
        'version' : 3}
    return np.asarray(d).view( dtype=dtype )
# END FROM

def shmem_slice_nparray(arr,r,c,h,w):
    sliced_array=arr[r:r + h, c:c + w]
    return sliced_array

