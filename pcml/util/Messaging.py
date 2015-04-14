"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
class PCMLUserInformation(UserWarning):
    def __init__(self,value):
        self.value=" [ INFORMATION ] This is informational : %s" % value 
        print(self.value)
    def __str__(self):
        return repr(self.value)

class PCMLNotImplemented(UserWarning):
    def __init__(self,value):
        self.value=" [ WARNING ] This is not implemented : %s" % value 
        #print(self.value)
    def __str__(self):
        return repr(self.value)

class PCMLNotSupported(Exception):
    def __init__(self,value):
        self.value=" [ ERROR ] This is not supported : %s" % value 
        print(self.value)
        exit(1)

    def __str__(self):
        return repr(self.value)

class PCMLInvalidInput(Exception):
    def __init__(self,msg,value):
        self.value=" [ ERROR ] The input is invalid (%s) : %s" % (msg,value)
        print(self.value)
        #exit(1)

    def __str__(self):
        return repr(self.value)

class PCMLOperationError(Exception):
    def __init__(self,value):
        self.value=" [ ERROR ] Operation error : %s" % value
        print(self.value)
    def __str__(self):
        return repr(self.value)

class PCMLException(Exception):
    def __init__(self,value):
        self.value=" [ ERROR ] Exception : %s" % value 
        print(self.value)
    def __str__(self):
        return repr(self.value)


def PCMLTODO(msg):
    printtodo=False
    if printtodo:
        print(" [ TODO ] %s" % msg)


