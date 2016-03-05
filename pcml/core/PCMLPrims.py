"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu)
"""


# Define data structures supported in PCML
class Datastructure():
    """ Enumeration class. Defines the type of data structure.
    """
    array = 1
    pointlist = 2


# Define the type of Executor to use in the Scheduler
class ExecutorType():
    """ Enumeration class. It defines which executor to apply for an operation
    """
    serialpython = 1
    parallelpythonqueue = 2


class OpClass():
    """ Enumeration class. Classifies operations as local, focal, zonal, or global.
    """
    localclass = 1
    focalclass = 2
    zonalclass = 3
    globalclass = 4
