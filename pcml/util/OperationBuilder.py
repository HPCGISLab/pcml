"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from ..core.Operation import *
from ..core.Scheduler import *
from ..util.Messaging import PCMLOperationError
import types


class OperationDecorator:
    def __init__(self, opclass=OpClass.localclass, override='function', **kwargs):
        """ When user write "@operation(opclass=OpClass.focalclass, override='exectuor')",
        this function will be called first.

        :param opclass: Operation type, supported types are in OpClass .

        :param override: The override instance method in Operation class that will be overrided.
        Currently we support 'function' and 'executor', but the override function can be anything.

        :param kwargs: Extra mappings. Users may want to override other methods of Operation class,
        e.g.:
        > @operation(opclass=OpClass.localclass, getOutputLayers=_getOutputLayers)
        then the *instance method* getOutputLayers() will be replaced by _getOutputLayers.
        Using this feature, users can create customized Operation without writting a whole Operation class
        definition.
        Notice: Overriding class methods and static methods may cause problems.
        If the override function(either 'function' or 'executor') is overrided using kwargs, undefined behavior may happen, the order of overriding
        is not guaranteed and it is subject to change in future versions.
        :return: OperationDecorator object
        """
        self.opclass = opclass
        self.override = override
        self.mapping = kwargs

    def __call__(self, func):
        """ When users use @operation, this method will be called, user defined function will be passed as argument func.
        """
        def _func(*layers, **kwargs):
            """ This is the function that is actually exposed to users. kwargs are extra arguments, they will pass down
            to Operation.__init__, so one can write 'buffersize=20', and the Operation will construct an operation object
            with buffersize of 20.
            This function will be renamed to func.__name__, and will be added into current module's global name space.
            :param layers: Layer objects to process.
            :param kwargs: A list of arguments passed down to Operation.__init__
            """
            opclass = getattr(_func, 'opclass', None) or self.opclass
            override = getattr(_func, 'override', None) or self.override
            op = Operation(func.__name__, opclass=opclass, layers=layers, **kwargs)
            # Replace the override function
            setattr(op, override, types.MethodType(func, op, Operation))
            # Traverse the method mapping dictionary and replace them.
            for t, m in self.mapping.iteritems():
                setattr(op, t, type(m)((m, op, Operation)))
            return scheduler(op)
        # Mark _func as the function created by OperationDecorator
        _func._PCML_exported = True
        # Rename _func and add it to the current global name space
        _func.__name__ = func.__name__
        globals()[_func.__name__] = _func
        return _func

operation = OperationDecorator

# Helper decorators


def globaloperation(fn):
    if getattr(fn, '_PCML_exported', False):
        exported_func = fn
        exported_func.opclass = OpClass.globalclass
    else:
        exported_func = operation(opclass=OpClass.globalclass)(fn)
    return exported_func


def localoperation(fn):
    if getattr(fn, '_PCML_exported', False):
        exported_func = fn
        exported_func.opclass = OpClass.localclass
    else:
        exported_func = operation(opclass=OpClass.localclass)(fn)
    return exported_func


def focaloperation(fn):
    if getattr(fn, '_PCML_exported', False):
        exported_func = fn
        exported_func.opclass = OpClass.focalclass
    else:
        exported_func = operation(opclass=OpClass.focalclass)(fn)
    return exported_func


def zonaloperation(fn):
    if getattr(fn, '_PCML_exported', False):
        exported_func = fn
        exported_func.opclass = OpClass.zonalclass
    else:
        exported_func = operation(opclass=OpClass.zonalclass)(fn)
    return exported_func


def executor(fn):
    if not isinstance(fn, types.FunctionType):
        raise PCMLOperationError("Decorator @executor can only be used on functions")
    if not getattr(fn, '_PCML_exported', False):
        raise PCMLOperationError("Function %s should not be decorated directed by @executor" % fn.__name__)
    fn.override = 'executor'
    return fn
