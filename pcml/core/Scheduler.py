"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Zhengliang Feng (odayfans@gmail.com, zfeng2@kent.edu)
"""
from ..util.Messaging import *
from .PCMLPrims import *
import pcml.core.PCMLConfig as PCMLConfig

import numpy as np
import multiprocessing as mp
import sys
import time

class PoolProcess( mp.Process ):
    def __init__(self, rank, numproc, lock, queue, subdomainlists,operation):
        mp.Process.__init__(self)
        self.rank = rank
        self.numproc = numproc
        self.lock = lock
        self.queue=queue
        self.subdomainlists=subdomainlists
        self.operation=operation
    def run(self):
        loop=True
        while loop:
            if self.queue.empty():
                loop=False
                continue
            try:
                subdomainsindex=self.queue.get(block=True, timeout=1)
            except:
                loop=False
                continue
            try:
                subdomains=self.subdomainlists[subdomainsindex]
                self.operation.executor(subdomains)
            except:
                PCMLException("Exception in process %i message %s " % (self.rank,sys.exc_info()[0]))

# This is a work-in-progress that will handle applying a function to a layer
# by dividing the layer into subdomains and applying the operation function to each location in a subdomain
# the function will return a value that will need to be saved in a new subdomain
# then we will need to combine the new subdomains back into a new layer or new layers to be returned
def scheduler(op):
    print("scheduling operation for execution %s"%op)

    # First decompose layers into multiple subdomains based on operation
    subdomainlists=op._decompositionrun()

    exectype=PCMLConfig.exectype
    #boolean array to check whether subdomains are processed
    processed=np.zeros(len(subdomainlists),dtype=bool)
    if exectype==ExecutorType.serialpython:
        print("Executing in serial python")

        # Call the executor for each group of subdomains in the subdomainlists
        for subdomains in subdomainlists:
            op.executor(subdomains)
        #iterate till all subdomains are processed
        while True:
            for i in xrange(len(subdomainlists)):
                if subdomainlists[i][0].get_itercount()==0:
                    processed[i]=True
            #check if all subdomains are processed
            if np.all(processed):
                break
            else:
                for i in xrange(len(subdomainlists)):
                    if processed[i]==True:
                        op.executor(subdomainlists[i])

    elif exectype==ExecutorType.parallelpythonqueue: # Parallel python version
        print("Executing in parallel python (Queue)")

        queue=mp.Queue()
        for i in xrange(len(subdomainlists)):
            queue.put(i)
            time.sleep(0.05) # Inserted to solve a bug in Python that is acknowledged online
            # https://pythonadventures.wordpress.com/tag/processes/
            # http://stackoverflow.com/questions/10607553/python-multiprocessing-queue-what-to-do-when-the-receiving-process-quits

        lock = mp.Lock()

        num_procs=PCMLConfig.num_procs

        print("Starting",num_procs,"processes to apply",op,"to",len(subdomainlists),"subdomains")
        pool = [PoolProcess(rank, num_procs, lock, queue, subdomainlists,op) for rank in range(num_procs)]
        for p in pool: p.start()
        for p in pool: p.join()
        #iterate till all subdomains are processed
        while True:
            for i in xrange(len(subdomainlists)):
                if subdomainlists[i][0].get_itercount()==0:
                    processed[i]=True
            if np.all(processed):
                break
            else:
                for i in xrange(len(subdomainlists)):
                    if processed[i]==True:
                        queue.put(i)
                pool = [PoolProcess(rank, num_procs, lock, queue, subdomainlists,op) for rank in range(num_procs)]
                for p in pool: p.start()
                for p in pool: p.join()
        if subdomainlists[0][0].data_structure==Datastructure.pointlist and not op.isnonlayeroutput:
            op.writepointdatatooutputlayer(subdomainlists)
    else:
        PCMLNotSupported("Scheduler does not support this exectype -"+str(exectype))
    return op.getOutputLayers()
