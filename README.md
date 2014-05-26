Parallelized Cartographic Modeling Language - pCML
===

Introduction
---
The parallel cartographic modeling language (pCML) is a multi-institutional collaborative project 
aiming to create a computing language for cyberGIScientists that is designed for
(1) usability, (2) programmability, and (3) scalability. pCML provides multi-core parallel
processing for spatial operations while hiding the implementation complexities of parallelism.

Installation
---
1. Make sure that pip and setuptools are installed in your current Python environment (global or a virtualenv).
2. Install GDAL library

    $ sudo apt-get install libgdal-dev

Change the command accordingly to your system.

3. run

    $ pip install -r requirements.txt # Need to use pip to install numpy or the next line will fail
    $ python setup.py install

Example
---

    from pCML import *
    PCMLConfig.num_procs = 4 # Run computation in 4 processes (this is the default configuration)
    layer1 = ReadASCIIGrid("layer1.asc")
    layer2 = ReadGeoTIFF("layer2.tiff")
    layer_out = layer1 + layer2
    layer_out.print_data()


    Please also see test.py for additional working examples

Build status
---
[![wercker status](https://app.wercker.com/status/99dd16339b190c2ab04db505fa7af57a/m "wercker status")](https://app.wercker.com/project/bykey/99dd16339b190c2ab04db505fa7af57a)
