Parallel Cartographic Modeling Language - PCML
==============================================


Introduction
------------

The parallel cartographic modeling language (PCML) is a multi-institutional 
collaborative project aiming to create a computing language for 
cyberGIScientists that is designed for (1) usability, (2) programmability, and 
(3) scalability. PCML provides multi-core parallel processing for spatial 
operations while hiding the implementation complexities of parallelism. 

Example
-------

PCML is designed to be simple to use for cartographic modeling.  The example below demonstrates the usability of PCML.  First users import the PCML files.  The second line configures PCML to use 4 processing cores for parallel processing, which is the default number.  Next the users will read an ASCII grid file and a GeoTIFF file as layer1 and layer2, respectively.  The two layers are added together using a LocalSum raster operation, which generates an output layers (layer_out).


    from pcml import *
    PCMLConfig.num_procs = 4 # Run computation in 4 processes (default)
    layer1 = ReadASCIIGrid("layer1.asc")
    layer2 = ReadGeoTIFF("layer2.tiff")
    layer_out = layer1 + layer2
    layer_out.print_data()


Please also see `test.py` for additional examples.



Installation
------------

### 1. Make sure that `pip` and `setuptools` are installed in your current Python environment (global or a virtualenv).

### 2. Install GDAL library

    $ sudo apt-get install libgdal-dev

or

    $ su -c 'yum install gdal-devel gdal-libs'

*The command may vary according to package manager and system.*

### 3. Install the required dependencies

    $ pip install -r requirements.txt

### 4. Finally, install

    $ python setup.py install

<!-- TODO: platform/distribution specific troubleshooting. -->



Build status
------------

[![wercker status](https://app.wercker.com/status/99dd16339b190c2ab04db505fa7af57a/m "wercker status")](https://app.wercker.com/project/bykey/99dd16339b190c2ab04db505fa7af57a)

