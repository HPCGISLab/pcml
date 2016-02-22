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

PCML is easy to use for parallel spatial data processing.  First, users import PCML.  The second line tells PCML to use 4 processing cores for parallel processing.  Next, the users will read an ASCII grid file and a GeoTIFF file as layer1 and layer2, respectively.  The two layers are added together using a LocalSum raster operation. The output output layer (layer_out) is printed to the screen.


    from pcml import *
    PCMLConfig.num_procs = 4 # Run computation in 4 processes (default)
    layer1 = ReadASCIIGrid("layer1.asc")
    layer2 = ReadGeoTIFF("layer2.tiff")
    layer_out = layer1 + layer2
    layer_out.print_data()


Please also see `test.py` or `examples.py` for additional examples.


Windows Installation
------------

#### 0. Download and install Ananconda (Optional)

   The Anaconda Python distribution package is available at http://continuum.io/downloads.
   It provides Python and several Python packages, but can be skipped if you have Python already.

#### 1. Download PCML
   Navigate to the PCML page on GitHub at https://github.com/hpcgislab/pcml

   Click "Download ZIP"

#### 2.	Unzip PCML

   A file named “PCML-master.zip” will be downloaded to your computer. Unzip it.

#### 3. Open ‘demo_30min.py’ 

   Open 'demo_30min.py' by clicking File->Open in the Spyder editor (included in the Ananconda package)

#### 4. Create your own LocalSum function in under 30 minutes!

#### 5. Run PCML and start processing spatial data in parallel
