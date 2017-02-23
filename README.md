Parallel Cartographic Modeling Language - PCML
==============================================


Introduction
------------

The parallel cartographic modeling language (PCML) is a multi-institutional 
collaborative project aiming to create a computing language for 
cyberGIScientists that is designed for (1) usability, (2) programmability, and 
(3) scalability. PCML provides multi-core parallel processing for spatial 
operations while hiding the implementation complexities of parallelism. 

**Note** PCML has been making progress behind-the-scenes.  One of the primary developers (Eric Shook), moved to the University of Minnesota, which disrupted some of the public-facing PCML development.  However, several students and Dr. Shook continue to advance PCML behind this Github repository so stay tuned for some very exciting updates very soon. 

PCML in action
--------------
PCML has been used in several graduate and undergraduate courses taught at Kent State University and the University of Minnesota including 'Introduction to CyberGIS' and 'Advanced Geocomputing.' It was featured in a workshop (see the 30 minutes PCML demo below).  We also have plans in place to teach it in several courses this Spring 2017 semester at different institutions.

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
For detailed instructions please see https://github.com/HPCGISLab/pcml/blob/master/demo_30min_install.pdf

#### 0. Download and install Ananconda (Optional)

   The Anaconda Python distribution package is available at http://continuum.io/downloads.
   It provides Python and several Python packages, but can be skipped if you have Python already.

#### 1. Download PCML
   Navigate to this page (i.e., the PCML page on GitHub at https://github.com/hpcgislab/pcml)

   Scroll up and click the "Download ZIP" button on the right side.

#### 2.	Unzip PCML

   A file named “PCML-master.zip” will be downloaded to your computer. Unzip it.

#### 3. Open ‘demo_30min.py’ 

   Open 'demo_30min.py' by clicking File->Open in the Spyder editor (included in the Ananconda package)

#### 4. Create your own LocalSum function in under 30 minutes!

   Open 'demo_30min_walkthrough.pdf' and follow the 30 minute tutorial.

#### 5. Run PCML and start processing spatial data in parallel

   Success!
   
#### 6. Run into a problem or a bug? Or have a great idea?

   Let us know!  We are trying to make PCML as easy-to-use and as useful as possible.
   We live to collaborate and work on interesting problems (and we will debug your PCML code for you).
   
   Submit a 'New issue' on github (https://github.com/HPCGISLab/pcml/issues) or email Eric Shook (eshook@gmail.com).
