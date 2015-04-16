#!/usr/bin/python
"""
Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); 
"""

# Import the PCML package to start developing in the parallel cartographic modeling language (PCML)
from pcml import *

# Import the operation system (os) package path
# so that it works with MS Windows, Mac OS X, or Linux 
import os.path as path

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 
# 1. In this section, you will create your own myLocalSum operation
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #







# This is the 'main' section that will execute our code
if __name__ == '__main__':

    print("In 30 minutes or less you can begin developing in the parallel cartographic modeling language.")
    print("You will need to create and call your own operation, and if you are bold you can compare your")
    print("results with our results to make sure your operation is correct.")
    print("You can do this in 5 easy steps and we walk you through every step of the way.")

    # Default data directory (where our data -- dataa.asc and datab.asc -- are located)
    datadir="data"

    # Read 2 layers that we will add together using PCML 
    layera=ReadASCIIGrid(path.join(datadir,"dataa.asc"))
    layerb=ReadASCIIGrid(path.join(datadir,"datab.asc"))
    
    # Print out the layers so that we can see the data 
    print("layera",layera) # This prints the layer information and bounding box
    layera.print_data()    # This prints the data 

    print("layerb",layerb)
    layerb.print_data()

    # Add the two layers together using PCML-style map algebra
    sumresult=layera+layerb

    # Print the sum results
    print("sumresult",sumresult) # This prints the layer information and bounding box
    sumresult.print_data()       # This prints the data

    # Write the sum results to a file named 'sumresult.asc' in the data directory
    WriteASCIIGrid(path.join(datadir,"sumresult.asc"), sumresult)

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # 
    # 2. In this section, call your myLocalSum operation and save the results to a 'myresult' layer
    #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #




    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # 
    # 3. In this section, print your 'myresult' layer (layer information and the data)
    #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    #print("myresult",myresult)
    #myresult.print_data()

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # 
    # 4. Write your 'myresult' layer to a file named 'myresult.asc' by uncommenting the line of code
    #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    #WriteASCIIGrid(path.join(datadir,"myresult.asc"), myresult)

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # 
    # 5. If you are bold, then uncomment the next three lines of code to see if the difference is zero! 
    #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    #differencelayer=myresult-sumresult
    #print("\n\nDifference (should be all zeroes)",differencelayer)
    #differencelayer.print_data()

