"""Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Jayakrishnan Ajayakumar (jajayaku@kent.edu)  
"""
from ..core.Operation import *
from ..core.Scheduler import *
from ..util.OperationBuilder import *
import numpy as np
import types
import math
from numba import jit
#from scipy import stats

# For LeastCostDistance
def sortlist(unsorted_list):
    # Added -x to reverse the sort
    sort_on = lambda pos: lambda x: -x[pos]
    return sorted(unsorted_list, key=sort_on(2))

# For LeastCostDistance
def boundscheck(r,c,arr):
    if r < 0 or c < 0:
        return False
    if r >= len(arr) or c >= len(arr[0]):
        return False
    return True


@executor
@focaloperation
def LeastCostDistance(self, subdomains):

    # Get the array from the output subdomain as well as the output subdomain
    # This will be the source (which will be copied from input)
    # and used in the computation
    outsubdomain = subdomains[0]
    source = outsubdomain.get_nparray() 

    # Get the cost subdomain
    costsubdomain = subdomains[1]
    cost = costsubdomain.get_nparray() 

    # Get the source subdomain
    sourcesubdomain = subdomains[2]
    sourceinput = sourcesubdomain.get_nparray() 

    '''
    print "out nrows,ncols",outsubdomain.nrows,outsubdomain.ncols
    print "cst nrows,ncols",costsubdomain.nrows,costsubdomain.ncols
    print "src nrows,ncols",sourcesubdomain.nrows,sourcesubdomain.ncols

    print "out r,c",outsubdomain.r,outsubdomain.c
    print "cst r,c",costsubdomain.r,costsubdomain.c
    print "src r,c",sourcesubdomain.r,sourcesubdomain.c
    '''


    # Track the number of locations that were processed in this iteration
    locations_processed = 0

    # Start with an empty list of new locations with values 
    new_list_val=[]

    # If this is the first time applying least-cost-path to this subdomain
    if outsubdomain.iteration_number == 0:
        # Copy input source to output source
        outsubdomain.copy_boundingbox_data(sourcesubdomain)

        # Find the start (i.e. 0 value) locations
        new_list_inv=np.where(source==0)
        #print "nli=",new_list_inv

        # Extract the two arrays of X and Y dimension and
        # zip (invert) them into coordinate sets
        new_list=zip(new_list_inv[0], new_list_inv[1])
        #print "new_list=",new_list

        for loc in new_list:
            #new_list_val.append([loc[0],loc[1],cost[loc]])
            new_list_val.append([loc[0], loc[1], 0]) # Initial cost is zero

    # If this is the n'th time applying last-cost-path
    # then we only look at the edges.
    # Notice this is not terribly efficient, but it will be correct
    else:
        # Get the ghost zone edges for the subdomain (N, S, E, W)
        # Use xor hash on edges to detect change for optimization
        # If there were changes to last time

        # Add the edge locations to a list of new locations with values
        outlayer = outsubdomain.layer
        outraster = outlayer.get_nparray()

        # ASSUMPTION: Row decomposition only, need to add row check for column/grid        
        # Check the upper and lower ghost zone 
        for r in [outsubdomain.r-1,outsubdomain.r+outsubdomain.h]:
            for c in xrange(outsubdomain.ncols):
                if not boundscheck(r,c,outraster): # Then this is not a valid location in the layer
                    continue                        # So skip it
                if outraster[r,c] == -9999:         # Then this is not set yet
                    continue                        # So skip it
                # Calculate the offset to switch location from global to local for grabbing the cost 
                #ro = r - (outsubdomain.r - sourcesubdomain.r)
                #co = c - (outsubdomain.c - sourcesubdomain.c)
                #loccost = cost[ro,co]

                #point = [r-outsubdomain.r,c-outsubdomain.c,loccost]

                point = [r-outsubdomain.r,c-outsubdomain.c,outraster[r,c]]

                #print "ghost zone possible point",point," global",point[0]+outsubdomain.r,point[1]+outsubdomain.c

                # Now calculate the offset for the outsubdomain (I know it is crazy)
                new_list_val.append(point)


        #print "ROW",outsubdomain.r,"spl",new_list_val
        #for p in new_list_val:
        #    print "   offset potential",p[0]+outsubdomain.r,p[1]+outsubdomain.c
        #return
        #return

    # Sort the list
    sorted_potential_list=sortlist(new_list_val)

    #print "ROW",outsubdomain.r,"spl",sorted_potential_list
    #for p in sorted_potential_list:
    #    print "   offset potential",p[0]+outsubdomain.r,p[1]+outsubdomain.c

    # While potential locations are available, process them
    while sorted_potential_list:
        # Take the shortest location to process
        location = sorted_potential_list.pop()
        #print "PROCESSING location=",location,"abs",location[0]+outsubdomain.r,location[1]+outsubdomain.c

        # Get neighboring locations (3x3 window to location)
        for r in xrange(location[0]-1,location[0]+2):
            for c in xrange(location[1]-1,location[1]+2):
                # Calculate offset for sourceinput/cost to source
                # These must be used whenever referencing sourceinput/cost to function properly
                roff = - (outsubdomain.r - sourcesubdomain.r)
                coff = - (outsubdomain.c - sourcesubdomain.c)

                if r == location[0] and c == location[1]: # Skip self
                    continue
                # Boundary check gets tricky with ghost zones, because the checks go through with r,c offsets
                #if (not boundscheck(r,c,source)) and (not boundscheck(r-roff,c-coff,source)):
                if not boundscheck(r,c,source):
                     continue
                pt = [r,c,source[r,c]]

                if cost[r-roff,c-coff] == 9999: # Don't write skipped cells
                    #print "\t\t\t\tfound 9999 location = ",pt[0]+outsubdomain.r,pt[1]+outsubdomain.c
                    source[r,c]=9999 # Copy the skipped cell to source and skip processing the point
                    continue

                # Calculate cost
                costval = float(cost[location[0]-roff,location[1]-coff] + cost[r-roff,c-coff]) / 2
                if location[0] != r and location[1] != c: # Then we are diagonal
                     costval*=1.41421356237 # Multiply by sqrt(2)
                costval = costval + location[2]

                # If the new cost val will be the first or smaller than the current costval
                if pt[2] == -9999 or pt[2] > costval:
                    locations_processed += 1

                    # Save cost in source and point
                    source[r,c] = costval
                    pt[2] = costval
                    # Add this point as a potential for new neighbors
                    sorted_potential_list.append(pt)

            #print "source=\n",np.array_str(source,precision=2,suppress_small=True)

        # Resort the list
        sorted_potential_list = sortlist(sorted_potential_list)

    print "locations processed", locations_processed

    print "r,c=",costsubdomain.r,",",costsubdomain.c,"cost=\n", cost
    print "source=\n", np.array_str(source, precision=1, suppress_small=True)

    if locations_processed == 0: # If no updates
        # If we are already over-processing, mark as more negative
        # otherwise start with smallest negative (-1)
        if outsubdomain.iteration_number < 0:
            outsubdomain.iteration_number -= 1
            print "It was negative so going in the hole"
        else:
            outsubdomain.iteration_number = -1
    else: # If updates, mark how many
        outsubdomain.iteration_number = locations_processed

    print "iteration_number",outsubdomain.iteration_number

@executor
@focaloperation
def FocalMean_np_exec(self, subdomains):
    # Get the array from the output subdomain as well as the output subdomain
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray() 
    # Get the input subdomain
    insubdomain = subdomains[1]
    inarr = insubdomain.get_nparray() 

    # Get the buffersize for the focal operation
    buffersize=self.buffersize

    locind={}

    for roworig in xrange(outsubdomain.nrows):
        for colorig in xrange(outsubdomain.ncols):
            # Calculate the row and column positions (global for the study area)
            # which is why we add outsubdomain.{r,c}
            locind['r']=roworig+outsubdomain.r
            locind['c']=colorig+outsubdomain.c

            # The following 2 lines of code are the easy-to-read version of the code below
            # Essentially we get the location for the input array and then average it
            # Get the array from the location
            #arr=insubdomain.bufferedlocgetarr(locind,buffersize)
            # Calculate the mean and set it
            #outarr[roworig][colorig]=np.average(arr)

            # Add the offsets for the subdomains (difference etween out/input subdomains) 
            row=locind['r']-insubdomain.r
            col=locind['c']-insubdomain.c
            # Calculate the starting row and column accounting for buffersize and the edge of layer
            r=max(0,row-buffersize)
            c=max(0,col-buffersize)

            # Calculate the h accounting for buffersize and edge of layer
            h=buffersize+(row-r)+1
            if (r+h > insubdomain.nrows):
                h=insubdomain.nrows-r

            # Same for width 
            w=buffersize+(col-c)+1
            if (c+w > insubdomain.ncols):
                w=insubdomain.ncols-c

            # Take a slice of the array (inarr) based on those calculations, average the results
            outarr[roworig][colorig]=np.average(inarr[r:r + h, c:c + w])

#Helper function to claculate shade using Numba
@jit(nopython=True)
def numbashadecalculator_new(outarray,outarrdim,inputarray,inputarrdim,parameterdata,nodata_value):
    buffersize=1
    for i in xrange(outarrdim[0]):
        for j in xrange(outarrdim[1]):
            rout,cout=i+outarrdim[2],j+outarrdim[3]
            rin,cin=rout-inputarrdim[2],cout-inputarrdim[3]
            r,c=rin-buffersize,cin-buffersize
            if r<0:
                r=0
            if c<0:
                c=0
            h=buffersize+(rin-r)+1
            if r+h > inputarrdim[0]:
                h=inputarrdim[0]-r
            w=buffersize+(cin-c)+1
            if c+w > inputarrdim[1]:
                w=inputarrdim[1]-c
            arr=inputarray[r:r+h,c:c+w]
            arraysize=arr.size
            if(arraysize!= 9):
                outarray[i][j]=nodata_value
                continue
            containsnodata=False
            for ii in xrange(3):
                for jj in xrange(3):
                    if arr[ii][jj]==nodata_value:
                        containsnodata=True
                        break
                if (containsnodata):
                    break
            if (containsnodata):
                outarray[i][j]=nodata_value
                continue
            dzdx=((arr[2][2]+2*arr[1][2]+arr[0][2]) - (arr[2][0]+2*arr[1][0]+arr[0][0]))/ parameterdata[0]
            dzdy=((arr[2][0]+2*arr[2][1]+arr[2][2]) - (arr[0][0]+2*arr[0][1]+arr[0][2])) / parameterdata[0]
            xx_plus_yy = (dzdx*dzdx) + (dzdy*dzdy)
            aspect=np.arctan2(dzdy,-dzdx)
            shade=(parameterdata[2] - parameterdata[4] * np.sqrt(xx_plus_yy) * np.sin(aspect - parameterdata[3]))/np.sqrt(1+parameterdata[5] * xx_plus_yy)
            if (shade <= 0):
                shade=1.0
            else:
                shade=1.0 + (254.0 * shade)
            outarray[i][j]=np.ceil(shade)


@executor
@focaloperation
def HillShade_Exec_Numba(self,subdomains): # Experimental
    start=timeit.default_timer()
    dataarray=GetDataForBoundingBox(subdomains[1])
    altitude=45
    #altitude=60 # Override standard
    azimuth=315
    rtod=3.1415926535897932384626433832795/180.0 # ( pi / 180.0 )
    #zenith_rad=(90-altitude)*rtod
    #azimuth_rad=(360.0-azimuth+90)*rtod
    outarray=np.zeros((subdomains[0].nrows,subdomains[0].ncols))
    #constval=8*subdomains[1].cellsize
    #constval=np.array([subdomains[1].nsres,subdomains[1].ewres])
    #new stuff
=======
            outarray[i][j]=np.ceil(shade)   

#Hillshade calculation using numba
@executor
@focaloperation
def HillShade_Exec_Numba(self,subdomains): # Experimental
    dataarray=subdomains[1].get_nparray()
    altitude=45
    azimuth=315
    rtod=3.1415926535897932384626433832795/180.0 # ( pi / 180.0 )
    outarray=subdomains[0].get_nparray()

    sin_altRadians = np.sin(altitude*rtod)
    azRadians=azimuth*rtod
    z_scale_factor = 1/8.0
    cos_altRadians_mul_z_scale_factor=np.cos(altitude * rtod) * z_scale_factor
    square_z_scale_factor = z_scale_factor * z_scale_factor
    parameterdata=np.array([subdomains[1].nsres,subdomains[1].ewres,sin_altRadians,azRadians,cos_altRadians_mul_z_scale_factor,square_z_scale_factor])
    outarrdim=np.array([subdomains[0].nrows,subdomains[0].ncols,subdomains[0].r,subdomains[0].c])
    datarrdim=np.array([subdomains[1].nrows,subdomains[1].ncols,subdomains[1].r,subdomains[1].c])
    numbashadecalculator_new(outarray,outarrdim,dataarray,datarrdim,parameterdata,subdomains[1].nodata_value)
    end=timeit.default_timer()
    print "time taken on numba based hillshade calculation alone :",end-start
    AddBoundingBoxData(subdomains[0],outarray,self.lock)


#Helper function to calculate contour lines using numba
def countour_line_calc(outarray,outarrdim,inputarray,inputarrdim,buffersize,nodata_value):
    cntrline=5.0
    for i in xrange(outarrdim[0]):
        for j in xrange(outarrdim[1]):
            rout,cout=i+outarrdim[2],j+outarrdim[3]
            rin,cin=rout-inputarrdim[2],cout-inputarrdim[3]
            r,c=rin-buffersize,cin-buffersize
            if r<0:
                r=0
            if c<0:
                c=0
            h=buffersize+(rin-r)+1
            if r+h > inputarrdim[0]:
                h=inputarrdim[0]-r
            w=buffersize+(cin-c)+1
            if c+w > inputarrdim[1]:
                w=inputarrdim[1]-c
            arr=inputarray[r:r+h,c:c+w]
            arraysize=arr.size
            if(arraysize!= 9):
                outarray[i][j]=nodata_value
                continue
            containsnodata=False
            for ii in xrange(3):
                for jj in xrange(3):
                    if arr[ii][jj]==nodata_value:
                        containsnodata=True
                        break
                if (containsnodata):
                    break
            if (containsnodata):
                outarray[i][j]=nodata_value
                continue
            q=np.floor(arr[1][1]/cntrline)
            r=arr[1][1]%cntrline
            cutoff=(q+1)*cntrline
            if(arr[2][1]>=cutoff or arr[1][2]>=cutoff or arr[0][1]>=cutoff or arr[1][0]>=cutoff):
                outarray[i][j]=1
                continue
            sum=0
            for tt in xrange(3):
                for ss in xrange(3):
                    sum=sum+arr[tt][ss]
            if(np.floor(sum/arraysize)%10 != 0):
                outarray[i][j]=1
            else:
                outarray[i][j]=0

#Function to calculate contour lines using Numba
@executor
@focaloperation
def Contour_Lines_Numba(self,subdomains): # Experimental
    dataarray=subdomains[1].get_nparray()
    outarray=subdomains[0].get_nparray()
    outarrdim=np.array([subdomains[0].nrows,subdomains[0].ncols,subdomains[0].r,subdomains[0].c])
    datarrdim=np.array([subdomains[1].nrows,subdomains[1].ncols,subdomains[1].r,subdomains[1].c])
    countour_line_calc(outarray,outarrdim,dataarray,datarrdim,self.buffersize,subdomains[1].nodata_value)

