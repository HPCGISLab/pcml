"""Copyright (c) 2014 High-Performance Computing and GIS (HPCGIS) Laboratory. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
Authors and contributors: Eric Shook (eshook@kent.edu); Jayakrishnan Ajayakumar (jajayaku@kent.edu)  
"""
from ..core.Operation import *
from ..core.Scheduler import *
from ..util.OperationBuilder import *
import pcml.core.PCMLConfig as PCMLConfig
import numpy as np
import types
import math

try:
   PCMLConfig.numbaenabled = 1
   from numba import jit
except ImportError as e:
   PCMLConfig.numbaenabled = 0

try:
   PCMLConfig.scipyenabled = 1
   from scipy import ndimage
except ImportError as e:
   PCMLConfig.scipyenabled = 0


#from scipy import stats

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

if PCMLConfig.numbaenabled == 1:

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
    parameterdata=np.array([subdomains[1].cellsize,-(subdomains[1].cellsize),sin_altRadians,azRadians,cos_altRadians_mul_z_scale_factor,square_z_scale_factor])
    outarrdim=np.array([subdomains[0].nrows,subdomains[0].ncols,subdomains[0].r,subdomains[0].c])
    datarrdim=np.array([subdomains[1].nrows,subdomains[1].ncols,subdomains[1].r,subdomains[1].c])
    numbashadecalculator_new(outarray,outarrdim,dataarray,datarrdim,parameterdata,subdomains[1].nodata_value)

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

# Makes a footprint array for ndimage filters (can be used by Buffer for example)
def _makefp(dist,cellsize):
    # Calculate the array size for footprint
    fpsize=int(math.ceil(dist/cellsize)*2)+1
    # Calculate the "raster distance" in number of cells from center
    rdist=math.ceil(dist/cellsize)
    # Create footprint array
    fp=np.zeros((fpsize,fpsize))
    for i in xrange(fpsize):
        for j in xrange(fpsize):
            xd=(rdist-i)*cellsize
            yd=(rdist-j)*cellsize
            d=math.sqrt(xd*xd + yd*yd) # Calculate distance
            if(d<=dist):               # if less than dist, 1 else 0
                fp[i][j]=1
            else:
                fp[i][j]=0
    return fp

# Function for nearby classify used in scipy.ndimage.generic_filter
def _nearbyclassifyfunct(values,classifyval):
    if classifyval in values: # If we find classifyval in the search space
        return 1              # Return true (1)
    else:
        return 0

# Function for buffer to be used in scipy.ndimage.generic_filter
def _bufferfunct(values,classifyval):
    if classifyval in values: # If we find classifyval in the search space
        return classifyval    # Return it (which acts as a buffer)
    else:
        return values[len(values)/2]  # Otherwise return the original value (in the middle of values array)

@executor
@focaloperation
def FocalBuffer(self,subdomains):
    if PCMLConfig.scipyenabled == 0:
        PCMLNotSupported("SciPy module required for getneighbors()")
 
    # Get the array from the output subdomain as well as the output subdomain
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()
    # Get the input subdomain
    insubdomain = subdomains[1]
    inarr = insubdomain.get_nparray()

    classifyval = self.kwargs.get('classtype',0)

    # Create the footprint array for ndimage filter
    fp = _makefp(self.buffersize,outsubdomain.cellsize)

    arr=ndimage.generic_filter(inarr,function=_bufferfunct,footprint=fp,extra_arguments=(classifyval,))

    # Copy values to outarr (outsubdomain)
    roffset=outsubdomain.r-insubdomain.r
    coffset=outsubdomain.c-insubdomain.c
    outarr[:,:]=arr[roffset:outsubdomain.nrows+roffset,coffset:outsubdomain.ncols+coffset]
    

@executor
@focaloperation
def BufferedClassify(self,subdomains):
    # Get the array from the output subdomain as well as the output subdomain
    outsubdomain = subdomains[0]
    outarr = outsubdomain.get_nparray()
    # Get the input subdomain
    insubdomain = subdomains[1]
    inarr = insubdomain.get_nparray()

    classifyval = self.kwargs.get('classtype',0)

    # Create the footprint array for ndimage filter
    fp = _makefp(self.buffersize,outsubdomain.cellsize)

    arr=ndimage.generic_filter(inarr,function=_nearbyclassifyfunct,footprint=fp,extra_arguments=(classifyval,))

    # Only get the inside portion of the insubdomain array (this complex thing trims out ghost zones)
    roffset=outsubdomain.r-insubdomain.r
    coffset=outsubdomain.c-insubdomain.c
    outarr[:,:]=arr[roffset:outsubdomain.nrows+roffset,coffset:outsubdomain.ncols+coffset]
    #outarr[:,:]=arr[outsubdomain.r-insubdomain.r:outsubdomain.nrows,outsubdomain.c-insubdomain.c:outsubdomain.ncols]
    


