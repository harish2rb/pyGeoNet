# pyGeoNet_readGeotiff
#import sys
#import os

from osgeo import gdal
#from string import *
import numpy as np
from time import clock
import pygeonet_defaults as defaults
import pygeonet_prepare as Parameters
from math import modf, floor
#from scipy.stats.mstats import mquantiles


def read_dem_from_geotiff(demFileName,demFilePath):
    # Open the GeoTIFF format DEM
    fullFilePath = demFilePath + demFileName
    #fullFilePath = "G:\\HarishLaptop_Backup\\TI102782W0E\\PythonScripts\\pyGeoNet1.0\\data\\skunk.tif"
    print fullFilePath
    ary = []
    ds = gdal.Open(fullFilePath, gdal.GA_ReadOnly)
    geotransform = ds.GetGeoTransform()
    '''
    print 'Driver: ', ds.GetDriver().ShortName,'/', \
          ds.GetDriver().LongName
    print 'Size is ',ds.RasterXSize,'x',ds.RasterYSize, \
          'x',ds.RasterCount
    print 'Projection is ',ds.GetProjection()
    
    
    if not geotransform is None:
        print 'Origin = (',geotransform[0], ',',geotransform[3],')'
        print 'Pixel Size = (',geotransform[1], ',',geotransform[5],')'
    '''
    ary = ds.GetRasterBand(1).ReadAsArray()
    #Parameters.geospatialReferenceArray
    #Parameters.geoReferencingMatrix
    #Parameters.geoBoundingBox
    Parameters.demPixelScale = geotransform[1]
    Parameters.xLowerLeftCoord = geotransform[0]
    Parameters.yLowerLeftCoord = geotransform[3]
    return ary


def quantile(x, q,  qtype = 7, issorted = False):
	"""
	Args:
	   x - input data
	   q - quantile
	   qtype - algorithm
	   issorted- True if x already sorted.
 
	Compute quantiles from input array x given q.For median,
	specify q=0.5.
 
	References:
	   http://reference.wolfram.com/mathematica/ref/Quantile.html
	   http://wiki.r-project.org/rwiki/doku.php?id=rdoc:stats:quantile
 
	Author:
	Ernesto P.Adorio Ph.D.
	UP Extension Program in Pampanga, Clark Field.
	"""
	if not issorted:
		y = sorted(x)
	else:
		y = x
	if not (1 <= qtype <= 9):
	   return None  # error!
 
	# Parameters for the Hyndman and Fan algorithm
	abcd = [(0,   0, 1, 0), # inverse empirical distrib.function., R type 1
			(0.5, 0, 1, 0), # similar to type 1, averaged, R type 2
			(0.5, 0, 0, 0), # nearest order statistic,(SAS) R type 3
 
			(0,   0, 0, 1), # California linear interpolation, R type 4
			(0.5, 0, 0, 1), # hydrologists method, R type 5
			(0,   1, 0, 1), # mean-based estimate(Weibull method), (SPSS,Minitab), type 6
			(1,  -1, 0, 1), # mode-based method,(S, S-Plus), R type 7
			(1.0/3, 1.0/3, 0, 1), # median-unbiased ,  R type 8
			(3/8.0, 0.25, 0, 1)   # normal-unbiased, R type 9.
		   ]
 
	a, b, c, d = abcd[qtype-1]
	n = len(x)
	g, j = modf( a + (n+b) * q -1)
	if j < 0:
		return y[0]
	elif j >= n:
		return y[n-1]   # oct. 8, 2010 y[n]???!! uncaught  off by 1 error!!!
 
	j = int(floor(j))
	if g ==  0:
	   return y[j]
	else:
	   return y[j] + (y[j+1]- y[j])* (c + d * g)

def main():    
    #demFileName = "skunk.tif"
    #demFilePath = "G:\\HarishLaptop_Backup\\TI102782W0E\\PythonScripts\\pyGeoNet1.0\\data\\"
    print "Reading input file path :",Parameters.demDataFilePath
    print "Reading input file :",Parameters.demFileName
    rawDemArray = read_dem_from_geotiff(Parameters.demFileName,Parameters.demDataFilePath)
    
    nanDemArray=rawDemArray
    nanDemArray[nanDemArray < defaults.demNanFlag]= np.NAN
    Parameters.minDemValue= np.min(nanDemArray[:])
    Parameters.maxDemValue= np.max(nanDemArray[:])
    # Area of analysis
    Parameters.xDemSize=np.size(rawDemArray,0)
    Parameters.yDemSize=np.size(rawDemArray,1)
    # Calculate pixel length scale and assume square
    Parameters.maxLowerLeftCoord = np.max([Parameters.xDemSize, Parameters.yDemSize])
    print 'DTM size: ',Parameters.xDemSize, 'x' ,Parameters.yDemSize
    #-----------------------------------------------------------------------------

    # Compute slope magnitude for raw and filtered DEMs
    print 'Computing slope of raw DTM'
    slopeMagnitudeDemArray = np.gradient(nanDemArray,Parameters.demPixelScale)
    print slopeMagnitudeDemArray

    # Computation of the threshold lambda used in Perona-Malik nonlinear
    # filtering. The value of lambda (=edgeThresholdValue) is given by the 90th
    # quantile of the absolute value of the gradient.
    print'Computing lambda = q-q-based nonlinear filtering threshold'
    mult = Parameters.xDemSize * Parameters.yDemSize
    print np.size(slopeMagnitudeDemArray,0)
    edgeThresholdValue = quantile(np.reshape(slopeMagnitudeDemArray,mult),defaults.demSmoothingQuantile)    
    print edgeThresholdValue

  


if __name__ == '__main__':
    t0 = clock()
    main()
    t1 = clock()
    print "time taken to complete the script is::",t1-t0," seconds"
    print "script complete"



