# pyGeoNet_readGeotiff
import sys
import os

from osgeo import gdal,osr,ogr
import statsmodels.api as sm
import numpy as np
from time import clock
import pygeonet_defaults as defaults
import pygeonet_prepare as Parameters
from math import modf, floor
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.stats.mstats import mquantiles
import bottleneck as bn
import scipy.signal as conv2
from scipy import stats
import skfmm
from scipy import ndimage
import numpy.ma as npma

# ----------- GRASS GIS SETUP ------------------
#setting up the environment for grass gis access
"""
# The below grass script will work assuming you have installed
# Grass GIS 7 on your machine and that the required environment
# variables are set on windows /linux as required

"""
sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))
import grass.script as g
import grass.script.setup as gsetup

# -------------------------- FUNCTIONS ----------------------------
def read_dem_from_geotiff(demFileName,demFilePath):
    # Open the GeoTIFF format DEM
    fullFilePath = demFilePath + demFileName
    print fullFilePath
    ary = []
    gdal.UseExceptions()
    ds = gdal.Open(fullFilePath, gdal.GA_ReadOnly)
    Parameters.driver = ds.GetDriver()
    geotransform = ds.GetGeoTransform()
    Parameters.geotransform = geotransform
    #'''
    print 'Driver: ', ds.GetDriver().ShortName,'/', \
          ds.GetDriver().LongName
    print 'Size is ',ds.RasterXSize,'x',ds.RasterYSize, \
          'x',ds.RasterCount
    print 'Projection is ',ds.GetProjection()
    print 'geotransform', geotransform,type(geotransform[0])
    
    if not geotransform is None:
        print 'Origin = (',geotransform[0], ',',geotransform[3],')'
        print 'Pixel Size = (',geotransform[1], ',',geotransform[5],')'
    #'''
    ary = ds.GetRasterBand(1).ReadAsArray()
    #Parameters.geospatialReferenceArray
    #Parameters.geoReferencingMatrix
    #Parameters.geoBoundingBox
    Parameters.demPixelScale = float(geotransform[1])
    Parameters.xLowerLeftCoord = float(geotransform[0])
    Parameters.yLowerLeftCoord = float(geotransform[3])
    Parameters.inputwktInfo = ds.GetProjection()
    
    return ary

def quantileasmatlab(a, prob):
    """
    Estimates the prob'th quantile of the values in a data array.

    Uses the algorithm of matlab's quantile(), namely:
        - Remove any nan values
        - Take the sorted data as the (.5/n), (1.5/n), ..., (1-.5/n) quantiles.
        - Use linear interpolation for values between (.5/n) and (1 - .5/n).
        - Use the minimum or maximum for quantiles outside that range.

    See also: scipy.stats.mstats.mquantiles
    """
    a = np.asanyarray(a)
    a = a[np.logical_not(np.isnan(a))].ravel()
    n = a.size

    if prob >= 1 - .5/n:
        return a.max()
    elif prob <= .5 / n:
        return a.min()

    # find the two bounds we're interpreting between:
    # that is, find i such that (i+.5) / n <= prob <= (i+1.5)/n
    t = n * prob - .5
    i = np.floor(t)

    # partial sort so that the ith element is at position i, with bigger ones
    # to the right and smaller to the left
    a = bn.partsort(a, i)

    if i == t: # did we luck out and get an integer index?
        return a[i]
    else:
        # we'll linearly interpolate between this and the next index
        smaller = a[i]
        larger = a[i+1:].min()
        if np.isinf(smaller):
            return smaller # avoid inf - inf
        return smaller + (larger - smaller) * (t - i)


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

def simple_gaussian_smoothing(inputDemArray,kernelWidth,diffusionSigmaSquared):
    """
    smoothing input array with gaussian
    """
    [Ny,Nx]=inputDemArray.shape;
    #print Ny,Nx
    halfKernelWidth=(kernelWidth-1)/2;
    #print "halfKernelWidth",halfKernelWidth
    # Make a ramp array with 5 rows each containing [-2, -1, 0, 1, 2]
    x = np.linspace(-halfKernelWidth, halfKernelWidth, kernelWidth)
    #x= np.ones((kernelWidth,1))* range(-halfKernelWidth,halfKernelWidth+1)
    #y=x.T
    y = x
    xv,yv = np.meshgrid(x,y)
    gaussianFilter = np.exp(-(xv**2+yv**2)/(2*diffusionSigmaSquared))  # 2D Gaussian
    #print "gaussianFilter", gaussianFilter
    gaussianFilter=gaussianFilter/np.sum(gaussianFilter[:]) # Normalize
    #print inputDemArray[:,0:halfKernelWidth]
    xL= np.nanmean(inputDemArray[:,0:halfKernelWidth],axis=1)
    #xL = np.matrix(xL)
    xR= np.nanmean(inputDemArray[:,Nx-halfKernelWidth:Nx],axis =1)
    #xR = np.matrix(xR)
    #print "xR",xR.shape
    part1T = np.vstack((xL,xL))
    part1 = part1T.T
    #part1 = xL * np.matrix(np.ones((1,halfKernelWidth)))
    part2T = np.vstack((xR,xR))
    part2 = part2T.T
    #part2 = xR * np.matrix(np.ones((1,halfKernelWidth)))
    #print part1.shape, part2.shape
    eI = np.hstack((part1,inputDemArray,part2))
    #eI= np.matrix(np.concatenate((part1,inputDemArray,part2),1))
    #print 'eI',eI.shape
    xU= np.nanmean(eI[0:halfKernelWidth,:],axis=0)
    #xU = np.matrix(xU)
    #print "xU",xU.shape
    xD= np.nanmean(eI[Ny-halfKernelWidth:Ny,:],axis =0)
    #xD = np.matrix(xD)
    #print "xD",xD.shape
    part3 = np.vstack((xU,xU))
    part4 = np.vstack((xD,xD))
    #print part3.shape, part4.shape
    #part3 = np.matrix(np.ones((halfKernelWidth,1)))*xU
    #part4 = np.matrix(np.ones((halfKernelWidth,1)))*xD
    # Generate the expanded DTM array, 4 pixels wider in both x,y directions
    eI = np.vstack((part3,eI,part4))
    #eI= np.matrix(np.concatenate((part3, eI, part4),0))
    #print 'eI',eI.shape
    # The 'valid' option forces the 2d convolution to clip 2 pixels off the edges
    # NaNs spread from one pixel to a 5x5 set centered on the NaN
    #smoothedDemArray=conv2.convolve2d(eI,gaussianFilter,'valid'); # original
    fillvalue = np.nanmean(inputDemArray[:])
    smoothedDemArray=conv2.convolve2d(eI,gaussianFilter,'valid')
    return smoothedDemArray


def geonet_diffusion(demArray, diffusionMethod, nFilterIterations,\
    edgeThreshold, diffusionTimeIncrement, diffusionSigmaSquared, pixelSize):
    """
	References:
	   Based on diffusion() by Guy Gilboa
	   Code imported from GeoNet2.1 by Harish Sangireddy June 2014
    """
    print 'Performing Perona Malik'
    print diffusionMethod, nFilterIterations,edgeThreshold, \
    diffusionTimeIncrement, diffusionSigmaSquared, pixelSize
    # DTM dimensions
    [Ny,Nx]=demArray.shape;
    print Ny,Nx
    for i in range(0,nFilterIterations):
        print "iteration",i
        # Gaussian filter the DTM using a 5x5 kernel (Catte et al)
        if diffusionSigmaSquared>0:
            originalDemArray = demArray   # Save original DTM array
            demArrayout = simple_gaussian_smoothing(demArray,5,diffusionSigmaSquared)
        del demArray
        demArray = demArrayout
        #print 'demArray after gaussian smoothing',demArray.shape
        # Now calculate gradient in all directions (N,S,E,W) by simple differencing
        # - with repeat padding in each direction.
        # This step will propagate NaNs one pixel inward in each dirn.
        demArrayMatrix = np.matrix(demArray)
        In = (np.concatenate((demArrayMatrix[0,:],demArrayMatrix[0:Ny-1,:]),0)\
              - demArrayMatrix)/pixelSize
        #print In.shape
        #stop
        Is = (np.concatenate((demArrayMatrix[1:Ny,:],demArrayMatrix[Ny-1,:]),0)\
              - demArrayMatrix)/pixelSize
        #print Isgr.shape, Is.shape
        #print np.nanmean(Isgr - Is)
        #stop
        Ie = (np.concatenate((demArrayMatrix[:,1:Nx],demArrayMatrix[:,Nx-1]),1)\
              - demArrayMatrix)/pixelSize
        #print Ie.shape
        
        Iw = (np.concatenate((demArrayMatrix[:,0],demArrayMatrix[:,0:Nx-1]),1)\
              - demArrayMatrix)/pixelSize
        #print Iw.shape
        In[np.isnan(In)] = 0
        Is[np.isnan(Is)] = 0
        Ie[np.isnan(Ie)] = 0
        Iw[np.isnan(Iw)] = 0
        
        # Calculate diffusion coefficients in all dirns according to diffusionMethod
        if diffusionMethod =='linear':
            Cn=edgeThreshold
            Cs=edgeThreshold
            Ce=edgeThreshold
            Cw=edgeThreshold
        elif diffusionMethod == 'PeronaMalik1':
            Cn=np.exp(-(np.abs(In)/edgeThreshold)**2)
            Cs=np.exp(-(np.abs(Is)/edgeThreshold)**2)
            Ce=np.exp(-(np.abs(Ie)/edgeThreshold)**2)
            Cw=np.exp(-(np.abs(Iw)/edgeThreshold)**2)
        elif diffusionMethod == 'PeronaMalik2':
            Cn=1/(1+ np.array(np.abs(In)/edgeThreshold)**2)
            Cs=1/(1+ np.array(np.abs(Is)/edgeThreshold)**2)
            Ce=1/(1+ np.array(np.abs(Ie)/edgeThreshold)**2)
            Cw=1/(1+ np.array(np.abs(Iw)/edgeThreshold)**2)
        else:
            print 'Unknown smoothing method', diffusionMethod

        if diffusionSigmaSquared>0:
            #Calculate real gradients (not smoothed) - with repeat padding in each
            # direction.  This step will propagate NaNs one pixel inward in each dirn.
            originalDemArrayMatrix = np.matrix(originalDemArray)
            In = (np.concatenate((originalDemArrayMatrix[0,:],originalDemArrayMatrix[0:Ny-1,:]),0)\
                  - originalDemArrayMatrix)/pixelSize
            Is = (np.concatenate((originalDemArrayMatrix[1:Ny,:],originalDemArrayMatrix[Ny-1,:]),0)\
                  - originalDemArrayMatrix)/pixelSize
            Ie = (np.concatenate((originalDemArrayMatrix[:,1:Nx],originalDemArrayMatrix[:,Nx-1]),1)\
                  - originalDemArrayMatrix)/pixelSize
            Iw = (np.concatenate((originalDemArrayMatrix[:,0],originalDemArrayMatrix[:,0:Nx-1]),1)\
                  - originalDemArrayMatrix)/pixelSize
            In[np.isnan(In)] = 0
            Is[np.isnan(Is)] = 0
            Ie[np.isnan(Ie)] = 0
            Iw[np.isnan(Iw)] = 0
            demArray=originalDemArray

        part6 = np.array(Cn)*np.array(In) + np.array(Cs)*np.array(Is) +\
                np.array(Ce)*np.array(In) + np.array(Cw)*np.array(Iw)
        demArrayMatrix = demArrayMatrix + diffusionTimeIncrement*(part6)
        #print demArrayMatrix.shape
        
        demArrayMatrix = demArray
        #print demArrayMatrix.shape
        
    return demArrayMatrix


def compute_dem_curvature(demArray,pixelDemScale,curvatureCalcMethod):
    print 'computing DTM curvature'
    #demArray[demArray<0]=np.nan
    gradXArray,gradYArray = np.gradient(demArray,pixelDemScale)
    slopeArrayT = np.sqrt(gradXArray**2 + gradYArray**2)
    if curvatureCalcMethod=='geometric':
        #Geometric curvature
        print 'using geometric curvature'
        gradXArrayT = np.divide(gradXArray,slopeArrayT)
        gradYArrayT = np.divide(gradYArray,slopeArrayT)
        #gradXArrayT[slopeArrayT==0.0]=0.0
        #gradYArrayT[slopeArrayT==0.0]=0.0
    elif curvatureCalcMethod=='laplacian':
        # do nothing..
        print 'using laplacian curvature'
        gradXArrayT = gradXArray
        gradYArrayT = gradYArray
    
    gradGradXArray,tmpy = np.gradient(gradXArrayT,pixelDemScale)
    tmpX,gradGradYArray = np.gradient(gradYArrayT,pixelDemScale)
    curvatureDemArray = gradGradXArray + gradGradYArray
    del tmpy, tmpX
    return curvatureDemArray

def compute_quantile_quantile_curve(x):
    print 'getting qqplot estimate'
    defaults.figureNumber = defaults.figureNumber + 1
    plt.figure(defaults.figureNumber)
    res = stats.probplot(x, plot=plt)
    res1 = sm.ProbPlot(x, stats.t, fit=True)
    print res1
    return res

# Write geotif to file on a disk
def write_geotif_generic(inputArray,outfilepath,outfilename):
    print 'writing geotiff', outfilename
    output_fileName = outfilepath + outfilename
    # Get shape
    ncols = inputArray.shape[0]
    nrows = inputArray.shape[1]
    # create the output image
    driver = Parameters.driver
    #print driver
    outDs = driver.Create(output_fileName, nrows, ncols, 1, gdal.GDT_Float64)
    if outDs is None:
        print 'Could not create reclass_40.tif'
        sys.exit(1)
    outBand = outDs.GetRasterBand(1)
    #outData = inputArray
    #outBand.SetNoDataValue(-3.402823e+038)
    
    # set the reference info
    geotransform = Parameters.geotransform
    cc = (geotransform[0],geotransform[1],geotransform[2],\
          geotransform[3],geotransform[4],geotransform[5])
    outDs.SetGeoTransform(cc)
    outDs.SetProjection(Parameters.inputwktInfo)

    # write the band
    PMarray=np.array(inputArray)
    outBand.WriteArray(PMarray)
    # flush data to disk, set the NoData value and calculate stats
    outBand.FlushCache()
    
    del PMarray, outDs, outBand, driver

    

# Write geotiff to disk
def write_geotif_filteredDEM(filteredDemArray,filepath,filename):
    print 'writing filtered DEM'
    output_fileName = "C:\\Mystuff\\grassgisdatabase\\PM_filtered_grassgis.tif"
    # Create gtif
    ncols = filteredDemArray.shape[0]
    nrows = filteredDemArray.shape[1]
    print ncols, nrows
    # create the output image
    driver = gdal.GetDriverByName('GTiff')
    #print driver
    outDs = driver.Create(output_fileName,nrows,ncols,1, gdal.GDT_Float32)
    if outDs is None:
        print 'Could not create tif file'
        sys.exit(1)
    # set the reference info
    geotransform = Parameters.geotransform
    outDs.SetGeoTransform(geotransform)
    outDs.SetProjection(Parameters.inputwktInfo)
    # write the band
    #PMarray=np.array(filteredDemArray)
    outband = outDs.GetRasterBand(1)
    PMarray = np.array(filteredDemArray)
    print type(PMarray)
    outband.WriteArray(PMarray)
    outRasterSRS = osr.SpatialReference(wkt=Parameters.inputwktInfo)
    authoritycode = outRasterSRS.GetAuthorityCode("PROJCS")
    outRasterSRS.ImportFromEPSG(int(authoritycode))
    outDs.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()
    #stop
    

    """
    # Sine the above raster is not north up, we need to write it
    # again so that grass gis can read it properly
    # this needs to be done for the filtered DEM only
    fullFilePath = output_fileName
    dataset1 = gdal.Open(fullFilePath, gdal.GA_ReadOnly)
    output_fileName1 = "C:\\Mystuff\\grassgisdatabase\\PM_filtered_grassgis.tif"
    ncols = dataset1.RasterXSize
    nrows = dataset1.RasterYSize
    # create the output image
    driver = dataset1.GetDriver()

    #print driver
    outDs = driver.Create(output_fileName1,ncols,nrows,1, gdal.GDT_Float32)
    if outDs is None:
        print 'Could not create reclass_40.tif'
        sys.exit(1)
    outBand = outDs.GetRasterBand(1)
    outData = dataset1.GetRasterBand(1).ReadAsArray()
    nanDemArray1=np.array(outData.T)
    # flush data to disk, set the NoData value and calculate stats
    outBand.FlushCache()
    outBand.SetNoDataValue(-9999)

    # georeference the image and set the projection
    geotransform = dataset1.GetGeoTransform()
    cc1 = (geotransform[0],geotransform[1],geotransform[2],\
          geotransform[3],geotransform[4],-geotransform[5])
    outDs.SetGeoTransform(cc1)
    outDs.SetProjection(dataset1.GetProjection())

    # write the data as per grass gis
    outBand.WriteArray(nanDemArray1)
    #"""
    # finishing the writing of filtered DEM
    del outDs, outband, driver,outRasterSRS
    #del dataset, outData,outDs,outBand

    
# Flow accumulation is computed by calling GRASS GIS functions.
def flowaccumulation(filteredDemArray):
    ncols = filteredDemArray.shape[0]
    nrows = filteredDemArray.shape[1]
    print ncols,nrows
    gisbase = os.environ['GISBASE']
    gisdbdir = 'C:\\Users\\Harish\\Documents\\grassdata'
    originalGeotiff = Parameters.demDataFilePath + Parameters.demFileName
    #geotiff = originalGeotiff
    geotiff = 'c:\\Mystuff\\grassgisdatabase\\PM_filtered_grassgis.tif'
    locationGeonet = 'geonet'
    mapsetGeonet = 'geonetuser'
    print 'Making the geonet location'
    print g.run_command('g.proj', georef=geotiff,\
                    location= locationGeonet)
    location = locationGeonet  
    mapset = mapsetGeonet
    print 'Number of Mapsets after making locations'
    print g.read_command('g.mapsets', flags = 'l')
    print 'Setting environ'
    print gsetup.init(gisbase, gisdbdir, locationGeonet, 'PERMANENT')
    print g.gisenv()
    print 'Making mapset now'
    print g.run_command('g.mapset', flags = 'c', mapset = mapsetGeonet,\
                    location = locationGeonet, gisdb = gisdbdir)
    #after adding new mapset
    print 'mapsets after making new'
    print g.read_command('g.mapsets', flags = 'l')
    # gsetup initialization 
    print gsetup.init(gisbase, gisdbdir, locationGeonet, mapsetGeonet)
    # Read the filtered DEM
    print 'r.in.gdal'
    tmpfile = Parameters.demFileName # this reads something like skunkroi.tif
    geotiffmapraster = tmpfile.split('.')[0]
    print 'geotiffmapraster: ',geotiffmapraster
    print g.run_command('r.in.gdal', input=geotiff, \
                        output=geotiffmapraster,overwrite=True)
    gtf = Parameters.geotransform
    print gtf
    """
    # May be not required!
    # Set up grass region
    print 'g.region'
    north = gtf[3]
    south = gtf[3]-nrows
    west = gtf[0]
    east = gtf[0]+ncols   #flags = 'p' 
    print g.run_command('g.region', \
              n = north ,s = south, e = east, w = west,\
              res = 1, rows = nrows ,cols = ncols)
    #"""
    #Flow computation for massive grids (float version)
    print "Calling the r.watershed command from GRASS GIS"
    subbasinThreshold = defaults.thresholdAreaSubBasinIndexing
    """
    print g.run_command('r.watershed',overwrite=True,\
        elevation=geotiffmapraster,\
        threshold=subbasinThreshold, accumulation='acc1v23',basin = 'bas1v23',\
        drainage = 'dra1v23')
        flags ='a',
    """
    print g.run_command('r.watershed',flags ='a',overwrite=True,\
        elevation=geotiffmapraster,\
        threshold=subbasinThreshold, accumulation='acc1v23',\
        drainage = 'dra1v23')
    
    """
    if the size of input data > 5000:
        then call r.watershed with only one input multiple
        times
    
    print g.run_command('r.watershed',overwrite=True,\
        elevation=geotiffmapraster,\
        threshold=subbasinThreshold, accumulation='acc1v23',basin = 'bas1v23',\
        drainage = 'dra1v23')
    
    print g.run_command('r.watershed',overwrite=True,\
        elevation=geotiffmapraster,\
        threshold=subbasinThreshold, accumulation='acc1v23',basin = 'bas1v23',\
        drainage = 'dra1v23')
    """
    print 'r.maplac'
    print g.run_command('r.mapcalc',overwrite=True,\
        expression='outletmap = dra1v23 < 0')
    outlet_filename = geotiffmapraster + '_outlets.tif'
    print g.run_command('r.out.gdal',overwrite=True,\
                        input='outletmap', type='Float32',\
                        output='C:\\Mystuff\\grassgisdatabase\\' +\
                        outlet_filename,\
                        format='GTiff')
    
    # Number of rasters after computation
    print g.read_command('g.list', _type = 'rast')

    # Save the outputs as TIFs
    outputFAC_filename = geotiffmapraster + '_fac.tif'
    print g.run_command('r.out.gdal',overwrite=True,\
                        input='acc1v23', type='Float64',\
                        output='C:\\Mystuff\\grassgisdatabase\\' +\
                        outputFAC_filename,\
                        format='GTiff')
    
    outputFDR_filename = geotiffmapraster + '_fdr.tif'
    print g.run_command('r.out.gdal',overwrite=True,\
                        input = "dra1v23", type='Float64',\
                        output='C:\\Mystuff\\grassgisdatabase\\'+\
                        outputFDR_filename,\
                        format='GTiff')
    """
    outputBAS_filename = geotiffmapraster + '_basins.tif'
    print g.run_command('r.out.gdal',overwrite=True,\
                        input = "bas1v23", type='Float64',\
                        output='C:\\Mystuff\\grassgisdatabase\\'+\
                        outputBAS_filename,\
                        format='GTiff')
    """
    # plot the flow directions
    fdrtif = 'C:\\Mystuff\\grassgisdatabase\\'+outputFDR_filename
    dsfdr = gdal.Open(fdrtif, gdal.GA_ReadOnly)
    aryfdr = dsfdr.GetRasterBand(1).ReadAsArray()
    nanDemArrayfdr=np.array(aryfdr)
    nanDemArrayfdrT = nanDemArrayfdr.T
    del dsfdr,aryfdr

    outlettif = 'C:\\Mystuff\\grassgisdatabase\\'+outlet_filename
    dsout = gdal.Open(outlettif, gdal.GA_ReadOnly)
    aryfdrout = dsout.GetRasterBand(1).ReadAsArray()
    nanDemArrayfdrout=np.array(aryfdrout)
    del dsout,aryfdrout

    outletfromtif = np.where(nanDemArrayfdrout==1)
    print 'outletfromtif'
    print outletfromtif
    """
    Output drainage raster map contains drainage direction.
    Provides the "aspect" for each cell measured CCW from East.
    Multiplying positive values by 45 will give the direction
    in degrees that the surface runoff will travel from that cell.
    The value 0 (zero) indicates that the cell is a depression area
    (defined by the depression input map).
    
    Negative values indicate that surface runoff is leaving the boundaries
    of the current geographic region. The absolute value of these
    negative cells indicates the direction of flow.
    
    """ 
    #outlets = np.where((nanDemArrayfdr<0) & (nanDemArrayfdr!=-3.402823e+038))
    outlets = np.where(nanDemArrayfdrT<0)
    print "Number of outlets :", str(len(outlets[0]))    
    defaults.figureNumber = defaults.figureNumber + 1
    plt.figure(defaults.figureNumber)
    plt.imshow(nanDemArrayfdr)
    plt.xlabel('X[m]')
    plt.ylabel('Y[m]')
    plt.title('Flow directions DEM')
    if defaults.doPlot==1:
        plt.show()

    # plot the flow accumulation
    factif = 'C:\\Mystuff\\grassgisdatabase\\'+outputFAC_filename
    dsfac = gdal.Open(factif, gdal.GA_ReadOnly)
    aryfac = dsfac.GetRasterBand(1).ReadAsArray()
    nanDemArrayfac=np.array(aryfac)
    del dsfac,aryfac
    defaults.figureNumber = defaults.figureNumber + 1
    plt.figure(defaults.figureNumber)
    plt.imshow(nanDemArrayfac,cmap=cm.BrBG)
    plt.xlabel('X[m]')
    plt.ylabel('Y[m]')
    plt.title('Flow accumulations DEM')
    if defaults.doPlot==1:
        plt.show()

    
    # outlets locations in projection of the input dataset
    print outlets
    outletsxx = outlets[0]
    outletsxxfloat = [float(x)+0.5 for x in outletsxx]
    outletsyy = outlets[1]
    outletsyyfloat = [float(x)+0.5 for x in outletsyy]
    """
    outletsxxProj = np.array(outletsxxfloat) * Parameters.demPixelScale + \
                    Parameters.xLowerLeftCoord + float(0.0164)
    outletsyyProj = Parameters.yLowerLeftCoord - np.array(outletsyyfloat) * \
                    Parameters.demPixelScale + \
                    Parameters.yDemSize * Parameters.demPixelScale + float(0.0155)
    
    # The extra decimal digits is essentially a hack into
    # Grass GIS r.water.outlet routine, which only, works
    # with atleast 4 significant digits
    """
    print gtf
    outletsxxProj = float(gtf[0])+ \
                    float(gtf[1]) * np.array(outletsxxfloat)
    
    outletsyyProj = float(gtf[3])+ \
                    float(gtf[5])*np.array(outletsyy)
    
    #"""
    #print outletsxxProj,outletsyyProj
    
    # Call the watershed outlet grass gis function to find the
    # basin Index that will be used for FM marching
    #"""
    #outletFac = []
    print 'using: r.water.outlet'
    nanDemArrayfacT = nanDemArrayfac.T
    for op in range(0,len(outletsxxProj)):
        east = outletsxxProj[op]
        north = outletsyyProj[op]
        outletfac = nanDemArrayfacT[outletsxx[op],outletsyy[op]]
        print '#:',str(op),'east :',east,' north :',north,' outletfac: ',outletfac
        print g.run_command('r.water.outlet',overwrite=True,\
                    input = 'dra1v23', output = 'oneoutletbasin'+str(op+ 1),\
                    coordinates=[east,north])
        # done making big basins
        outputoneBAS_filename = geotiffmapraster + '_onebasins_'+str(op+ 1)+'.tif'
        print g.run_command('r.out.gdal',overwrite=True,\
                        input = "oneoutletbasin"+str(op+ 1), type='Float32',\
                        output='C:\\Mystuff\\grassgisdatabase\\basinTiffs\\'+\
                        outputoneBAS_filename,\
                        format='GTiff')

    # print g.read_command('g.list', _type='rast')
    allbasins  = np.zeros((nanDemArrayfacT.shape))
    #print outletFac    
    
    # Read big basin files
    for op in range(1,len(outletsxxProj)):
        bigbasintif = 'C:\\Mystuff\\grassgisdatabase\\basinTiffs\\'+\
                      geotiffmapraster+'_onebasins_'+str(op)+'.tif'
        dsbasin = gdal.Open(bigbasintif, gdal.GA_ReadOnly)
        arybasin = dsbasin.GetRasterBand(1).ReadAsArray()
        nanDemArrayBasins=np.array(arybasin)
        nanDemArrayBasinsT = nanDemArrayBasins.T
        allbasins[nanDemArrayBasinsT==1]=op
        dsbasin = None
    
    #"""
    return {'outlets':outlets, 'fac':nanDemArrayfac ,\
            'fdr':nanDemArrayfdr,\
            'outletsxxProj':outletsxxProj, 'outletsyyProj':outletsyyProj,\
            'bigbasins':allbasins}
    """
    # return back the dictionary of flow routing
    return {'outlets':outlets, 'fac':nanDemArrayfac ,\
            'fdr':nanDemArrayfdr ,'basins':nanDemArraybasins,\
            'outletsxxProj':outletsxxProj, 'outletsyyProj':outletsyyProj,\
            'bigbasins':allbasins}
    """
    # end of flow accumulation


# Skeleton by thresholding one grid measure e.g. flow or curvature
def compute_skeleton_by_single_threshold(inputArray, threshold):
    skeletonArray = np.zeros((inputArray.shape))
    #skeletonArray = skeletonArray.T
    skeletonArray[np.where(inputArray> threshold)] = 1
    return skeletonArray

# Skeleton by thresholding two grid measures e.g. flow and curvature
def compute_skeleton_by_dual_threshold(inputArray1, inputArray2, threshold1, threshold2):
    skeletonArray = np.zeros((inputArray1.shape))
    mask1 = np.where(inputArray1> threshold1,1,False)
    mask2 = np.where(inputArray2>threshold2,1,False)
    skeletonArray= mask1*mask2
    return skeletonArray

# Normalize curvature
def normalize(inputArray):
    normalizedArray = inputArray- np.min(inputArray[~np.isnan(inputArray)])
    normalizedArrayR = normalizedArray/ np.max(normalizedArray[~np.isnan(normalizedArray)])
    return normalizedArrayR

# Compute discrete geodesics
def compute_discrete_geodesic(geodesicDistanceArray,skeletonEndPoint,doTrueGradientDescent ):
    #print 'computing discrete geodesics'
    # Extract a discrete geodesic path in 2D
    # D = geodesic distance matrix
    # x = channel head or start point
    # path = variable that stores the pixel values of the stream line.
    skeletonEndPoint = skeletonEndPoint[:]
    #print skeletonEndPoint[:]
    streamPathPixelList = skeletonEndPoint[:]
    print 'skeletonEndPoint',skeletonEndPoint
    # Creating the 8 cell neighbor moves
    tempArrayDxMoves = [1, -1, 0, 0, 1, -1, 1, -1]
    tempArrayDyMoves = [0, 0, 1, -1, 1, -1, -1, 1]
    tempArray = [tempArrayDxMoves,tempArrayDyMoves]
    # Get the geodesic value for the channel head
    channelHeadGeodesicDistance = geodesicDistanceArray[skeletonEndPoint[0],skeletonEndPoint[1]]
    #print 'channelHeadGeodesicDistance',channelHeadGeodesicDistance
    # Get the size of the geodesic distance
    geodesicDistanceArraySize = geodesicDistanceArray.shape
    #print geodesicDistanceArraySize
    # While we find a geodesic distance less then previous value
    while True:
        cardinalDxMoves = [1, -1, 0, 0]
        cardinalDyMoves = [0, 0, 1, -1]
        diagonalDxMoves = [1, -1, 1, -1]
        diagonalDyMoves = [1, -1, -1, 1]
        cardinalAllPossibleMoves = [cardinalDxMoves,cardinalDyMoves]
        diagonalAllPossibleMoves = [diagonalDxMoves,diagonalDyMoves]
        tempStreamPathPixelList = streamPathPixelList[:,-1]
        #print tempStreamPathPixelList
        tempStreamPathPixelListA = np.array([[tempStreamPathPixelList[0]],\
                                             [tempStreamPathPixelList[1]]])
        cardinalSkeletonEndPoint = np.repeat(tempStreamPathPixelListA,4,axis=1)+\
                                      cardinalAllPossibleMoves
        diagonalSkeletonEndPoint = np.repeat(tempStreamPathPixelListA,4,axis=1)+\
                                   diagonalAllPossibleMoves
        r1 = cardinalSkeletonEndPoint.tolist()[0]
        r2 = cardinalSkeletonEndPoint.tolist()[1]
        r3 = diagonalSkeletonEndPoint.tolist()[0]
        r4 = diagonalSkeletonEndPoint.tolist()[1]

        neighborPixelSkeletonEndPointList = np.array([r1 + r3,r2 + r4])

        r5 = neighborPixelSkeletonEndPointList.tolist()[0]
        r6 = neighborPixelSkeletonEndPointList.tolist()[1]

        # Get the indices which are not on boundary
        cardinalAllowedIndex0 = np.array([cardinalSkeletonEndPoint[0,:] > 0] and\
                               [cardinalSkeletonEndPoint[0,:] < geodesicDistanceArraySize[0]])
        cardinalAllowedIndex1 = np.array([cardinalSkeletonEndPoint[1,:] > 0] and\
                               [cardinalSkeletonEndPoint[1,:] < geodesicDistanceArraySize[1]])        
        cardinalAllowedIndex = cardinalAllowedIndex0 * cardinalAllowedIndex1
        
        diagonalAllowedIndex0 = np.array([diagonalSkeletonEndPoint[0,:] > 0] and \
                               [diagonalSkeletonEndPoint[0,:] < geodesicDistanceArraySize[0]])
        diagonalAllowedIndex1 = np.array([diagonalSkeletonEndPoint[1,:] > 0] and\
                               [diagonalSkeletonEndPoint[1,:] < geodesicDistanceArraySize[1]])
        diagonalAllowedIndex = diagonalAllowedIndex0 * diagonalAllowedIndex1
        
        allAllowedIndex0 = np.array([neighborPixelSkeletonEndPointList[0,:] > 0] and\
                    [neighborPixelSkeletonEndPointList[0,:] < geodesicDistanceArraySize[0]])        
        allAllowedIndex1= np.array([neighborPixelSkeletonEndPointList[1,:] > 0] and\
                    [neighborPixelSkeletonEndPointList[1,:] < geodesicDistanceArraySize[1]])
        allAllowedIndex = allAllowedIndex0 * allAllowedIndex1
        
        #print cardinalAllowedIndex[0]
        #print diagonalAllowedIndex[0]
        #print allAllowedIndex[0]

        # Now remove neighbors that are no boundary
        # build the true false array
        tfCarray = np.array([cardinalAllowedIndex[0],cardinalAllowedIndex[0]])
        tfCarrayMask = np.zeros((tfCarray.shape))
        tfCarrayMask[tfCarray==False]=1
        popinfC = np.where(tfCarray[0,:]==False)
        #print popinfC
        
        tfDarray = np.array([diagonalAllowedIndex[0],diagonalAllowedIndex[0]])
        tfDarrayMask = np.zeros((tfDarray.shape))
        tfDarrayMask[tfDarray==False]=1
        popinfD = np.where(tfDarray[0,:]==False)
        #print popinfD
        
        tfAarray = np.array([allAllowedIndex[0],allAllowedIndex[0]])
        tfAarrayMask = np.zeros((tfAarray.shape))
        tfAarrayMask[tfAarray==False]=1
        popinfA = np.where(tfAarray[0,:]==False)
        #print popinfA
        
        # Now remove the false indices from our neighborhood matrix
        # Now arrange the arrays above
        cardinalSkeletonEndPointAllowed = npma.masked_array(cardinalSkeletonEndPoint,\
                                                            mask=tfCarrayMask)
        diagonalSkeletonEndPointAllowed = npma.masked_array(diagonalSkeletonEndPoint,\
                                                            mask=tfDarrayMask)
        neighborPixelSkeletonEndPointListAllowed=npma.masked_array(neighborPixelSkeletonEndPointList,\
                                                            mask=tfAarrayMask)

        rw1 = neighborPixelSkeletonEndPointListAllowed[0,:]
        rw2 = neighborPixelSkeletonEndPointListAllowed[1,:]
        rw3 = cardinalSkeletonEndPointAllowed[0,:]
        rw4 = cardinalSkeletonEndPointAllowed[1,:]
        rw5 = diagonalSkeletonEndPointAllowed[0,:]
        rw6 = diagonalSkeletonEndPointAllowed[1,:]
        # Get the minimum value of geodesic distance in the 8 cell neighbor
        # Get the values of D(I) and adjust values for diagonal elements
        try:
            allGeodesicDistanceList = np.array(geodesicDistanceArray[rw1[~rw1.mask],\
                rw2[~rw2.mask]])
            # new line   
            cardinalPixelGeodesicDistanceList = np.array(geodesicDistanceArray[rw3[~rw3.mask],\
                rw4[~rw4.mask]])
            diagonalPixelGeodesicDistanceList= np.array(geodesicDistanceArray[rw5[~rw5.mask],\
                rw6[~rw6.mask]])
        except:
            print neighborPixelSkeletonEndPointList
            print allAllowedIndex
            print allGeodesicDistanceList
            print popinfC
            print popinfD
            print popinfA
            print rw1,rw2,rw3,rw4,rw5,rw6
            print rw1[~rw1.mask]
            print rw2[~rw2.mask]
        
        #print allGeodesicDistanceList
        #print cardinalPixelGeodesicDistanceList
        #print diagonalPixelGeodesicDistanceList
        # We have to insert np.nan values for masked values
        allFinal = np.zeros((1,8))
        #print popinfA
        allFinal[0,popinfA[0]]= np.nan
        aF = 0
        cardinalFinal = np.zeros((1,4))
        #print popinfC
        cardinalFinal[0,popinfC[0]]= np.nan
        cF = 0
        diagonalFinal = np.zeros((1,4))
        #print popinfD
        diagonalFinal[0,popinfD[0]] = np.nan
        dF = 0

        #print allFinal,cardinalFinal,diagonalFinal
        for aFi in xrange(0,8):
            if ~np.isnan(allFinal[0,aFi]):
                allFinal[0,aFi] = allGeodesicDistanceList[aF]
                aF = aF+1
        #--------
        for cFi in xrange(0,4):
            if ~np.isnan(cardinalFinal[0,cFi]):
                cardinalFinal[0,cFi] = cardinalPixelGeodesicDistanceList[cF]
                cF = cF+1
        #--------
        for dFi in xrange(0,4):
            if ~np.isnan(diagonalFinal[0,dFi]):
                diagonalFinal[0,dFi] = diagonalPixelGeodesicDistanceList[dF]
                dF = dF+1
        #--------
        del allGeodesicDistanceList,  cardinalPixelGeodesicDistanceList,\
            diagonalPixelGeodesicDistanceList

        allGeodesicDistanceList = allFinal
        cardinalPixelGeodesicDistanceList = cardinalFinal
        diagonalPixelGeodesicDistanceList = diagonalFinal
        #print allGeodesicDistanceList
        #print cardinalPixelGeodesicDistanceList
        #print diagonalPixelGeodesicDistanceList
        #stop
        # for cells in horizontal and vertical positions to the
        # current cell
        cardinalPixelGeodesicDistanceList = channelHeadGeodesicDistance - \
                                            cardinalPixelGeodesicDistanceList
        # for cells in the diagonal position to the current cell
        diagonalPixelGeodesicDistanceList = (channelHeadGeodesicDistance - \
                                            diagonalPixelGeodesicDistanceList)/np.sqrt(2)
        tcL = cardinalPixelGeodesicDistanceList.tolist()
        tdL = diagonalPixelGeodesicDistanceList.tolist()
        neighborPixelGeodesicDistanceList = np.array(tcL[0]+tdL[0])
        
        #print type(cardinalPixelGeodesicDistanceList)
        #print tcL[0] + tdL[0]
        #print diagonalPixelGeodesicDistanceList
        #print neighborPixelGeodesicDistanceList
        
        # get the index of the maximum geodesic array
        chosenGeodesicIndex = np.argmax(neighborPixelGeodesicDistanceList)
        #print 'chosenGeodesicIndex',chosenGeodesicIndex
        # This is required to break out of the while loop
        chosenGeodesicDistanceFromAll = np.amin(allGeodesicDistanceList)
        #print 'neighborPixelSkeletonEndPointList',neighborPixelSkeletonEndPointList
        neighborPixelSkeletonEndPointList = neighborPixelSkeletonEndPointList[:,chosenGeodesicIndex]
        #print neighborPixelSkeletonEndPointList
        #stop
        if chosenGeodesicDistanceFromAll > channelHeadGeodesicDistance:
            #print "greater geo distance"
            #print channelHeadGeodesicDistance
            break
        elif chosenGeodesicDistanceFromAll == 0:
            print "equal zzero"
            break
        #print 'before assig:',channelHeadGeodesicDistance
        channelHeadGeodesicDistance = chosenGeodesicDistanceFromAll
        #print 'afetr assig:',channelHeadGeodesicDistance
        #print channelHeadGeodesicDistance
        # Finally add the value of neighborPixelSkeletonEndPointList
        # to path list
        b = np.array([[neighborPixelSkeletonEndPointList[0]],\
                      [neighborPixelSkeletonEndPointList[1]]])
        #print 'b',b
        streamPathPixelList = np.hstack((streamPathPixelList,b))
    #print 'streamPathPixelList',streamPathPixelList
    #stop
    return streamPathPixelList
        
# Writing channel head shapefiles
def write_channel_heads(xx,yy):
    print "Writing Channel Heads shapefile"
    
    # set up the shapefile driver
    driver = ogr.GetDriverByName(Parameters.driverName)
    # This will delete and assist in overwrite of the shape files
    if os.path.exists(Parameters.FileName):
        driver.DeleteDataSource(Parameters.FileName)
    
    # create the data source
    data_source = driver.CreateDataSource(Parameters.FileName)
    
    # create the spatial reference, same as the input dataset
    srs = osr.SpatialReference()
    gtf = Parameters.geotransform
    georef = Parameters.inputwktInfo
    
    srs.ImportFromWkt(georef)
    
    # Project the xx, and yy points
    xxProj = float(gtf[0])+ \
                    float(gtf[1]) * np.array(xx)
    yyProj = float(gtf[3])+ \
                    float(gtf[5])*np.array(yy)
    
    
    # create the layer
    layer = data_source.CreateLayer(Parameters.shapefileName,\
                                    srs, ogr.wkbPoint)
    # Add the fields we're interested in
    field_name = ogr.FieldDefn("Name", ogr.OFTString)
    field_name.SetWidth(24)
    layer.CreateField(field_name)
    
    field_region = ogr.FieldDefn("Region", ogr.OFTString)
    field_region.SetWidth(24)
    layer.CreateField(field_region)
    
    layer.CreateField(ogr.FieldDefn("Latitude", ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn("Longitude", ogr.OFTReal))
    tmpfname = Parameters.demFileName
    # Now add the channel heads as features to the layer
    for i in xrange(0,len(xxProj)):
        # create the feature
        feature = ogr.Feature(layer.GetLayerDefn())
        # Set the attributes using the values
        feature.SetField("Name", 'ChannelHead')        
        feature.SetField("Region", tmpfname.split(".")[0])
        feature.SetField("Latitude", xxProj[i])
        feature.SetField("Longitude", yyProj[i])
        # create the WKT for the feature using Python string formatting
        wkt = "POINT(%f %f)" %  (float(xxProj[i]) , float(yyProj[i]))
        # Create the point from the Well Known Txt
        point = ogr.CreateGeometryFromWkt(wkt)
        # Set the feature geometry using the point
        feature.SetGeometry(point)
        # Create the feature in the layer (shapefile)
        layer.CreateFeature(feature)
        # Destroy the feature to free resources
        feature.Destroy()

    # Destroy the data source to free resources
    data_source.Destroy()



#---------------------------------------------------------------------------------
#------------------- MAIN FUNCTION--------------------------------------------------
#---------------------------------------------------------------------------------

def main():
    print "current working directory", os.getcwd()
    print "Reading input file path :",Parameters.demDataFilePath
    print "Reading input file :",Parameters.demFileName
    defaults.figureNumber = 0

    rawDemArray = read_dem_from_geotiff(Parameters.demFileName,\
                                        Parameters.demDataFilePath)

    nanDemArraylr=np.array(rawDemArray)
    nanDemArray = nanDemArraylr
    nanDemArray[nanDemArray < defaults.demNanFlag]= np.nan
    Parameters.minDemValue= np.min(nanDemArray[:])
    Parameters.maxDemValue= np.max(nanDemArray[:])

    defaults.figureNumber = defaults.figureNumber + 1
    plt.figure(defaults.figureNumber)
    plt.imshow(nanDemArray,cmap=cm.coolwarm)
    plt.xlabel('X[m]')
    plt.ylabel('Y[m]')
    plt.title('Input DEM')
    if defaults.doPlot==1:
        plt.show()

    # Area of analysis
    Parameters.xDemSize=np.size(nanDemArray,0)
    Parameters.yDemSize=np.size(nanDemArray,1)

    # Calculate pixel length scale and assume square
    Parameters.maxLowerLeftCoord = np.max([Parameters.xDemSize, Parameters.yDemSize])
    print 'DTM size: ',Parameters.xDemSize, 'x' ,Parameters.yDemSize
    #-----------------------------------------------------------------------------

    # Compute slope magnitude for raw and filtered DEMs
    print 'Computing slope of raw DTM'
    print 'DEM pixel scale:',Parameters.demPixelScale
    print np.array(nanDemArray).shape
    slopeXArray,slopeYArray = np.gradient(np.array(nanDemArray),Parameters.demPixelScale)
    slopeMagnitudeDemArray = np.sqrt(slopeXArray**2 + slopeYArray**2)

    # plot the slope DEM array
    slopeMagnitudeDemArrayNp = np.array(slopeMagnitudeDemArray)
    print slopeMagnitudeDemArrayNp.shape

    # plotting the slope DEM of non filtered DEM
    defaults.figureNumber = defaults.figureNumber + 1
    plt.figure(defaults.figureNumber)
    plt.imshow(slopeMagnitudeDemArrayNp,cmap=cm.coolwarm)
    plt.xlabel('X[m]')
    plt.ylabel('Y[m]')
    plt.title('Slope of unfiltered DEM')
    if defaults.doPlot==1:
        plt.show()

    # Computation of the threshold lambda used in Perona-Malik nonlinear
    # filtering. The value of lambda (=edgeThresholdValue) is given by the 90th
    # quantile of the absolute value of the gradient.
    print'Computing lambda = q-q-based nonlinear filtering threshold'
    slopeMagnitudeDemArrayQ = slopeMagnitudeDemArrayNp
    slopeMagnitudeDemArrayQ = np.reshape(slopeMagnitudeDemArrayQ,\
                                         np.size(slopeMagnitudeDemArrayQ))
    slopeMagnitudeDemArrayQ = slopeMagnitudeDemArrayQ[~np.isnan(slopeMagnitudeDemArrayQ)]
    print 'dem smoothing Quantile',defaults.demSmoothingQuantile

    edgeThresholdValue = quantile(np.absolute(slopeMagnitudeDemArrayQ),\
                                  defaults.demSmoothingQuantile)
    print 'edgeThresholdValue :', edgeThresholdValue

    tempArray = np.reshape(slopeMagnitudeDemArrayNp,np.size(slopeMagnitudeDemArrayNp))
    edgeThresholdValuescipy = mquantiles(np.absolute(slopeMagnitudeDemArrayQ),\
                                         defaults.demSmoothingQuantile)
    print 'edgeThresholdValuescipy :', edgeThresholdValuescipy

    edgeThresholdValueasmatlab = quantileasmatlab(np.absolute(tempArray),\
                                                  defaults.demSmoothingQuantile)
    print 'edgeThresholdValueasmatlab :', edgeThresholdValueasmatlab

    # performing Perona-Malik filtering
    print 'Performing Perona-Malik nonlinear filtering'
    filteredDemArray = geonet_diffusion (nanDemArray,\
    defaults.diffusionMethod,\
    defaults.nFilterIterations, edgeThresholdValuescipy,\
    defaults.diffusionTimeIncrement,\
    defaults.diffusionSigmaSquared, 1+0*Parameters.demPixelScale);

    #np.savetxt('C:\\Mystuff\\grassgisdatabase\\fildem.txt', filteredDemArray, delimiter=',')
    
    # plotting the filtered DEM
    defaults.figureNumber = defaults.figureNumber + 1
    plt.figure(defaults.figureNumber)
    plt.imshow(filteredDemArray,cmap=cm.coolwarm)
    plt.xlabel('X[m]')
    plt.ylabel('Y[m]')
    plt.title('Filtered DEM')
    if defaults.doPlot==1:
        plt.show()
    
    # Writing the filtered DEM as a tif
    write_geotif_filteredDEM(filteredDemArray,Parameters.demDataFilePath,\
                             Parameters.demFileName)

    # Computing slope of filtered DEM
    print 'Computing slope of filtered DTM'
    filteredDemArraynp = filteredDemArray#np.gradient only takes an array as input
    slopeXArray,slopeYArray = np.gradient(filteredDemArraynp,Parameters.demPixelScale)
    slopeDemArray = np.sqrt(slopeXArray**2 + slopeYArray**2)
    #slopeDemArray[np.isnan(slopeDemArray)] = 0.0
    #slopeDemArray[np.isnan(filteredDemArraynp)] = np.nan
    
    slopeMagnitudeDemArrayQ = slopeDemArray
    slopeMagnitudeDemArrayQ = np.reshape(slopeMagnitudeDemArrayQ,\
                                         np.size(slopeMagnitudeDemArrayQ))
    slopeMagnitudeDemArrayQ = slopeMagnitudeDemArrayQ[~np.isnan(slopeMagnitudeDemArrayQ)]
    print ' angle min:', np.arctan(quantile(slopeMagnitudeDemArrayQ,0.001))*180/np.pi
    print ' angle max:', np.arctan(quantile(slopeMagnitudeDemArrayQ,0.999))*180/np.pi
    print 'mean slope:',np.nanmean(slopeDemArray[:])
    print 'stdev slope:',np.nanstd(slopeDemArray[:])
    
    #Computing curvature
    print 'computing curvature'
    curvatureDemArrayIn= filteredDemArraynp
    #curvatureDemArrayIn[curvatureDemArrayIn== defaults.demErrorFlag]=np.nan
    curvatureDemArray = compute_dem_curvature(curvatureDemArrayIn,\
                                              Parameters.demPixelScale,\
                                              defaults.curvatureCalcMethod)
    print curvatureDemArray.shape
    # Writing the curvature array
    # inputArray = curvatureDemArray
    outfilepath = 'C:\\Mystuff\\grassgisdatabase\\'
    outfilename = Parameters.demFileName
    outfilename = outfilename.split('.')[0]+'_curvature.tif'
    write_geotif_generic(curvatureDemArray,outfilepath,outfilename)
    
    #Computation of statistics of curvature
    print 'Computing curvature statistics'
    print curvatureDemArray.shape
    tt = curvatureDemArray[~np.isnan(curvatureDemArray[:])]
    print tt.shape
    finiteCurvatureDemList = curvatureDemArray[np.isfinite(curvatureDemArray[:])]
    print finiteCurvatureDemList.shape
    curvatureDemMean = np.nanmean(finiteCurvatureDemList)
    curvatureDemStdDevn = np.nanstd(finiteCurvatureDemList)
    print ' mean: ', curvatureDemMean
    print ' standard deviation: ', curvatureDemStdDevn


    # plotting only for testing purposes
    defaults.figureNumber = defaults.figureNumber + 1
    plt.figure(defaults.figureNumber)
    plt.imshow(curvatureDemArray,cmap=cm.coolwarm)
    plt.xlabel('X[m]')
    plt.ylabel('Y[m]')
    plt.title('Curvature DEM')
    plt.colorbar()
    if defaults.doPlot==1:
        plt.show()
    
    #stop
    #Compute curvature quantile-quantile curve
    # This seems to take a long time ... is commented for now
    print 'Computing curvature quantile-quantile curve'
    #osm,osr = compute_quantile_quantile_curve(finiteCurvatureDemList)
    #print osm[0]
    #print osr[0]
    thresholdCurvatureQQxx = 1
    # have to add method to automatically compute the thresold
    # .....
    # .....
    
   

    # Computing contributing areas
    print 'Computing upstream accumulation areas using MFD from GRASS GIS'
    """
    return {'outlets':outlets, 'fac':nanDemArrayfac ,\
            'fdr':nanDemArrayfdr ,'basins':nanDemArraybasins,\
            'outletsxxProj':outletsxxProj, 'outletsyyProj':outletsyyProj,\
            'bigbasins':allbasins}
    """
    # Call the flow accumulation function
    flowroutingresults = flowaccumulation(filteredDemArray)

    # Read out the flowroutingresults into appropriate variables
    outletPointsList = flowroutingresults['outlets']
    flowArray = flowroutingresults['fac']
    flowDirectionsArray = flowroutingresults['fdr']
    # These are actually not sub basins, if the basin threshold
    # is large, then you might have as nulls, so best
    # practice is to keep the basin threshold close to 1000
    # default value is 10,000
    #subBasinIndexArray = flowroutingresults['basins']

    
    #subBasinIndexArray[subBasinIndexArray==-9999]=np.nan
    basinIndexArray = flowroutingresults['bigbasins']

    flowArray[np.isnan(filteredDemArray)]=np.nan
    flowMean = np.mean(flowArray[~np.isnan(flowArray[:])])
    print 'Mean upstream flow: ', flowMean

    # plotting only for testing purposes
    defaults.figureNumber = defaults.figureNumber + 1
    plt.figure(defaults.figureNumber)
    drainageMeasure = -np.sqrt(np.log10(flowArray))
    plt.imshow(drainageMeasure,cmap=cm.coolwarm)
    plt.plot(outletPointsList[0],outletPointsList[1],'go')
    plt.xlabel('X[m]')
    plt.ylabel('Y[m]')
    plt.title('flowArray with outlets')
    plt.colorbar()
    if defaults.doPlot==1:
        plt.show()
    #stop
    # plotting only for testing purposes
    defaults.figureNumber = defaults.figureNumber + 1
    plt.figure(defaults.figureNumber)
    plt.imshow(basinIndexArray.T,cmap=cm.Dark2)
    plt.plot(outletPointsList[0],outletPointsList[1],'go')
    plt.xlabel('X[m]')
    plt.ylabel('Y[m]')
    plt.title('basinIndexArray with outlets')
    if defaults.doPlot==1:
        plt.show()


    # Define a skeleton based on flow alone
    skeletonFromFlowArray = \
    compute_skeleton_by_single_threshold(flowArray.T,\
        defaults.flowThresholdForSkeleton)
    
    # Define a skeleton based on curvature alone
    skeletonFromCurvatureArray =\
    compute_skeleton_by_single_threshold(curvatureDemArray.T,\
        curvatureDemMean+thresholdCurvatureQQxx*curvatureDemStdDevn)
    
    
    # Define a skeleton based on curvature and flow
    skeletonFromFlowAndCurvatureArray =\
    compute_skeleton_by_dual_threshold(curvatureDemArray.T, flowArray.T, \
        curvatureDemMean+thresholdCurvatureQQxx*curvatureDemStdDevn, \
        defaults.flowThresholdForSkeleton)

    # plotting only for testing purposes
    defaults.figureNumber = defaults.figureNumber + 1
    plt.figure(defaults.figureNumber)
    plt.imshow(skeletonFromFlowAndCurvatureArray.T,cmap=cm.binary)
    plt.plot(outletPointsList[0],outletPointsList[1],'go')
    plt.xlabel('X[m]')
    plt.ylabel('Y[m]')
    plt.title('Curvetaure with outlets')
    if defaults.doPlot==1:
        plt.show()
    
    # Writing the skeletonFromFlowAndCurvatureArray array
    outfilepath = 'C:\\Mystuff\\grassgisdatabase\\'
    outfilename = Parameters.demFileName
    outfilename = outfilename.split('.')[0]+'_skeleton.tif'
    write_geotif_generic(skeletonFromFlowAndCurvatureArray.T,\
                         outfilepath,outfilename)

    
    # Computing the percentage drainage areas
    print 'Computing percentage drainage area of each indexed basin'
    fastMarchingStartPointList = np.array(outletPointsList)
    fastMarchingStartPointListFMM = np.zeros((fastMarchingStartPointList.shape))
    basinsUsedIndexList = np.zeros((len(fastMarchingStartPointList[0]),1))
    nx = Parameters.xDemSize
    ny = Parameters.yDemSize
    nDempixels = float(nx*ny)
    for label in range(0,len(fastMarchingStartPointList[0])):
        outletbasinIndex = basinIndexArray[fastMarchingStartPointList[0,label],\
                                         fastMarchingStartPointList[1,label]]
        numelments = basinIndexArray[basinIndexArray==outletbasinIndex]
        #print type(numelments), len(numelments)
        percentBasinArea = float(len(numelments)) * 100/nDempixels
        print 'Basin: ',outletbasinIndex,\
              '@ : ',fastMarchingStartPointList[:,label],' #Elements ',len(numelments),\
              ' area ',percentBasinArea,' %'
        if percentBasinArea > defaults.thresholdPercentAreaForDelineation and\
           len(numelments) > 6:
            # Get the watersheds used
            basinsUsedIndexList[label]= label
            # Preparing the outlets used for fast marching in ROI
            fastMarchingStartPointListFMM[:,label] = fastMarchingStartPointList[:,label]
        # finishing Making outlets for FMM
    #Closing Basin area computation
    print fastMarchingStartPointListFMM

    # Computing the local cost function
    print 'Preparing to calculate cost function'
    # lets normalize the curvature first
    if defaults.doNormalizeCurvature ==1:
        curvatureDemArrayNor = normalize(curvatureDemArray)
    del curvatureDemArray
    curvatureDemArray = curvatureDemArrayNor
    del curvatureDemArrayNor
    
    defaults.figureNumber = defaults.figureNumber + 1
    plt.figure(defaults.figureNumber)
    plt.figure(defaults.figureNumber)
    plt.imshow(curvatureDemArray,cmap=cm.coolwarm)
    plt.title('Curvature after normalization')
    plt.colorbar()
    if defaults.doPlot==1:
        plt.show()
    
    
    print 'Curvature min: ' ,str(np.min(curvatureDemArray[~np.isnan(curvatureDemArray)])), \
          ' exp(min): ',str(np.exp(3*np.min(curvatureDemArray[~np.isnan(curvatureDemArray)])))
    print 'Curvature max: ' ,str(np.max(curvatureDemArray[~np.isnan(curvatureDemArray)])),\
          ' exp(max): ',str(np.exp(3*np.max(curvatureDemArray[~np.isnan(curvatureDemArray)])))
    
    # set all the nan's to zeros before cost function is computed
    curvatureDemArray[np.isnan(curvatureDemArray)] = 0
    
    print 'Computing cost function & geodesic distance'
    # Calculate the local reciprocal cost (weight, or propagation speed in the
    # eikonal equation sense).  If the cost function isn't defined, default to
    # old cost function.
    flowArray = flowArray.T
    curvatureDemArray = curvatureDemArray.T
    
    if hasattr(defaults, 'reciprocalLocalCostFn'):
        print 'Evaluating local cost func.'
        reciprocalLocalCostArray = eval(defaults.reciprocalLocalCostFn)
    else:
        print 'Evaluating local cost func. (default)'
        reciprocalLocalCostArray = flowArray + \
                                   (flowMean*skeletonFromFlowAndCurvatureArray)\
                                   + (flowMean*curvatureDemArray)
    del reciprocalLocalCostArray
    # Forcing the evaluations
    reciprocalLocalCostArray = flowArray + \
                                   (flowMean*skeletonFromFlowAndCurvatureArray)\
                                   + (flowMean*curvatureDemArray)
    if hasattr(defaults,'reciprocalLocalCostMinimum'):
        if defaults.reciprocalLocalCostMinimum != 'nan':
            reciprocalLocalCostArray[reciprocalLocalCostArray[:]\
                                 < defaults.reciprocalLocalCostMinimum]=1.0
    
    print '1/cost min: ', np.nanmin(reciprocalLocalCostArray[:]) 
    print '1/cost max: ', np.nanmax(reciprocalLocalCostArray[:])

    # Writing the reciprocal array
    outfilepath = 'C:\\Mystuff\\grassgisdatabase\\'
    outfilename = Parameters.demFileName
    outfilename = outfilename.split('.')[0]+'_costfunction.tif'
    write_geotif_generic(reciprocalLocalCostArray,outfilepath,outfilename)

    # Fast marching
    print 'Performing fast marching'
    print '# of unique basins:',np.size(np.unique(basinIndexArray))
    # Now access each unique basin and get the
    # outlets for it
    basinIndexList = np.unique(basinIndexArray)
    print 'basinIndexList:', str(basinIndexList)
    print reciprocalLocalCostArray.shape
    #stop
    
    # Do fast marching for each sub basin
    geodesicDistanceArray = np.zeros((basinIndexArray.shape))
    geodesicDistanceArray[geodesicDistanceArray==0]=np.Inf
    filteredDemArrayTr = filteredDemArray.T
    defaults.figureNumber = defaults.figureNumber + 1
    for i in range(0,len(fastMarchingStartPointListFMM[0])):
        if fastMarchingStartPointListFMM[0,i] != 0 and\
           fastMarchingStartPointListFMM[1,i] != 0:
            basinIndexList = basinIndexArray[fastMarchingStartPointListFMM[0,i],\
                                         fastMarchingStartPointListFMM[1,i]]
            print 'basin Index:',str(basinIndexList)
            maskedBasin = np.zeros((basinIndexArray.shape))
            maskedBasin[basinIndexArray==basinIndexList]=1
            # For the masked basin get the maximum accumulation are
            # location and use that as an outlet for the basin.
            maskedBasinFAC = np.zeros((basinIndexArray.shape))
            maskedBasinFAC[basinIndexArray==basinIndexList]=\
            flowArray[basinIndexArray==basinIndexList]
            maskedBasinFAC[maskedBasinFAC==0]=np.nan
            # Get the outlet of subbasin
            maskedBasinFAC[np.isnan(maskedBasinFAC)]=0
            #subBasinoutletindices = np.where(maskedBasinFAC==maskedBasinFAC.max())
            # print subBasinoutletindices
            # outlets locations in projection of the input dataset
            outletsxx = fastMarchingStartPointList[0,i]#subBasinoutletindices[0]
            outletsyy = fastMarchingStartPointList[1,i]#subBasinoutletindices[1]
            # call the fast marching here
            phi = -1*np.ones((reciprocalLocalCostArray.shape))
            DEMarray = -1*np.ones((reciprocalLocalCostArray.shape))
            phi[maskedBasinFAC!=0] = reciprocalLocalCostArray[maskedBasinFAC!=0]
            phi[phi==-1]=np.nan
            phi[outletsxx,outletsyy] =-1
            speed = phi
            distancearray = skfmm.distance(1/phi, dx=1)
            travelTimearray = skfmm.travel_time(distancearray, speed, dx=1)
            print travelTimearray.shape
            geodesicDistanceArray[maskedBasin ==1]= travelTimearray[maskedBasin ==1]

            #plt.figure(defaults.figureNumber)
            #plt.imshow(speed.T,cmap=cm.coolwarm)
            #plt.contour(travelTimearray.T,cmap=cm.coolwarm)
            #plt.title('basin Index'+str(basinIndexList))
            #plt.show()
            #stop

    # Plot the geodesic array
    defaults.figureNumber = defaults.figureNumber + 1
    plt.figure(defaults.figureNumber)
    #plt.imshow(np.log10(geodesciDistanceArray),cmap=cm.coolwarm)
    plt.contour(geodesicDistanceArray.T,140,cmap=cm.coolwarm)
    plt.title('Geodesic distance array (travel time)')
    plt.colorbar()
    if defaults.doPlot==0:
        plt.show()
    
    print geodesicDistanceArray.shape
    # Writing the geodesic distance array
    outfilepath = 'C:\\Mystuff\\grassgisdatabase\\'
    outfilename = Parameters.demFileName
    outfilename = outfilename.split('.')[0]+'_geodesicDistance.tif'
    write_geotif_generic(geodesicDistanceArray,outfilepath,outfilename)
    
    # Locating end points
    print 'Locating skeleton end points'
    xySkeletonSize = skeletonFromFlowAndCurvatureArray.shape
    print xySkeletonSize
    skeletonLabeledArray, skeletonNumConnectedComponentsList =\
                          ndimage.label(skeletonFromFlowAndCurvatureArray)
    #print skeletonLabeledArray, skeletonNumConnectedComponentsList
    #print np.unique(skeletonLabeledArray)
    """
     Through the histogram of skeletonNumElementsSortedList (skeletonNumElementsList minus the maximum value which
       corresponds to the largest connected element of the skeleton) we get the
       size of the smallest elements of the skeleton, which will likely
       correspond to small isolated convergent areas. These elements will be
       excluded from the search of end points.
    """
    print 'Counting the number of elements of each connected component'
    skeletonNumElementsList = ndimage.measurements.histogram(skeletonLabeledArray, 1, \
                skeletonNumConnectedComponentsList, \
                bins = skeletonNumConnectedComponentsList,\
                labels=None)
    skeletonNumElementsSortedList = np.sort(skeletonNumElementsList)
    
    histarray,skeletonNumElementsHistogramX=np.histogram(\
        skeletonNumElementsSortedList[1:len(skeletonNumElementsSortedList)-1],
        np.sqrt(len(skeletonNumElementsSortedList)))
    
    # Create skeleton gridded array
    skeletonNumElementsGriddedArray = np.zeros(xySkeletonSize)
    #"""
    for i in range(0,xySkeletonSize[0]):
        for j in range(0,xySkeletonSize[1]):
            #Gets the watershed label for this specified cell and checked in
            #subsequent if statement
            basinIndex = basinIndexArray[i,j]
            if skeletonLabeledArray[i, j] > 0:
                skeletonNumElementsGriddedArray[i,j] = \
                    skeletonNumElementsList[skeletonLabeledArray[i,j]-1]
    
    defaults.figureNumber = defaults.figureNumber + 1
    plt.figure(defaults.figureNumber)
    plt.imshow(skeletonNumElementsGriddedArray.T,cmap=cm.coolwarm)
    plt.title('Skeleton Num elements Array')
    plt.colorbar()
    if defaults.doPlot==1:
        plt.show()
    
    #"""
    # Elements smaller than skeletonNumElementsThreshold are not considered in the
    # skeletonEndPointsList detection
    print skeletonNumElementsHistogramX
    skeletonNumElementsThreshold = skeletonNumElementsHistogramX[3]
    
    print 'skeletonNumElementsThreshold',str(skeletonNumElementsThreshold)
    
    # Scan the array for finding the channel heads
    print 'Continuing to locate skeleton endpoints'
    #"""
    skeletonEndPointsList = []
    for i in range(0,xySkeletonSize[0]):
        for j in range(0,xySkeletonSize[1]):
            #print i,j
            # Skip this pixel if the current point is not a labeled or if the
            # number of connected skeleton elements is too small
            if skeletonLabeledArray[i,j]!=0 \
               and skeletonNumElementsGriddedArray[i,j]>=skeletonNumElementsThreshold:
                # Define search box and ensure it fits within the DTM bounds
                mx = i-1
                px = xySkeletonSize[0]-i
                my = j-1
                py = xySkeletonSize[1]-j
                xMinus = np.min([defaults.endPointSearchBoxSize, mx])
                xPlus  = np.min([defaults.endPointSearchBoxSize, px])
                yMinus = np.min([defaults.endPointSearchBoxSize, my])
                yPlus  = np.min([defaults.endPointSearchBoxSize, py])
                # Extract the geodesic distances geodesicDistanceArray for pixels within the search box
                searchGeodesicDistanceBox = geodesicDistanceArray[i-xMinus:i+xPlus, j-yMinus:j+yPlus]
                # Extract the skeleton labels for pixels within the search box
                searchLabeledSkeletonBox = skeletonLabeledArray[i-xMinus:i+xPlus, j-yMinus:j+yPlus]
                # Look in the search box for skeleton points with the same label
                # and greater geodesic distance than the current pixel at (i,j)
                # - if there are none, then add the current point as a channel head
                v = searchLabeledSkeletonBox==skeletonLabeledArray[i,j]
                v1 = v * searchGeodesicDistanceBox > geodesicDistanceArray[i,j]
                v3 = np.where(np.any(v1==True,axis=0))
                if len(v3[0])==0:
                    skeletonEndPointsList.append([i,j])
    
    # For loop ends here
    skeletonEndPointsListArray = np.array(skeletonEndPointsList)
    xx = skeletonEndPointsListArray[0:len(skeletonEndPointsListArray),0]
    yy = skeletonEndPointsListArray[0:len(skeletonEndPointsListArray),1]
    
    defaults.figureNumber = defaults.figureNumber + 1
    plt.figure(defaults.figureNumber)
    plt.imshow(skeletonFromFlowAndCurvatureArray.T,cmap=cm.binary)
    plt.plot(xx,yy,'or')
    plt.title('Skeleton Num elements Array with channel heads')
    plt.colorbar()
    if defaults.doPlot==1:
        plt.show()             

    defaults.figureNumber = defaults.figureNumber + 1
    plt.figure(defaults.figureNumber)
    plt.imshow(np.log(geodesicDistanceArray.T),cmap=cm.coolwarm)
    plt.plot(xx,yy,'or')
    plt.title('Geodesic distance Array with channel heads')
    plt.colorbar()
    if defaults.doPlot==1:
        plt.show()             
           
    #"""
    # Do compute discrete geodesics
    print 'Computing discrete geodesics'
    geodesicPathsCellList = []
    numberOfEndPoints = len(xx)
    outerbounds = geodesicDistanceArray.shape
    for i in range(0,numberOfEndPoints):
        print 'EndPoint# ',i,'/',numberOfEndPoints
        xEndPoint = xx[i]
        yEndPoint = yy[i]
        if xEndPoint > 0 and xEndPoint < outerbounds[0] and \
           yEndPoint > 0 and yEndPoint < outerbounds[1]:
            skeletonEndPoint = np.array([[xEndPoint],[yEndPoint]]) 
            watershedLabel = basinIndexArray[xEndPoint,yEndPoint]
            print 'watershedLabel',watershedLabel
            watershedIndexList = basinIndexArray == watershedLabel
            geodesicDistanceArrayMask = np.zeros((geodesicDistanceArray.shape))
            geodesicDistanceArrayMask[watershedIndexList]= geodesicDistanceArray[watershedIndexList]
            geodesicDistanceArrayMask[geodesicDistanceArrayMask == 0]= np.Inf
            #print geodesicDistanceArrayMask.shape
            geodesicPathsCellList.append(compute_discrete_geodesic(geodesicDistanceArrayMask,\
                                    skeletonEndPoint,defaults.doTrueGradientDescent))
    #
    #print 'geodesicPathsCellList',geodesicPathsCellList
    defaults.figureNumber = defaults.figureNumber + 1
    plt.figure(defaults.figureNumber)
    plt.imshow(flowDirectionsArray,cmap=cm.coolwarm)
    for pp in range(0,len(geodesicPathsCellList)):
        plt.plot(geodesicPathsCellList[pp][0,:],geodesicPathsCellList[pp][1,:],'-r')
    plt.plot(xx,yy,'og')
    plt.title('Geodesic Array with channel heads and streams')
    plt.colorbar()
    if defaults.doPlot==1:
        plt.show()
    
    # Write shapefiles of channel heads
    write_channel_heads(xx,yy)
    plt.show()

    # Write stream network as shapefiles
    
    

    
if __name__ == '__main__':
    t0 = clock()
    main()
    t1 = clock()
    print "time taken to complete the script is::",t1-t0," seconds"
    print "script complete"



