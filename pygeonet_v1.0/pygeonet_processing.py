#! /usr/bin/env python
# pygeonet_processing.py
# Run this file after setting up folder structure
# in pygeonet_parameters.py

import sys
import os
from osgeo import gdal,osr,ogr
import statsmodels.api as sm
import numpy as np
from time import clock
import pygeonet_prepare as Parameters
import pygeonet_defaults as defaults
from math import modf, floor
import matplotlib.pyplot as plt
from matplotlib import cm 
from scipy.stats.mstats import mquantiles
import scipy.signal as conv2
from scipy import stats
import skfmm
from scipy import ndimage
import numpy.ma as npma
import shutil

# ----------- GRASS GIS SETUP ------------------
#setting up the environment for grass gis access
"""
# The below grass script will work assuming you have installed
# Grass GIS 7 on your machine and that the required environment
# variables are set on windows as required.
This has not been tested on linux machines yet.

"""
sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))
import grass.script as g
import grass.script.setup as gsetup

# -------------------------- FUNCTIONS ----------------------------
def read_dem_from_geotiff(demFileName,demFilePath):
    # Open the GeoTIFF format DEM
    fullFilePath = demFilePath + demFileName
    #print fullFilePath
    ary = []
    gdal.UseExceptions()
    ds = gdal.Open(fullFilePath, gdal.GA_ReadOnly)
    Parameters.driver = ds.GetDriver()
    geotransform = ds.GetGeoTransform()
    Parameters.geotransform = geotransform
    ary = ds.GetRasterBand(1).ReadAsArray()
    Parameters.demPixelScale = float(geotransform[1])
    Parameters.xLowerLeftCoord = float(geotransform[0])
    Parameters.yLowerLeftCoord = float(geotransform[3])
    Parameters.inputwktInfo = ds.GetProjection()
    return ary

def anisodiff(img,niter,kappa,gamma,step=(1.,1.),option=2):
    # initialize output array
    img = img.astype('float32')
    imgout = img.copy()
    # initialize some internal variables
    deltaS = np.zeros_like(imgout)
    deltaE = deltaS.copy()
    NS = deltaS.copy()
    EW = deltaS.copy()
    gS = np.ones_like(imgout)
    gE = gS.copy()
    for ii in xrange(niter):
        # do a simple gaussian smoothing
        #imgout = simple_gaussian_smoothing(imgout,5,\
        #                                   defaults.diffusionSigmaSquared)
        # calculate the diffs
        deltaS[:-1,: ] = np.diff(imgout,axis=0)
        deltaE[: ,:-1] = np.diff(imgout,axis=1)
        if option == 2:
            gS = 1./(1.+(deltaS/kappa)**2.)/step[0]
            gE = 1./(1.+(deltaE/kappa)**2.)/step[1]
        elif option == 1:
            gS = np.exp(-(deltaS/kappa)**2.)/step[0]
            gE = np.exp(-(deltaE/kappa)**2.)/step[1]
        # update matrices
        E = gE*deltaE
        S = gS*deltaS
        # subtract a copy that has been shifted 'North/West' by one
        # pixel. don't as questions. just do it. trust me.
        NS[:] = S
        EW[:] = E
        NS[1:,:] -= S[:-1,:]
        EW[:,1:] -= E[:,:-1]
        # update the image
        imgout += gamma*(NS+EW)
    
    return imgout


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
    outDs = driver.Create(output_fileName, nrows, ncols, 1, gdal.GDT_Float32)
    if outDs is None:
        print 'Could not create DemName.tif'
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

    

# Write filtered geotiff to disk to be used by GRASS GIS
def write_geotif_filteredDEM(filteredDemArray,filepath,filename):
    print 'writing filtered DEM'
    output_fileName = Parameters.pmGrassGISfileName
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
    # finishing the writing of filtered DEM
    del outDs, outband, driver,outRasterSRS
    #del dataset, outData,outDs,outBand

    
# Flow accumulation is computed by calling GRASS GIS functions.
def flowaccumulation(filteredDemArray):
    ncols = filteredDemArray.shape[0]
    nrows = filteredDemArray.shape[1]
    print ncols,nrows
    gisbase = os.environ['GISBASE']
    gisdbdir = Parameters.gisdbdir
    originalGeotiff = Parameters.demDataFilePath + Parameters.demFileName
    geotiff = Parameters.pmGrassGISfileName
    print gsetup.init(gisbase, gisdbdir, 'demolocation', 'PERMANENT')
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
                    location = locationGeonet, dbase = gisdbdir)
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
    
    #Flow computation for massive grids (float version)
    print "Calling the r.watershed command from GRASS GIS"
    subbasinThreshold = defaults.thresholdAreaSubBasinIndexing
    if Parameters.xDemSize > 4000 or Parameters.yDemSize > 4000:
        print ('using swap memory option')
        print g.run_command('r.watershed',flags ='am',overwrite=True,\
                            elevation=geotiffmapraster, \
                            threshold=subbasinThreshold, \
                            drainage = 'dra1v23')
        print g.run_command('r.watershed',flags ='am',overwrite=True,\
                            elevation=geotiffmapraster, \
                            threshold=subbasinThreshold, \
                            accumulation='acc1v23')
    else :
        print g.run_command('r.watershed',flags ='a',overwrite=True,\
                            elevation=geotiffmapraster, \
                            threshold=subbasinThreshold, \
                            accumulation='acc1v23',\
                            drainage = 'dra1v23')
    
    print 'r.maplac'
    print g.run_command('r.mapcalc',overwrite=True,\
        expression='outletmap = if(dra1v23 >= 0,null(),1)')
    print 'r.to.vector'
    print g.run_command('r.to.vect',overwrite=True,\
                        input = 'outletmap', output = 'outletsmapvec',\
                        type='point')
    print "r.stream.basins"
    print g.run_command('r.stream.basins',overwrite=True,\
                        direction='dra1v23',points='outletsmapvec',\
                        basins = 'outletbains')
    # Number of rasters after computation
    print g.read_command('g.list', _type = 'rast')
    # Save the outputs as TIFs
    outlet_filename = geotiffmapraster + '_outlets.tif'
    print g.run_command('r.out.gdal',overwrite=True,\
                        input='outletmap', type='Float32',\
                        output=Parameters.geonetResultsDir +\
                        outlet_filename,\
                        format='GTiff')
    
    outputFAC_filename = geotiffmapraster + '_fac.tif'
    print g.run_command('r.out.gdal',overwrite=True,\
                        input='acc1v23', type='Float64',\
                        output=Parameters.geonetResultsDir +\
                        outputFAC_filename,\
                        format='GTiff')
    
    outputFDR_filename = geotiffmapraster + '_fdr.tif'
    print g.run_command('r.out.gdal',overwrite=True,\
                        input = "dra1v23", type='Float64',\
                        output=Parameters.geonetResultsDir+\
                        outputFDR_filename,\
                        format='GTiff')
    
    outputBAS_filename = geotiffmapraster + '_basins.tif'
    print g.run_command('r.out.gdal',overwrite=True,\
                        input = "outletbains", type='Int16',\
                        output=Parameters.geonetResultsDir+\
                        outputBAS_filename,\
                        format='GTiff')
    
    # plot the flow directions
    fdrtif = Parameters.geonetResultsDir+outputFDR_filename
    dsfdr = gdal.Open(fdrtif, gdal.GA_ReadOnly)
    aryfdr = dsfdr.GetRasterBand(1).ReadAsArray()
    nanDemArrayfdr=np.array(aryfdr)
    nanDemArrayfdrT = nanDemArrayfdr.T
    del dsfdr,aryfdr
    
    outlettif = Parameters.geonetResultsDir+outlet_filename
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
    outlets = np.where(nanDemArrayfdr<0)
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
    factif = Parameters.geonetResultsDir+outputFAC_filename
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

    # getting the bigbasins from the r.streams.basins modules
    basinstif = Parameters.geonetResultsDir+outputBAS_filename
    dsbasins = gdal.Open(basinstif, gdal.GA_ReadOnly)
    arybasins = dsbasins.GetRasterBand(1).ReadAsArray()
    nanDemArraybasins =np.array(arybasins)
    nanDemArraybasins = nanDemArraybasins.T
    del dsbasins,arybasins
    
    # outlets locations in projection of the input dataset
    print outlets
    outletsxx = outlets[0]
    outletsxxfloat = [float(x)+0.5 for x in outletsxx]
    outletsyy = outlets[1]
    outletsyyfloat = [float(x)+0.5 for x in outletsyy]
    """
    # The extra decimal digits is essentially a hack into
    # Grass GIS r.water.outlet routine, which only, works
    # with atleast 4 significant digits
    """
    print gtf
    outletsxxProj = float(gtf[0])+ \
                    float(gtf[1]) * np.array(outletsxxfloat)
    
    outletsyyProj = float(gtf[3])+ \
                    float(gtf[5])*np.array(outletsyy)
    
    return {'outlets':outlets, 'fac':nanDemArrayfac ,\
            'fdr':nanDemArrayfdr,\
            'outletsxxProj':outletsxxProj, 'outletsyyProj':outletsyyProj,\
            'bigbasins':nanDemArraybasins}
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


def compute_discrete_geodesic_v1():
    # this a new version using r.drain to extract discrete goedesics
    gisbase = os.environ['GISBASE']
    gisdbdir = Parameters.gisdbdir
    locationGeonet = 'geonet'
    mapsetGeonet = 'geonetuser'
    print gsetup.init(gisbase, gisdbdir, locationGeonet, mapsetGeonet)
    # Read the filtered DEM
    print 'r.in.gdal'
    outfilepathgeodesic = Parameters.geonetResultsDir
    outfilenamegeodesic = Parameters.demFileName
    outfilenamegeodesic = outfilenamegeodesic.split('.')[0]+'_geodesicDistance.tif'
    inputgeodesictifile = outfilepathgeodesic +'\\'+outfilenamegeodesic
    print 'importing goedesic tif: ',inputgeodesictifile
    print g.run_command('r.in.gdal', input=inputgeodesictifile, \
                        output=outfilenamegeodesic,overwrite=True)
    
    # The maximum number of points is 1024
    # --- have to add a check---
    # -- seems to run for large point shapefiles without fail.
    
    print 'importing channel heads shape file'
    channeheadsshapefileName = Parameters.pointshapefileName
    inputshapefilepath = Parameters.pointFileName
    print g.run_command('v.in.ogr',input = inputshapefilepath,\
                        layer=channeheadsshapefileName,output=channeheadsshapefileName,\
                        geometry='Point')
    
    print 'executing r.drain'
    print g.run_command('r.drain',input=outfilenamegeodesic,\
                        output='discretegeodesicsras',\
                        start_points=channeheadsshapefileName)
    print 'thining the discrete geodesic raster'
    print g.run_command('r.thin',input='discretegeodesicsras',\
                        output='discretegeodesicsrasthin')
    
    print 'converting the raster geodesic to vector map'
    print g.run_command('r.to.vect',input = 'discretegeodesicsrasthin',\
                        output='discretegeovec', type='line')
     
    
    print 'exporting the geodesics as shapefile'
    print g.run_command('v.out.ogr', input= 'discretegeovec',\
                        output=Parameters.drainagelineFileName,\
                        format='ESRI_Shapefile')
    print 'completed discrete geodesics'
    # ---draining algorithm finished


# Writing channel head shapefiles
def write_channel_heads(xx,yy):
    print "Writing Channel Heads shapefile"
    
    # set up the shapefile driver
    driver = ogr.GetDriverByName(Parameters.driverName)
    # This will delete and assist in overwrite of the shape files
    if os.path.exists(Parameters.pointFileName):
        driver.DeleteDataSource(Parameters.pointFileName)
    
    # create the data source
    data_source = driver.CreateDataSource(Parameters.pointFileName)
    
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
    layer = data_source.CreateLayer(Parameters.pointshapefileName,\
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
    
    # Now add the channel heads as features to the layer
    for i in xrange(0,len(xxProj)):
        # create the feature
        feature = ogr.Feature(layer.GetLayerDefn())
        # Set the attributes using the values
        feature.SetField("Name", 'ChannelHead')        
        feature.SetField("Region", Parameters.Region)
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

# Writing drainage paths as shapefile
def write_drainage_paths(geodesicPathsCellList):
    print 'Writing drainage paths'
    #print geodesicPathsCellList
    # set up the shapefile driver
    driver = ogr.GetDriverByName(Parameters.driverName)
    # This will delete and assist in overwrite of the shape files
    if os.path.exists(Parameters.drainagelineFileName):
        driver.DeleteDataSource(Parameters.drainagelineFileName)
    
    # create the data source
    data_source = driver.CreateDataSource(Parameters.drainagelineFileName)
    
    # create the spatial reference, same as the input dataset
    srs = osr.SpatialReference()
    gtf = Parameters.geotransform
    georef = Parameters.inputwktInfo
    
    srs.ImportFromWkt(georef)
    
    # create the layer
    layer = data_source.CreateLayer(Parameters.drainagelinefileName,\
                                    srs, ogr.wkbLineString)
    
    # Add the fields we're interested in
    field_name = ogr.FieldDefn("Name", ogr.OFTString)
    field_name.SetWidth(24)
    layer.CreateField(field_name)
    
    field_region = ogr.FieldDefn("Region", ogr.OFTString)
    field_region.SetWidth(24)
    layer.CreateField(field_region)
    
    layer.CreateField(ogr.FieldDefn("Latitude", ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn("Longitude", ogr.OFTReal))
    
    # Now add the channel heads as features to the layer
    print len(geodesicPathsCellList)
    for i in xrange(0,len(geodesicPathsCellList)):
        #print geodesicPathsCellList[i]
        # Project the linepoints to appropriate projection
        xx = geodesicPathsCellList[i][0]
        yy = geodesicPathsCellList[i][1]
        # Project the xx, and yy points
        xxProj = float(gtf[0])+ \
                    float(gtf[1]) * np.array(xx)
        yyProj = float(gtf[3])+ \
                    float(gtf[5])*np.array(yy)
        # create the feature
        feature = ogr.Feature(layer.GetLayerDefn())
        # Set the attributes using the values
        feature.SetField("Name", 'ChannelNetwork')
        feature.SetField("Region", Parameters.Region)
        # create the WKT for the feature using Python string formatting
        line = ogr.Geometry(ogr.wkbLineString)            
        for j in xrange(0,len(xxProj)):
            #print xxProj[j],yyProj[j]
            line.AddPoint(xxProj[j],yyProj[j])
        #print line
        # Create the point from the Well Known Txt
        #lineobject = line.ExportToWkt()
        # Set the feature geometry using the point
        feature.SetGeometryDirectly(line)
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
    Parameters.maxLowerLeftCoord = np.max([Parameters.xDemSize, \
                                           Parameters.yDemSize])
    print 'DTM size: ',Parameters.xDemSize, 'x' ,Parameters.yDemSize
    #-----------------------------------------------------------------------------

    # Compute slope magnitude for raw and filtered DEMs
    print 'Computing slope of raw DTM'
    print 'DEM pixel scale:',Parameters.demPixelScale
    print np.array(nanDemArray).shape
    slopeXArray,slopeYArray = np.gradient(np.array(nanDemArray),\
                                          Parameters.demPixelScale)
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

    edgeThresholdValuescipy = mquantiles(np.absolute(slopeMagnitudeDemArrayQ),\
                                         defaults.demSmoothingQuantile)
    print 'edgeThresholdValuescipy :', edgeThresholdValuescipy
    
    # performing PM filtering using the anisodiff
    print 'Performing Perona-Malik nonlinear filtering'
    filteredDemArray = anisodiff(nanDemArray, defaults.nFilterIterations, \
                                     edgeThresholdValuescipy,\
                                     defaults.diffusionTimeIncrement, \
                                     (Parameters.demPixelScale,\
                                      Parameters.demPixelScale),2)
    
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
    slopeMagnitudeDemArrayQ = slopeDemArray
    slopeMagnitudeDemArrayQ = np.reshape(slopeMagnitudeDemArrayQ,\
                                         np.size(slopeMagnitudeDemArrayQ))
    slopeMagnitudeDemArrayQ = slopeMagnitudeDemArrayQ[~np.isnan(slopeMagnitudeDemArrayQ)]
    print ' angle min:', np.arctan(np.percentile(slopeMagnitudeDemArrayQ,0.1))*180/np.pi
    print ' angle max:', np.arctan(np.percentile(slopeMagnitudeDemArrayQ,99.9))*180/np.pi
    print 'mean slope:',np.nanmean(slopeDemArray[:])
    print 'stdev slope:',np.nanstd(slopeDemArray[:])
    
    #Computing curvature
    print 'computing curvature'
    curvatureDemArrayIn= filteredDemArraynp
    #curvatureDemArrayIn[curvatureDemArrayIn== defaults.demErrorFlag]=np.nan
    curvatureDemArray = compute_dem_curvature(curvatureDemArrayIn,\
                                              Parameters.demPixelScale,\
                                              defaults.curvatureCalcMethod)
    #Writing the curvature array
    outfilepath = Parameters.geonetResultsDir
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
    
    #*************************************************
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
    #*************************************************
   

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
    plt.plot(outletPointsList[1],outletPointsList[0],'go')
    plt.xlabel('X[m]')
    plt.ylabel('Y[m]')
    plt.title('flowArray with outlets')
    plt.colorbar()
    if defaults.doPlot==1:
        plt.show()
    
    # plotting only for testing purposes
    defaults.figureNumber = defaults.figureNumber + 1
    plt.figure(defaults.figureNumber)
    plt.imshow(basinIndexArray.T,cmap=cm.Dark2)
    plt.plot(outletPointsList[1],outletPointsList[0],'go')
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
    plt.plot(outletPointsList[1],outletPointsList[0],'go')
    plt.xlabel('X[m]')
    plt.ylabel('Y[m]')
    plt.title('Curvature with outlets')
    if defaults.doPlot==1:
        plt.show()
    
    # Writing the skeletonFromFlowAndCurvatureArray array
    outfilepath = Parameters.geonetResultsDir
    outfilename = Parameters.demFileName
    outfilename = outfilename.split('.')[0]+'_skeleton.tif'
    write_geotif_generic(skeletonFromFlowAndCurvatureArray.T,\
                         outfilepath,outfilename)

    
    # Computing the percentage drainage areas
    print 'Computing percentage drainage area of each indexed basin'
    fastMarchingStartPointList = np.array(outletPointsList)
    print fastMarchingStartPointList
    #fastMarchingStartPointListFMM = np.zeros((fastMarchingStartPointList.shape))
    fastMarchingStartPointListFMMx = []
    fastMarchingStartPointListFMMy = []
    basinsUsedIndexList = np.zeros((len(fastMarchingStartPointList[0]),1))
    nx = Parameters.xDemSize
    ny = Parameters.yDemSize
    nDempixels = float(nx*ny)
    basinIndexArray = basinIndexArray.T
    for label in range(0,len(fastMarchingStartPointList[0])):        
        outletbasinIndex = basinIndexArray[fastMarchingStartPointList[0,label],\
                                         fastMarchingStartPointList[1,label]]
        print outletbasinIndex
        numelments = basinIndexArray[basinIndexArray==outletbasinIndex]
        #print type(numelments), len(numelments)
        percentBasinArea = float(len(numelments)) * 100/nDempixels
        print 'Basin: ',outletbasinIndex,\
              '@ : ',fastMarchingStartPointList[:,label],' #Elements ',len(numelments),\
              ' area ',percentBasinArea,' %'
        if percentBasinArea > defaults.thresholdPercentAreaForDelineation and\
           len(numelments) > Parameters.numBasinsElements:
            # Get the watersheds used
            basinsUsedIndexList[label]= label
            # Preparing the outlets used for fast marching in ROI
            #fastMarchingStartPointListFMM[:,label] = fastMarchingStartPointList[:,label]
            fastMarchingStartPointListFMMx.append(fastMarchingStartPointList[0,label])
            fastMarchingStartPointListFMMy.append(fastMarchingStartPointList[1,label])
        # finishing Making outlets for FMM
    #Closing Basin area computation

    fastMarchingStartPointListFMM = np.array([fastMarchingStartPointListFMMx,\
                                                  fastMarchingStartPointListFMMy])
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
    outfilepath = Parameters.geonetResultsDir
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
    geodesicDistanceArray = geodesicDistanceArray.T
    filteredDemArrayTr = filteredDemArray.T
    basinIndexArray = basinIndexArray.T
    # create a watershed outlet dictionary
    outletwatersheddict = {}
    defaults.figureNumber = defaults.figureNumber + 1
    for i in range(0,len(fastMarchingStartPointListFMM[0])):
        basinIndexList = basinIndexArray[fastMarchingStartPointListFMM[1,i],\
                                    fastMarchingStartPointListFMM[0,i]]
        print 'basin Index:',basinIndexList
        print 'start point :', fastMarchingStartPointListFMM[:,i]
        outletwatersheddict[basinIndexList]=fastMarchingStartPointListFMM[:,i]
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
        # print subBasinoutletindices
        # outlets locations in projection of the input dataset
        outletsxx = fastMarchingStartPointList[0,i]
        outletsyy = fastMarchingStartPointList[1,i]
        # call the fast marching here
        phi = np.nan * np.ones((reciprocalLocalCostArray.shape))
        speed = np.ones((reciprocalLocalCostArray.shape))* np.nan
        phi[maskedBasinFAC!=0] = 1
        speed[maskedBasinFAC!=0] = reciprocalLocalCostArray[maskedBasinFAC!=0]
        phi[fastMarchingStartPointListFMM[1,i],\
            fastMarchingStartPointListFMM[0,i]] =-1
        try:
            travelTimearray = skfmm.travel_time(phi,speed, dx=1)
        except IOError as e:            
            print 'Error in calculating skfmm travel time'
            print 'Error in catchment: ',basinIndexList
            # setting travel time to empty array
            travelTimearray = np.nan * np.zeros((reciprocalLocalCostArray.shape))
            plt.figure(defaults.figureNumber+1)
            plt.imshow(speed.T,cmap=cm.coolwarm)
            plt.plot(fastMarchingStartPointListFMM[1,i],\
                    fastMarchingStartPointListFMM[0,i],'ok')
            #plt.contour(speed.T,cmap=cm.coolwarm)
            plt.title('speed basin Index'+str(basinIndexList))
            plt.colorbar()
            plt.show()
            
            plt.figure(defaults.figureNumber+1)
            plt.imshow(phi.T,cmap=cm.coolwarm)
            plt.plot(fastMarchingStartPointListFMM[1,i],\
                    fastMarchingStartPointListFMM[0,i],'ok')
            #plt.contour(speed.T,cmap=cm.coolwarm)
            plt.title('phi basin Index'+str(basinIndexList))
            plt.colorbar()
            plt.show()
            
            print "I/O error({0}): {1}".format(e.errno, e.strerror)
            #stop
        
        #print travelTimearray.shape
        geodesicDistanceArray[maskedBasin ==1]= travelTimearray[maskedBasin ==1]

    #-----------------------------------
    #-----------------------------------
    # Plot the geodesic array
    defaults.figureNumber = defaults.figureNumber + 1
    plt.figure(defaults.figureNumber)
    plt.imshow(np.log10(geodesicDistanceArray.T),cmap=cm.coolwarm)
    plt.contour(geodesicDistanceArray.T,140,cmap=cm.coolwarm)
    plt.title('Geodesic distance array (travel time)')
    plt.colorbar()
    if defaults.doPlot==1:
        plt.show()
    
    print geodesicDistanceArray.shape
    # Writing the geodesic distance array
    outfilepath = Parameters.geonetResultsDir
    outfilename = Parameters.demFileName
    outfilename = outfilename.split('.')[0]+'_geodesicDistance.tif'
    write_geotif_generic(geodesicDistanceArray.T,outfilepath,outfilename)
    
    # Locating end points
    print 'Locating skeleton end points'
    xySkeletonSize = skeletonFromFlowAndCurvatureArray.shape
    skeletonLabeledArray, skeletonNumConnectedComponentsList =\
                          ndimage.label(skeletonFromFlowAndCurvatureArray)
    #print skeletonNumConnectedComponentsList
    """
     Through the histogram of skeletonNumElementsSortedList
     (skeletonNumElementsList minus the maximum value which
      corresponds to the largest connected element of the skeleton) we get the
      size of the smallest elements of the skeleton, which will likely
      correspond to small isolated convergent areas. These elements will be
      excluded from the search of end points.
    """
    print 'Counting the number of elements of each connected component'
    #print "ndimage.labeled_comprehension"
    lbls = np.arange(1, skeletonNumConnectedComponentsList+1)
    skeletonLabeledArrayNumtuple = ndimage.labeled_comprehension(skeletonFromFlowAndCurvatureArray,\
                                                                 skeletonLabeledArray,\
                                                                 lbls,np.count_nonzero,\
                                                                 int,0)
    skeletonNumElementsSortedList = np.sort(skeletonLabeledArrayNumtuple)
    print np.sqrt(len(skeletonNumElementsSortedList))
    histarray,skeletonNumElementsHistogramX=np.histogram(\
        skeletonNumElementsSortedList[0:len(skeletonNumElementsSortedList)-1],
        np.sqrt(len(skeletonNumElementsSortedList)))

    defaults.figureNumber = defaults.figureNumber + 1
    plt.figure(defaults.figureNumber)
    plt.imshow(skeletonLabeledArray.T,cmap=cm.coolwarm)
    plt.title('Skeleton Labeled Array elements Array')
    plt.colorbar()
    if defaults.doPlot==1:
        plt.show()

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
                    skeletonLabeledArrayNumtuple[skeletonLabeledArray[i,j]-1]
    
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
    skeletonNumElementsThreshold = skeletonNumElementsHistogramX[2]
    
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

    # Write shapefiles of channel heads
    write_channel_heads(xx,yy)
    
    # Do compute discrete geodesics
    print 'Computing discrete geodesics'
    compute_discrete_geodesic_v1()
    print 'Finished pyGeoNet'
    
if __name__ == '__main__':
    t0 = clock()
    main()
    t1 = clock()
    print "time taken to complete the script is::",t1-t0," seconds"
    #plt.show()
    print "script complete"
