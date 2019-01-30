# PyGeoNet functions for raster I/O
import os
import sys
import numpy as np
from osgeo import gdal
from osgeo import osr
import prepare_pygeonet_inputs as parameters


# Read dem information
def read_dem_from_geotiff(demFileName, demFilePath):
    """
    Read DEM from the input uri of the geotiff

    :param demFileName: The dem file name
    :param demFilePath: The dem full file path
    :return: The dem as a numpy array
    """
    # Open the GeoTIFF format DEM
    fullFilePath = os.path.join(demFilePath, demFileName)
    print 'reading geotiff', demFileName
    # Use GDAL functions to read the dem as a numpy array
    # and get the dem extent, resolution, and projection
    ary = []
    gdal.UseExceptions()
    ds = gdal.Open(fullFilePath, gdal.GA_ReadOnly)
    driver = ds.GetDriver()
    geotransform = ds.GetGeoTransform()
    parameters.geotransform = geotransform
    ary = ds.GetRasterBand(1).ReadAsArray()
    parameters.demPixelScale = float(geotransform[1])
    parameters.xLowerLeftCoord = float(geotransform[0])
    parameters.yLowerLeftCoord = float(geotransform[3])
    parameters.inputwktInfo = ds.GetProjection()
    # return the dem as a numpy array
    return ary


# Read geotif from file on a disk
def read_geotif_filteredDEM():
    """
    Read the Perona Malik filtered DEM

    :return: The filtered dem as a numpy array
    """
    intif = parameters.pmGrassGISfileName
    ds = gdal.Open(intif, gdal.GA_ReadOnly)
    driver = ds.GetDriver()
    ary = ds.GetRasterBand(1).ReadAsArray()
    geotransform = ds.GetGeoTransform()
    parameters.geotransform = geotransform
    parameters.demPixelScale = float(geotransform[1])
    parameters.inputwktInfo = ds.GetProjection()
    return ary


# Read geotif from file on a disk
def read_geotif_generic(intifpath, intifname):
    """
    Read a generic geotiff

    :param intifpath: The input dem file full path
    :param intifname: The input dem file name
    :return: The dem as a numpy array
    """
    intif = os.path.join(intifpath, intifname)
    ds = gdal.Open(intif, gdal.GA_ReadOnly)
    ary = ds.GetRasterBand(1).ReadAsArray()
    return ary


# Write geotif to file on a disk
def write_geotif_generic(inputArray, outfilepath, outfilename):
    """
    Write a geotiff to disk

    :param inputArray: The input array to be saved as geotif
    :param outfilepath: The output file path
    :param outfilename: The output file name
    :return: None
    """
    print 'writing geotiff', outfilename
    output_fileName = os.path.join(outfilepath, outfilename)
    # Get shape
    nrows = inputArray.shape[0]
    ncols = inputArray.shape[1]
    # create the output image
    driver = gdal.GetDriverByName('GTiff')
    outDs = driver.Create(output_fileName, ncols, nrows, 1, gdal.GDT_Float32)
    if outDs is None:
        print 'Could not create ' + outfilename
        sys.exit(1)
    outBand = outDs.GetRasterBand(1)
    # set the reference info
    geotransform = parameters.geotransform
    cc = (geotransform[0], geotransform[1], geotransform[2],
          geotransform[3], geotransform[4], geotransform[5])
    outDs.SetGeoTransform(cc)
    outDs.SetProjection(parameters.inputwktInfo)
    # write the band
    tmparray = np.array(inputArray)
    outBand.WriteArray(tmparray)
    # flush data to disk, set the NoData value and calculate stats
    outBand.FlushCache()
    del tmparray, outDs, outBand, driver


def write_geotif_filteredDEM(filteredDemArray, filepath, filename):
    """
    Write filtered geotiff to disk to be used by GRASS GIS

    :param filteredDemArray: The filtered DEM array
    :param filepath: The file path
    :param filename: The file name
    :return: None
    """
    print 'writing filtered DEM'
    output_fileName = parameters.pmGrassGISfileName
    # Create gtif
    nrows = filteredDemArray.shape[0]
    ncols = filteredDemArray.shape[1]
    print 'filtered DEM size:', str(nrows), 'rowsx', str(ncols), 'columns'
    # create the output image
    driver = gdal.GetDriverByName('GTiff')
    outDs = driver.Create(output_fileName, ncols, nrows, 1, gdal.GDT_Float32)
    if outDs is None:
        print 'Could not create tif file'
        sys.exit(1)
    # set the reference info
    geotransform = parameters.geotransform
    outDs.SetGeoTransform(geotransform)
    outDs.SetProjection(parameters.inputwktInfo)
    # write the band
    outband = outDs.GetRasterBand(1)
    outband.WriteArray(filteredDemArray)
    outRasterSRS = osr.SpatialReference(wkt=parameters.inputwktInfo)
    authoritycode = outRasterSRS.GetAuthorityCode("PROJCS")
    outRasterSRS.ImportFromEPSG(int(authoritycode))
    outDs.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()
    # finishing the writing of filtered DEM
    del outDs, outband, driver, outRasterSRS
