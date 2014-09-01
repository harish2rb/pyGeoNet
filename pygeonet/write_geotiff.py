from osgeo import gdal,osr
import numpy as np

#dem_2012_mission_v1
fullFilePath = 'C:\\Mystuff\\grassgisdatabase\\PM_filtered.tif'

dataset = gdal.Open(fullFilePath, gdal.GA_ReadOnly)
if not dataset is None:
    print 'Driver: ', dataset.GetDriver().ShortName,'/', \
          dataset.GetDriver().LongName
    print 'Size is ',dataset.RasterXSize,'x',dataset.RasterYSize, \
          'x',dataset.RasterCount
    print 'Projection is ',dataset.GetProjection()
    geotransform = dataset.GetGeoTransform()
    if not geotransform is None:
        print geotransform
        print 'Origin = (',geotransform[0], ',',geotransform[3],')'
        print 'Pixel Size = (',geotransform[1], ',',geotransform[5],')'
    
    band = dataset.GetRasterBand(1)

    print 'Band Type=',gdal.GetDataTypeName(band.DataType)

    min = band.GetMinimum()
    max = band.GetMaximum()
    if min is None or max is None:
        (min,max) = band.ComputeRasterMinMax(1)
        print 'Min=%.3f, Max=%.3f' % (min,max)

    if band.GetOverviewCount() > 0:
        print 'Band has ', band.GetOverviewCount(), ' overviews.'

    if not band.GetRasterColorTable() is None:
        print 'Band has a color table with ', \
        band.GetRasterColorTable().GetCount(), ' entries.'

    scanline = band.ReadRaster( 0, 0, band.XSize, 1, \
                                     band.XSize, 1, gdal.GDT_Float32 )


# Write the same tif as different name
output_fileName = "C:\\Mystuff\\grassgisdatabase\\PM_filtered_t.tif"
ncols = dataset.RasterXSize
nrows = dataset.RasterYSize
print 'ncols,nrows:',ncols, nrows

# create the output image
driver = dataset.GetDriver()

#print driver
outDs = driver.Create(output_fileName, ncols, nrows, 1, gdal.GDT_Float32)
if outDs is None:
    print 'Could not create reclass_40.tif'
    sys.exit(1)
outBand = outDs.GetRasterBand(1)
outData = dataset.GetRasterBand(1).ReadAsArray()
nanDemArray=np.array(outData.T)
print nanDemArray.shape[0],nanDemArray.shape[1]

# flush data to disk, set the NoData value and calculate stats
outBand.FlushCache()
outBand.SetNoDataValue(-99)

# georeference the image and set the projection
cc = (geotransform[0],geotransform[1],geotransform[2],geotransform[3],geotransform[4],-geotransform[5])
outDs.SetGeoTransform(cc)
outDs.SetProjection(dataset.GetProjection())

# write the data
outBand.WriteArray(nanDemArray.T)

del dataset, outData,outDs,outBand

