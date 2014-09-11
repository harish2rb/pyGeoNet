# write shapefile
import osgeo.ogr as ogr
import osgeo.osr as osr
from osgeo import gdal
import os
import numpy as np
from matplotlib import cm
from scipy import ndimage
import pylab as pl

# Set the variables
shapefilepath = "C:\\Users\\Harish\\Documents\\GitHub\\pyGeoNet\\test_results\\"
shapefileName = "volcanoes"
FileName = shapefilepath +shapefileName+".shp"
driverName = "ESRI Shapefile"

# Used to get the spatial reference
skeltif = 'C:\\Mystuff\\grassgisdatabase\\'+'skunkroi_skeleton.tif'
dsskel = gdal.Open(skeltif, gdal.GA_ReadOnly)
georef = dsskel.GetProjection()
print georef
gtf = dsskel.GetGeoTransform()
aryskel = dsskel.GetRasterBand(1).ReadAsArray()
nanDemArrayskel=np.array(aryskel.T)
#--------------

# --- create the channel heads demo---
geoDtif = 'C:\\Mystuff\\grassgisdatabase\\'+'skunkroi_geodesicDistance.tif'
dsgeoD = gdal.Open(geoDtif, gdal.GA_ReadOnly)
arygeoD = dsgeoD.GetRasterBand(1).ReadAsArray()
nanDemArraygeoD=np.array(arygeoD.T)

basintif = 'C:\\Mystuff\\grassgisdatabase\\'+'skunkroi_basins.tif'
dsbasins = gdal.Open(basintif, gdal.GA_ReadOnly)
arybasins = dsbasins.GetRasterBand(1).ReadAsArray()
nanDemArraybasins=np.array(arybasins.T)

skeletonLabeledArray, skeletonNumConnectedComponentsList =\
        ndimage.label(nanDemArrayskel)

skeletonNumElementsList = ndimage.measurements.histogram(skeletonLabeledArray, 1, \
                skeletonNumConnectedComponentsList, \
                bins = skeletonNumConnectedComponentsList,\
                labels=None)
print len(skeletonNumElementsList)
skeletonNumElementsSortedList = np.sort(skeletonNumElementsList)

histarray,bin_edges=np.histogram(\
    skeletonNumElementsSortedList[1:len(skeletonNumElementsSortedList)-1],
    np.sqrt(len(skeletonNumElementsSortedList)))

print bin_edges[1]
xySkeletonSize = skeletonLabeledArray.shape
skeletonNumElementsGriddedArray = np.zeros(xySkeletonSize)
#"""
for i in range(1,xySkeletonSize[0]):
    for j in range(1,xySkeletonSize[1]):
        #Gets the watershed label for this specified cell and checked in
        #subsequent if statement
        #basinIndex = subBasinIndexArray[i,j]
        if skeletonLabeledArray[i, j] > 0:
            skeletonNumElementsGriddedArray[i,j] = skeletonNumElementsList[skeletonLabeledArray[i,j]-1]

skeletonNumElementsThreshold = bin_edges[10]
endPointSearchBoxSize = 30
nanDemArraygeoD[np.isnan(nanDemArraygeoD)]==0
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
            xMinus = np.min([endPointSearchBoxSize, mx])
            xPlus  = np.min([endPointSearchBoxSize, px])
            yMinus = np.min([endPointSearchBoxSize, my])
            yPlus  = np.min([endPointSearchBoxSize, py])
            # Extract the geodesic distances geodesicDistanceArray for pixels within the search box
            searchGeodesicDistanceBox = nanDemArraygeoD[i-xMinus:i+xPlus, j-yMinus:j+yPlus]
            #print searchGeodesicDistanceBox
            # Extract the skeleton labels for pixels within the search box
            searchLabeledSkeletonBox = skeletonLabeledArray[i-xMinus:i+xPlus, j-yMinus:j+yPlus]
            # Look in the search box for skeleton points with the same label
            # and greater geodesic distance than the current pixel at (i,j)
            # if there are none, then add the current point as a channel head
            #if i > 200 and j >200:
            #print searchGeodesicDistanceBox
            #print searchLabeledSkeletonBox
            #print skeletonLabeledArray[i,j]
            v = searchLabeledSkeletonBox==skeletonLabeledArray[i,j]
            v1 = v * searchGeodesicDistanceBox > nanDemArraygeoD[i,j]
            #v2 = nanDemArraygeoD[i,j] > v1
            v3 = np.where(np.any(v1==True,axis=0))
            """
            if i > 200 and j >200:
                print v
                print v1
                print v3
                print len(v3[0])
                stop
            #"""
            #print v3
            if len(v3[0])==0:
                #print i,j
                skeletonEndPointsList.append([i,j])
            #stop
                       
            
# For loop ends here
skeletonEndPointsListArray = np.array(skeletonEndPointsList)
print type(skeletonEndPointsListArray)

xx = skeletonEndPointsListArray[0:len(skeletonEndPointsListArray),0]
yy = skeletonEndPointsListArray[0:len(skeletonEndPointsListArray),1]

xxProj = float(gtf[0])+ \
                    float(gtf[1]) * np.array(xx) + \
                    float(0.00964)
yyProj = float(gtf[3])+ \
                    float(gtf[5])*np.array(yy) + \
                    float(0.00155)

pl.imshow(nanDemArrayskel,cmap=cm.binary)
pl.plot(yy,xx,'or')
pl.title('Skeleton Num elements Array with channel heads')
pl.colorbar()
pl.show()
#-----------Actual shape file write---------------------------
# set up the shapefile driver
driver = ogr.GetDriverByName(driverName)

# This will delete and assist in overwrite the original shape files
if os.path.exists(FileName):
     driver.DeleteDataSource(FileName)

# create the data source
data_source = driver.CreateDataSource(FileName)

# create the spatial reference, WGS84
srs = osr.SpatialReference()
#srs.ImportFromEPSG(4326)
srs.ImportFromWkt(georef)


# create the layer
layer = data_source.CreateLayer(shapefileName, srs, ogr.wkbPoint)

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
    feature.SetField("Name", 'Head')
    feature.SetField("Region", 'Skunk')
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
