import os

# Some used demFileNames
#ikawa_roi1_nutm54_clipped
#dem_2012_mission_v1

# Prepare GeoNet parameters just prior to main code execution
currentWorkingDir = os.getcwd()
geoNetHomeDir = "C:\\Mystuff\\IO_Data\\"
demDataFilePath =  geoNetHomeDir +"data\\"
demFileName = "ikawa_roi1_nutm54_clipped.tif"
Region = demFileName.split(".")[0]
geonetResultsDir = 'C:\\Mystuff\\grassgisdatabase\\'
geonetResultsBasinDir='C:\\Mystuff\\grassgisdatabase\\basinTiffs\\'

# Write shapefile file paths
shapefilepath = "C:\\Users\\Harish\\Documents\\GitHub\\pyGeoNet\\test_results\\"
driverName = "ESRI Shapefile"

pointshapefileName = demFileName.split(".")[0]+"_channelHeads"
pointFileName = shapefilepath +pointshapefileName+".shp"


drainagelinefileName = demFileName.split(".")[0]+"_channelNetwork"
drainagelineFileName = shapefilepath +drainagelinefileName+".shp"

# Things to be changed
# PM Filtered DEM to be used in GRASS GIS for flow accumulation
pmGrassGISfileName = "C:\\Mystuff\\grassgisdatabase\\PM_filtered_grassgis.tif"

# GRASS GIS Parameters
# This is the grass database directory
gisdbdir = 'C:\\Users\\Harish\\Documents\\grassdata'

# Skfmm parameters
numBasinsElements = 6

