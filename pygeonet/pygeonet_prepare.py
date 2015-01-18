import os

"""
Folder structure for pyGeoNet is as follows
geoNetHomeDir : defines where files will be written
e.g.
geoNetHomeDir = "C:\\Mystuff\\IO_Data\\"
        --- \\data     (input lidar files will be read from this folder)
        --- \\results  (outputs from pygeonet will be written to this folder)
        --- \\basinTiffs (intermediate GRASS GIS files will be written
                          and deleted from this location. some times these
                          files could be huge, so have enough space)

pmGrassGISfileName -- this is an important intermediate GRASS GIS file name.
# Skfmm parameters
numBasinsElements = 6

# Some used demFileNames
#ikawa_roi1_nutm54_clipped
#dem_2012_mission_v1

#PLEASE DO NOT CHANGE VARIABLES,UNLESS YOU KNOW WHAT YOU ARE DOING

"""



# Prepare GeoNet parameters just prior to main code execution
currentWorkingDir = os.getcwd()
geoNetHomeDir = "C:\\Mystuff\\IO_Data\\"
demDataFilePath =  geoNetHomeDir +"data\\"
demFileName = "skunk.tif" #"ikawa_roi1_nutm54_clipped.tif"
Region = demFileName.split(".")[0]


geonetResultsDir = geoNetHomeDir +"results\\"
geonetResultsBasinDir= geoNetHomeDir +"basinTiffs\\"

# Write shapefile file paths
shapefilepath = geoNetHomeDir +"results\\"
driverName = "ESRI Shapefile"

pointshapefileName = demFileName.split(".")[0]+"_channelHeads"
pointFileName = shapefilepath +pointshapefileName+".shp"


drainagelinefileName = demFileName.split(".")[0]+"_channelNetwork"
drainagelineFileName = shapefilepath +drainagelinefileName+".shp"

# Things to be changed
# PM Filtered DEM to be used in GRASS GIS for flow accumulation
pmGrassGISfileName = geonetResultsDir +"PM_filtered_grassgis.tif"

# GRASS GIS Parameters
# This is the grass database directory
gisdbdir = 'C:\\Users\\Harish\\Documents\\grassdata'

# Skfmm parameters
numBasinsElements = 2

