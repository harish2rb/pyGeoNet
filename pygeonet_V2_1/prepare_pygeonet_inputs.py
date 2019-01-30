#! /usr/bin/env python
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
geoNetHomeDir = currentWorkingDir + "/Test"
demDataFilePath = os.path.join(geoNetHomeDir, "data")
demFileName = "skunk.tif"  # "Bath2015_A01_UTM55.tif"
# channelheadFileName = "channelhead.shp"
channelheadFileName = "Hou_weights.tif"
channeljunctionFileName = "junction.shp"

geonetResultsDir = os.path.join(geoNetHomeDir, "results", demFileName[:-4])
geonetResultsBasinDir = os.path.join(geoNetHomeDir, "basinTiffs")

# Write shapefile file paths
shapefilepath = os.path.join(geoNetHomeDir, "results", demFileName[:-4])
driverName = "ESRI Shapefile"

pointshapefileName = demFileName[:-4] + "_channelHeads"
pointFileName = os.path.join(shapefilepath, pointshapefileName + ".shp")

drainagelinefileName = demFileName[:-4] + "_channelNetwork"
drainagelineFileName = os.path.join(shapefilepath, drainagelinefileName + ".shp")

junctionshapefileName = demFileName[:-4] + "_channelJunctions"
junctionFileName = os.path.join(shapefilepath, junctionshapefileName + ".shp")

streamcellFileName = os.path.join(geonetResultsDir,
                                  demFileName[:-4] + "_streamcell.csv")

xsshapefileName = demFileName[:-4] + "_crossSections"
xsFileName = os.path.join(shapefilepath, xsshapefileName + ".shp")

banklinefileName = demFileName[:-4] + "_bankLines"
banklineFileName = os.path.join(shapefilepath, banklinefileName + ".shp")

# Things to be changed
# PM Filtered DEM to be used in GRASS GIS for flow accumulation
pmGrassGISfileName = os.path.join(geonetResultsDir, "PM_filtered_grassgis.tif")

# Skfmm parameters
numBasinsElements = 6

print("Making results directory")
if not os.path.exists(geonetResultsDir):
    os.mkdir(geonetResultsDir)
