import os
# Prepare GeoNet parameters just prior to main code execution
currentWorkingDir = os.getcwd()
geoNetHomeDir = "C:\\Mystuff\\IO_Data\\"
demDataFilePath =  geoNetHomeDir +"data\\"
demFileName = "dem_2012_mission_v1.tif"
#ikawa_roi1_nutm54_clipped
#dem_2012_mission_v1

print demDataFilePath
print demFileName

# Write shapefile file paths
shapefilepath = "C:\\Users\\Harish\\Documents\\GitHub\\pyGeoNet\\test_results\\"
shapefileName = demFileName.split(".")[0]+"_channelHeads"
FileName = shapefilepath +shapefileName+".shp"
driverName = "ESRI Shapefile"


