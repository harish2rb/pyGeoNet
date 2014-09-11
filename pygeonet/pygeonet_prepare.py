import os
# Prepare GeoNet parameters just prior to main code execution
currentWorkingDir = os.getcwd()
geoNetHomeDir = "C:\\Users\\Harish\\Dropbox\\pyGeoNet1.0"
demDataFilePath =  geoNetHomeDir +"\\data\\"
demFileName = "skunkroi.tif"

print demDataFilePath
print demFileName

# Write shapefile file paths
shapefilepath = "C:\\Users\\Harish\\Documents\\GitHub\\pyGeoNet\\test_results\\"
shapefileName = demFileName.split(".")[0]+"_channelHeads"
FileName = shapefilepath +shapefileName+".shp"
driverName = "ESRI Shapefile"
