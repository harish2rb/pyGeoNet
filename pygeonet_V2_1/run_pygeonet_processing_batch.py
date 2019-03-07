import os
import glob



# Prepare GeoNet parameters just prior to main code execution
currentWorkingDir = os.getcwd()
geoNetHomeDir = currentWorkingDir + "/IO"
demDataFilePath = os.path.join(geoNetHomeDir, "data")
print(demDataFilePath)
list_of_tifs = glob.glob(demDataFilePath + "/*.tif")
for demfilename in list_of_tifs:
    print(demfilename.split('/')[-1])

