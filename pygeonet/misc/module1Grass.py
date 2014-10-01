#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Harish
#
# Created:     21/08/2014
# Copyright:   (c) Harish 2014
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import sys, os
sys.path.append(os.path.join(os.environ['GISBASE'], "etc", "python"))

import grass.script as g
import grass.script.setup as gsetup

import grass.script.array as garray

from osgeo import gdal,osr
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

gisbase = os.environ['GISBASE']
# lets make a new location
gisdbdir = 'C:\\Users\\Harish\\Documents\\grassdata'
locationGeonet = 'geonet'
mapsetGeonet = 'geonetuser'
geotiff = 'C:\\Users\\Harish\\Documents\\GitHub\\pyGeoNet\\data\\ikawa_lidardem_z54_masked.tif'
print geotiff
ds = gdal.Open(geotiff, gdal.GA_ReadOnly)
geotransform = ds.GetGeoTransform()
prj_wkt = ds.GetProjection()
srs=osr.SpatialReference(wkt=prj_wkt)
print prj_wkt

if srs.IsProjected:
    print srs.GetAttrValue('EPSG')
    print srs.GetAttrValue('geogcs')

#print ds.GetSpatialRef()

geotiffPM = 'C:\\Mystuff\\grassgisdatabase\\PM_filtered_grassgis.tif'
print geotiffPM
dsPM = gdal.Open(geotiffPM, gdal.GA_ReadOnly)
print dsPM.GetProjection()
print dsPM.GetMetadata()

stop

#stop
print 'Making the geonet location'
print g.run_command('g.proj', georef=geotiff,\
                    location= locationGeonet)

#gisdbdir = 'C:\Users\Harish\Documents\grassdata'
location = locationGeonet  #'australia'
mapset = mapsetGeonet  #'missioncreek'
print 'Mapsets after making locations'
print g.read_command('g.mapsets', flags = 'l')

print 'Setting environ'
print gsetup.init(gisbase, gisdbdir, locationGeonet, 'PERMANENT')
print g.gisenv()

#create new mapset in the location
print 'Making mapset now'
#print g.run_command('g.mapset', flags = 'c', mapset = 'water', location = 'australia', gisdb = 'C:\Users\Harish\Documents\grassdata')
print g.run_command('g.mapset', flags = 'c', mapset = mapsetGeonet,\
                    location = locationGeonet, gisdb = gisdbdir)

#after adding new mapset
print 'mapsets after making new'
print g.read_command('g.mapsets', flags = 'l')

print gsetup.init(gisbase, gisdbdir, locationGeonet, mapsetGeonet)
print 'r.in.gdal'
#outfiletiff = 'C:\\Mystuff\\grassgisdatabase\\PM_rotated.tif'
#os.sys('gdalwarp '+geotiff+' '+outfiletiff+' -t_srs "+proj=longlat +ellps=WGS84"')
print g.run_command('r.in.gdal', input=geotiff, output='dem_2012_mission_v1', \
    overwrite=True)
#ds = gdal.Open(geotiff, gdal.GA_ReadOnly)
#geotransform = ds.GetGeoTransform()
#DEM_arr = 'dem_2012_mission_v1@water'#ds.GetRasterBand(1).ReadAsArray()
###****************************region adjustment***********************************
### We create a temporary region that is only valid in this python session
#g.use_temp_region()
#rows = 1423#7739#DEM_arr.shape[0]
#cols = 1429#8356#DEM_arr.shape[1]
#resolution = 1
#print rows,cols
#n = 4399029.000000002#5175067.02134887 #some arbitrary value
#s = 4397606.000000002#5167328.02134887
#e = 445876.99999999907#557587.86925264  #some arbitrary value
#w = 444447.99999999907#549231.86925264
#print "g.region command"
#print g.run_command('g.region', flags = 'p', \
#              n = n ,s = s, e = e, w = w,\
#              res = resolution, rows = rows ,cols = cols)

#Flow computation for massive grids (float version)
print "r.watershed command"
print g.run_command('r.watershed',overwrite=True,\
        elevation='dem_2012_mission_v1',\
        threshold=1000, accumulation='acc_v23',basin = 'bas1v23',\
        drainage = 'dra1v23')
print g.read_command('g.list', _type = 'rast')

oneoutletbasin = 'C:\Mystuff\grassgisdatabase\basin_out23.tif' 
print "r.water.outlet"

print g.run_command('r.water.outlet',overwrite=True,\
                    input = 'dra1v23', output= 'oneoutletbasin',\
                    coordinates=[444470.550695,4397299.50548])



# read map
a = garray.array()
a.read(acc_v23)

print a

print g.raster_info(acc_v23)['datatype']
stop

print g.run_command('r.out.gdal',overwrite=True,\
        input = "oneoutletbasin",type='Float32', \
        output='C:\Mystuff\grassgisdatabase\onebasinProg.tif',nodata=-9999)

print g.read_command('g.list', _type = 'rast')
##stop
print g.run_command('r.out.gdal',overwrite=True,\
 input='acc_v23',type='Float64', \
 output='C:\Mystuff\grassgisdatabase\missioncreek_fac23.tif', format='GTiff')

print 'drainages'
# We create a new array and read map2d_1 to modify it
print g.run_command('r.out.gdal',overwrite=True, type='Float32',\
        input = "dra1v23", nodata=-9999,\
        output='C:\Mystuff\grassgisdatabase\dra1v23.tif')

fdrtif = 'C:\\Mystuff\\grassgisdatabase\\dra1v23.tif'
dsfdr = gdal.Open(fdrtif, gdal.GA_ReadOnly)
dsfdrarry = dsfdr.GetRasterBand(1).ReadAsArray()
nanDemArrayfdr=np.array(dsfdrarry)
outlets = np.where((nanDemArrayfdr<0) & (nanDemArrayfdr!=-9999))
print outlets

factif = 'C:\\Mystuff\\grassgisdatabase\\missioncreek_fac23.tif'
ds = gdal.Open(factif, gdal.GA_ReadOnly)
geotransform = ds.GetGeoTransform()
ary = ds.GetRasterBand(1).ReadAsArray()

nanDemArray=np.array(ary.T)
#nanDemArray[nanDemArray < -999]= np.nan

figureNumber = 1
plt.figure(figureNumber)
plt.imshow(nanDemArray,cmap=cm.coolwarm)
plt.plot(outlets[0],outlets[1],'ro')
plt.xlabel('X[m]')
plt.ylabel('Y[m]')
plt.title('bigbasin DEM')
plt.show()

del ds,ary
