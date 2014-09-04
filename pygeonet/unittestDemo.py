import os
import numpy as np
import numpy.ma as npma
from osgeo import gdal,osr,ogr

skeltif = 'C:\\Mystuff\\grassgisdatabase\\'+'skunkroi_skeleton.tif'
dsskel = gdal.Open(skeltif, gdal.GA_ReadOnly)
aryskel = dsskel.GetRasterBand(1).ReadAsArray()
nanDemArrayskel=np.array(aryskel.T)

gtf =  dsskel.GetGeoTransform()

gtif = gdal.Open( skeltif )

geotransform = dsskel.GetGeoTransform()
originX = geotransform[0]
originY = geotransform[3]

print originX,originY

srs = osr.SpatialReference()
srs.ImportFromWkt(dsskel.GetProjection())

srsLatLong = srs.CloneGeogCS()
ct = osr.CoordinateTransformation(srs,srsLatLong)
print ct.TransformPoint(originX,originY)


stop
print 'xlowerCoord %8f' %gtf[0]

daShapefile = r"C:\Mystuff\grassgisdatabase\outletPoints\outletPoints.shp"
driver = ogr.GetDriverByName('ESRI Shapefile')
dataSource = driver.Open(daShapefile, 0)
if dataSource is None:
    print 'Could not open %s' % (daShapefile)
else:
    print 'Opened %s' % (daShapefile)
    layer = dataSource.GetLayer()
    featureCount = layer.GetFeatureCount()
    print "Number of features in %s: %d" \
          % (os.path.basename(daShapefile),featureCount)
    for feature in layer:
        geom = feature.GetGeometryRef()
        print geom.GetX(),geom.GetY()





a = []
a = np.array([[ 12 , 10 , 11,11],
 [944, 944 ,945, 943]])

b = np.array([[12 , 10 , 12, 10],
 [ 945, 943, 943 ,945]])

r1 = a.tolist()[0]
r2 = a.tolist()[1]
r3 = b.tolist()[0]
r4 = b.tolist()[1]


print a
print b

c = np.array([r1+r3,r2+r4])

print c
sizeA = [1423, 944]
lg1 = [a[0,:]>0] and [a[1,:]>0] and\
      [a[0,:] <= sizeA[0]] and [a[1,:] <= sizeA[1]]

print lg1[0]

truefalseArray = np.array([lg1[0],lg1[0]])

print truefalseArray

bf = np.array(a[truefalseArray])

print bf

truefalseArrayMask = np.zeros((truefalseArray.shape))

truefalseArrayMask[truefalseArray==False]=1
print truefalseArrayMask

af = npma.masked_array(a, mask=truefalseArrayMask)

print 'af'
print af

print np.array([bf[0:len(bf)/2],bf[3:len(bf)]])


startpoint = np.array([[12],[23]])
print 'startpoint\n',startpoint

b = np.array([[21],[32]])
startpoint = np.hstack((startpoint,b))

print 'startpoint\n',startpoint

outlet = [   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
          0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
          0,    0,    0,    0,    0,    0,    0,   10,   20,   63,   73,
         93,  158,  184,  222,  226,  241,  250,  253,  282,  284,  354,
        382,  483,  503,  513,  571,  600,  731,  731,  810,  859,  875,
        896,  902,  915,  929,  945,  953,  986,  993, 1011, 1014, 1019,
       1025, 1034, 1060, 1122, 1253, 1293, 1315, 1322, 1428, 1428, 1428,
       1428, 1428, 1428, 1428, 1428, 1428, 1428, 1428, 1428, 1428, 1428,
       1428, 1428, 1428, 1428, 1428, 1428, 1428, 1428, 1428, 1428, 1428,
       1428, 1428]
print outlet
outletfloat = [float(x) for x in outlet]

print np.array(outletfloat) * 1.0
print outletfloat
