import numpy as np
import pylab as pl
import skfmm
from osgeo import gdal,osr
from matplotlib import cm
from scipy import ndimage
import numpy.ma as npma

skeltif = 'C:\\Mystuff\\grassgisdatabase\\'+'skunkroi_skeleton.tif'
dsskel = gdal.Open(skeltif, gdal.GA_ReadOnly)
aryskel = dsskel.GetRasterBand(1).ReadAsArray()
nanDemArrayskel=np.array(aryskel.T)
allbasins  = np.zeros((nanDemArrayskel.shape))

# read big basin files
for op in range(1,101):
    bigbasintif = 'C:\\Mystuff\\grassgisdatabase\\'+'skunkroi_onebasins_'+str(op)+'.tif'
    dsskel = gdal.Open(bigbasintif, gdal.GA_ReadOnly)
    aryskel = dsskel.GetRasterBand(1).ReadAsArray()
    nanDemArrayskel=np.array(aryskel.T)
    allbasins[nanDemArrayskel==1]=op
    dsskel = None
    

pl.imshow(allbasins,cmap=cm.Set1)
pl.colorbar()
pl.show()




