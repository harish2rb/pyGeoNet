import numpy as np
import pylab as pl
import skfmm
from osgeo import gdal,osr
from matplotlib import cm
from scipy import ndimage
import numpy.ma as npma

# read big basin files
for op in range(1,101):
    bigbasintif = 'C:\\Mystuff\\grassgisdatabase\\'+'skunkroi_onebasins_'+str(op)+'.tif'
    dsskel = gdal.Open(bigbasintif, gdal.GA_ReadOnly)
    aryskel = dsskel.GetRasterBand(1).ReadAsArray()
    nanDemArrayskel=np.array(aryskel.T)

    pl.imshow(nanDemArrayskel)
    pl.show()




