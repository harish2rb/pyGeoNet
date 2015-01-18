from osgeo import gdal,osr
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import pygeonet_defaults as defaults

#plt.ion()

# plot the flow directions
geotiffmapraster = 'skunkroi'
outputFDR_filename = geotiffmapraster + '_fdr.tif'
fdrtif = 'C:\\Mystuff\\grassgisdatabase\\'+outputFDR_filename
dsfdr = gdal.Open(fdrtif, gdal.GA_ReadOnly)
aryfdr = dsfdr.GetRasterBand(1).ReadAsArray()
nanDemArrayfdr=np.array(aryfdr.T)

"""
Output drainage raster map contains drainage direction.
Provides the "aspect" for each cell measured CCW from East.
Multiplying positive values by 45 will give the direction
in degrees that the surface runoff will travel from that cell.
The value 0 (zero) indicates that the cell is a depression area
(defined by the depression input map).

Negative values indicate that surface runoff is leaving the boundaries
of the current geographic region. The absolute value of these
negative cells indicates the direction of flow.

"""
outlets = np.where((nanDemArrayfdr<0) & (nanDemArrayfdr!=-9999))
#print outlets
#print np.transpose(outlets)

factif = 'C:\\Mystuff\\grassgisdatabase\\skunkroi_fac.tif'
dsfac = gdal.Open(factif, gdal.GA_ReadOnly)
aryfac = dsfac.GetRasterBand(1).ReadAsArray()
nanDemArrayfac=np.array(aryfac.T)

nanDemArrayfac[np.isnan(nanDemArrayfac)]=np.nan
flowMean = np.mean(nanDemArrayfac[~np.isnan(nanDemArrayfac[:])])
print 'Mean upstream flow: ', flowMean

skeletonArray = np.zeros((nanDemArrayfac.shape))
print nanDemArrayfac.shape

figureNumber = 1
plt.figure(figureNumber)
plt.imshow(np.log10(nanDemArrayfac))
plt.plot(outlets[0],outlets[1],'go')
plt.xlabel('X[m]')
plt.ylabel('Y[m]')
plt.title('Flow directions DEM')
plt.show()

print 'still going on..'
mask = np.where(nanDemArrayfac> 3000)
skeletonArray[mask] = 1

figureNumber = 2
plt.figure(figureNumber)
plt.imshow(skeletonArray,cmap=cm.binary)
plt.plot(outlets[0],outlets[1],'go')
plt.xlabel('X[m]')
plt.ylabel('Y[m]')
plt.title('skeleton from flow')
plt.show()

#plt.ioff()

if hasattr(defaults, 'reciprocalLocalCostFn'):
    print 'pass'
else:
    print 'failed'
