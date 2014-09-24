# np mask testing
import numpy as np
import numpy.ma as npma
import matplotlib.pyplot as plt
from matplotlib import cm

from osgeo import gdal

outlettif = 'C:\\Mystuff\\IO_Data\\data\\ikawa_roi1_nutm54_clipped.tif'
dsout = gdal.Open(outlettif, gdal.GA_ReadOnly)
aryfdrout = dsout.GetRasterBand(1).ReadAsArray()
nanDemArrayfdrout=np.array(aryfdrout)

print nanDemArrayfdrout.shape

np.savetxt('C:\\Mystuff\\grassgisdatabase\\testikawapy.txt', nanDemArrayfdrout, delimiter=',')

nanDemArrayfdrout[nanDemArrayfdrout<0]=np.nan

slice1 = nanDemArrayfdrout#np.array(nanDemArrayfdrout[0:99,0:99])
print slice1



#np.savetxt('C:\\Mystuff\\grassgisdatabase\\testikawaMatlb.txt', nanDemArrayfdrout, delimiter=',')   # X is an array

plt.figure(1)
plt.imshow(nanDemArrayfdrout,cmap=cm.coolwarm)
plt.title('DEM')
plt.show()


plt.figure(2)
plt.imshow(slice1,cmap=cm.coolwarm)
plt.title('Slice')
#plt.show()

sliceslopeX,sliceslopeY = np.gradient(slice1,1)
sliceslope = np.sqrt(sliceslopeX**2+sliceslopeY**2)
gradXArrayT = np.divide(sliceslopeX,sliceslope)
gradYArrayT = np.divide(sliceslopeY,sliceslope)

gradGradXArray,tmpy = np.gradient(gradXArrayT,1)
tmpX,gradGradYArray = np.gradient(gradYArrayT,1)

curvatureDemArray = gradGradXArray + gradGradYArray

print curvatureDemArray[75,50]

plt.figure(3)
plt.imshow(sliceslope,cmap=cm.coolwarm)
plt.colorbar()
plt.title('Slice slope')
#plt.show()

plt.figure(3)
plt.imshow(curvatureDemArray,cmap=cm.coolwarm)
plt.title('curvatureDemArray slice')
plt.show()

np.savetxt('C:\\Mystuff\\grassgisdatabase\\curvslice.txt', curvatureDemArray, delimiter=',')


print 'mean curv:',str(np.nanmean(curvatureDemArray[:]))
print 'stdev curv:',str(np.nanstd(curvatureDemArray[:]))
# X is an array



"""
outlettif = 'C:\\Mystuff\\grassgisdatabase\\ikawa_roi1_nutm54_clipped_geodesicDistance.tif'
dsout = gdal.Open(outlettif, gdal.GA_ReadOnly)
aryfdrout = dsout.GetRasterBand(1).ReadAsArray()
nanDemArrayfdroutPy=np.array(aryfdrout)

nanDemArrayfdroutPy[nanDemArrayfdroutPy<0]=np.nan

plt.figure(2)
plt.imshow(nanDemArrayfdroutPy.T,cmap=cm.coolwarm)
plt.title('Geodesic Distance Py')
plt.show()

del dsout,aryfdrout

stop
inarray = np.array([[92,90,91,91,92,90,92,90,1429,1427,1428,1428,1429,1427,1429,1427],\
                    [1423, 1423, 1424 ,1422, 1424, 1422, 1422,1424,1423, 1423, 1424, 1422, 1424, 1422, 1422 ,1424]])

print inarray.T

tfarray = np.array([inarray[0]>0] and [inarray[0] < 1429])
tfarray1 =  np.array([inarray[1] > 0] and [inarray[1] < 1423])

print tfarray
print tfarray1

tfarray2 = tfarray * tfarray1

print tfarray2[0]

tfarray =tfarray2

tfarraynp = np.array([tfarray[0],tfarray[0]])
tfcarray = np.zeros((tfarraynp.shape))
tfcarray[tfarraynp==False]=1

print tfcarray

popin = np.where(tfarraynp[0,:]==False)

print len(popin[0])

#maskedinarray = npma.masked_array(inarray,mask=tfcarray)

maskedinarray = npma.array(inarray,mask=tfcarray)

print maskedinarray,type(maskedinarray)



r1 = maskedinarray[0,:]
r2 = maskedinarray[1,:]

print r1[~r1.mask], type(r1)
print r2[~r2.mask], type(r1)

print r1
print r2

interme = nanDemArrayfdrout[r1[~r1.mask],r2[~r2.mask]]

# create the maskless array
# we would have total of 8 cells
finalary = np.zeros((1,8))
print finalary

finalary[0,popin[0]]=np.nan
j = 0
for i in xrange(0,8):
    if ~np.isnan(finalary[0,i]):
        finalary[0,i] = interme[j]
        j = j+1
    print finalary

        

print finalary
"""
# end


