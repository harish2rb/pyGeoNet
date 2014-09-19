# np mask testing
import numpy as np
import numpy.ma as npma
import matplotlib.pyplot as plt
from matplotlib import cm

from osgeo import gdal

outlettif = 'C:\\Mystuff\\grassgisdatabase\\skunkroi_geodesicDistance.tif'
dsout = gdal.Open(outlettif, gdal.GA_ReadOnly)
aryfdrout = dsout.GetRasterBand(1).ReadAsArray()
nanDemArrayfdrout=np.array(aryfdrout)

print nanDemArrayfdrout.shape


plt.figure(1)
plt.imshow(nanDemArrayfdrout.T,cmap=cm.coolwarm)
plt.title('Geodesic Distance')
#plt.show()

del dsout,aryfdrout



inarray = np.array([[92,90,91,91,92,90,92,90],\
                    [1422,1422,1423,1421,1423,1421,1421,1423]])

print inarray

tfarray = [inarray[0]>1] and [inarray[1] < 1423]

print tfarray

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



