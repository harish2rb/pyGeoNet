# np mask testing
import numpy as np
import numpy.ma as npma
import matplotlib.pyplot as plt
from matplotlib import cm

testA = np.array([[910,910,910,910,909,909,909,909,909,909,910,911,912,913,914,915,916,916\
,916,916,917,917,917,917,917,917,917,918,918,918,918,919,919,919,919,919\
,919,919,919,919,919,919,919,919,919,919,920,920,920,919,919,919,919,919\
,919,919,919,919,919,920,920,920,920,920,920,920,919,919,918,918,918,918\
,917]\
,[537,536,535,534,533,532,531,530,529,528,527,526,527,528,529,530,531,532\
,533,534,535,536,537,538,539,540,541,542,543,544,545,546,547,548,549,550\
,551,552,553,554,555,556,557,558,559,560,561,562,563,564,565,566,567,568\
,569,570,571,572,573,574,575,576,577,578,579,580,581,582,583,584,585,586\
,587]])


print testA.T


from osgeo import gdal

outlettif = 'C:\\Mystuff\\grassgisdatabase\\ikawa_roi1_nutm54_clipped_geodesic_distanceMatlb.tif'
dsout = gdal.Open(outlettif, gdal.GA_ReadOnly)
aryfdrout = dsout.GetRasterBand(1).ReadAsArray()
nanDemArrayfdrout=np.array(aryfdrout)

print nanDemArrayfdrout.shape

nanDemArrayfdrout[nanDemArrayfdrout<0]=np.nan

#np.savetxt('C:\\Mystuff\\grassgisdatabase\\testikawaMatlb.txt', nanDemArrayfdrout, delimiter=',')   # X is an array

plt.figure(1)
plt.imshow(nanDemArrayfdrout,cmap=cm.coolwarm)
plt.title('Geodesic Distance')
#plt.show()


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



