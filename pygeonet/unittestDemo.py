import os
import numpy as np
import numpy.ma as npma

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

