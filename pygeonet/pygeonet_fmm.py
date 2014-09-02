import numpy as np
import pylab as pl
import skfmm
from osgeo import gdal,osr
from matplotlib import cm
from scipy import ndimage


skeltif = 'C:\\Mystuff\\grassgisdatabase\\'+'skunkroi_skeleton.tif'
dsskel = gdal.Open(skeltif, gdal.GA_ReadOnly)
aryskel = dsskel.GetRasterBand(1).ReadAsArray()
nanDemArrayskel=np.array(aryskel.T)

geoDtif = 'C:\\Mystuff\\grassgisdatabase\\'+'skunkroi_geodesicDistance.tif'
dsgeoD = gdal.Open(geoDtif, gdal.GA_ReadOnly)
arygeoD = dsgeoD.GetRasterBand(1).ReadAsArray()
nanDemArraygeoD=np.array(arygeoD.T)

pl.imshow(nanDemArrayskel,cmap=cm.binary)
pl.show()

skeletonLabeledArray, skeletonNumConnectedComponentsList =\
        ndimage.label(nanDemArrayskel)


pl.imshow(skeletonLabeledArray)
pl.colorbar()
pl.show()

skeletonNumElementsList = ndimage.measurements.histogram(skeletonLabeledArray, 1, \
                skeletonNumConnectedComponentsList, \
                bins = skeletonNumConnectedComponentsList,\
                labels=None)
print len(skeletonNumElementsList)
skeletonNumElementsSortedList = np.sort(skeletonNumElementsList)

histarray,bin_edges=np.histogram(\
    skeletonNumElementsSortedList[1:len(skeletonNumElementsSortedList)-1],
    np.sqrt(len(skeletonNumElementsSortedList)))

print bin_edges[1]

pl.bar(bin_edges[1:],histarray)
pl.show()

xySkeletonSize = skeletonLabeledArray.shape
print xySkeletonSize

skeletonNumElementsGriddedArray = np.zeros(xySkeletonSize)
#"""
for i in range(1,xySkeletonSize[0]):
    for j in range(1,xySkeletonSize[1]):
        #Gets the watershed label for this specified cell and checked in
        #subsequent if statement
        #basinIndex = subBasinIndexArray[i,j]
        if skeletonLabeledArray[i, j] > 0:
            skeletonNumElementsGriddedArray[i,j] = skeletonNumElementsList[skeletonLabeledArray[i,j]-1]


pl.imshow(skeletonNumElementsGriddedArray,cmap=cm.coolwarm)
pl.title('Skeleton Num elements Array')
pl.colorbar()
pl.show()

skeletonNumElementsThreshold = bin_edges[10]
endPointSearchBoxSize = 30
nanDemArraygeoD[np.isnan(nanDemArraygeoD)]==0
skeletonEndPointsList = []
for i in range(0,xySkeletonSize[0]):
    for j in range(0,xySkeletonSize[1]):
        #print i,j
        # Skip this pixel if the current point is not a labeled or if the
        # number of connected skeleton elements is too small
        if skeletonLabeledArray[i,j]!=0 \
            and skeletonNumElementsGriddedArray[i,j]>=skeletonNumElementsThreshold:
            # Define search box and ensure it fits within the DTM bounds
            mx = i-1
            px = xySkeletonSize[0]-i
            my = j-1
            py = xySkeletonSize[1]-j
            xMinus = np.min([endPointSearchBoxSize, mx])
            xPlus  = np.min([endPointSearchBoxSize, px])
            yMinus = np.min([endPointSearchBoxSize, my])
            yPlus  = np.min([endPointSearchBoxSize, py])
            # Extract the geodesic distances geodesicDistanceArray for pixels within the search box
            searchGeodesicDistanceBox = nanDemArraygeoD[i-xMinus:i+xPlus, j-yMinus:j+yPlus]
            #print searchGeodesicDistanceBox
            # Extract the skeleton labels for pixels within the search box
            searchLabeledSkeletonBox = skeletonLabeledArray[i-xMinus:i+xPlus, j-yMinus:j+yPlus]
            # Look in the search box for skeleton points with the same label
            # and greater geodesic distance than the current pixel at (i,j)
            # if there are none, then add the current point as a channel head
            #if i > 200 and j >200:
            #print searchGeodesicDistanceBox
            #print searchLabeledSkeletonBox
            #print skeletonLabeledArray[i,j]
            v = searchLabeledSkeletonBox==skeletonLabeledArray[i,j]
            v1 = v * searchGeodesicDistanceBox > nanDemArraygeoD[i,j]
            #v2 = nanDemArraygeoD[i,j] > v1
            v3 = np.where(np.any(v1==True,axis=0))
            """
            if i > 200 and j >200:
                print v
                print v1
                print v3
                print len(v3[0])
                stop
            #"""
            #print v3
            if len(v3[0])==0:
                #print i,j
                skeletonEndPointsList.append([i,j])
            #stop
                       
            
# For loop ends here
skeletonEndPointsListArray = np.array(skeletonEndPointsList)
print type(skeletonEndPointsListArray)
#print skeletonEndPointsList
#stop
#print len(skeletonEndPointsList[0:][0:])
xx = skeletonEndPointsListArray[0:len(skeletonEndPointsListArray),0]
yy = skeletonEndPointsListArray[0:len(skeletonEndPointsListArray),1]
pl.imshow(nanDemArrayskel,cmap=cm.binary)
pl.plot(yy,xx,'or')
pl.title('Skeleton Num elements Array with channel heads')
pl.colorbar()
pl.show() 

"""
stop
skeletonNumElementsList = np.zeros((1,skeletonNumConnectedComponentsList))
for i in range(1,skeletonNumConnectedComponentsList):
    tmp = np.where(skeletonLabeledArray == i)
    skeletonNumElementsList[0][i] = len(tmp)
skeletonNumElementsSortedList = np.sort(skeletonNumElementsList)
a = skeletonNumElementsSortedList[1:(skeletonNumElementsSortedList.size)-1]
bins = np.sqrt(skeletonNumElementsSortedList.size)
binrange = None
normed = False
weigths = None
density = None
[skeletonNumElementsHistogramX,binEdges] = np.histogram(a,bins,binrange,\
                                        normed,weigths,density)
print [skeletonNumElementsHistogramX,binEdges]

pl.bar(index,histarray)
pl.show()



stop

#-----------------
factif = 'C:\\Mystuff\\grassgisdatabase\\'+'skunkroi_fac.tif'
dsfac = gdal.Open(factif, gdal.GA_ReadOnly)
aryfac = dsfac.GetRasterBand(1).ReadAsArray()
nanDemArrayfac=np.array(aryfac.T)

basintif = 'C:\\Mystuff\\grassgisdatabase\\'+'skunkroi_basins.tif'
dsbasins = gdal.Open(basintif, gdal.GA_ReadOnly)
arybasins = dsbasins.GetRasterBand(1).ReadAsArray()
nanDemArraybasins=np.array(arybasins.T)

curvtif = 'C:\\Mystuff\\grassgisdatabase\\'+'skunkroi_curvature.tif'
dscurv = gdal.Open(curvtif, gdal.GA_ReadOnly)
arycurv = dscurv.GetRasterBand(1).ReadAsArray()
nanDemArraycurv=np.array(arycurv.T)

costtif = 'C:\\Mystuff\\grassgisdatabase\\'+'skunkroi_costfunction.tif'
dscost = gdal.Open(costtif, gdal.GA_ReadOnly)
arycost = dscost.GetRasterBand(1).ReadAsArray()
nanDemArraycost=np.array(arycost.T)

print '# of unique basins:',np.size(np.unique(nanDemArraybasins))
# Now access each unique basin and get the
# outlets for it
basinIndexList = np.unique(nanDemArraybasins)
print 'basinIndexList:', str(basinIndexList)

geodesciDistanceArray = np.zeros((nanDemArraybasins.shape))
geodesciDistanceArray[geodesciDistanceArray==0]=np.nan
for i in range(1,basinIndexList.size):
    print 'basin index',str(i)
    maskedBasin = np.zeros((nanDemArraybasins.shape))
    maskedBasin[nanDemArraybasins==basinIndexList[i]]=1
    #pl.imshow(maskedBasin,cmap=cm.binary)
    maskedBasinFAC = np.zeros((nanDemArraybasins.shape))
    maskedBasinFAC[maskedBasin==1]=\
        nanDemArrayfac[maskedBasin==1]
    maskedBasinFAC[maskedBasinFAC==0]=np.nan
    # Get the outlet of subbasin
    maskedBasinFAC[np.isnan(maskedBasinFAC)]=-9999
    subBasinoutletindices = np.where(maskedBasinFAC==maskedBasinFAC.max())
    # outlets locations in projection of the input dataset
    outletsxx = subBasinoutletindices[0]
    outletsyy = subBasinoutletindices[1]
    # call the fast marching here
    phi = -1*np.ones((nanDemArraycost.shape))
    phi[maskedBasinFAC!=-9999] = nanDemArraycost[maskedBasinFAC!=-9999]
    phi[phi==-1]=np.nan
    phi[outletsxx,outletsyy] =-1 # define your start points
    speed = phi # this is my speed == cost function
    d = skfmm.distance(1/phi, dx=1) # the geodesic distance
    t = skfmm.travel_time(d, speed, dx=1) # travel time based on geodesic distance
    d[d<0]=np.nan
    #t[t==0]=np.nan
    tlog  = t #np.sqrt(t)
    geodesciDistanceArray[maskedBasin ==1]= tlog[maskedBasin ==1]
    pl.imshow(tlog,cmap=cm.coolwarm)
    pl.plot(outletsyy,outletsxx,'ro')
    #pl.contour(tlog,100,cmap=cm.coolwarm)
    pl.colorbar()
    pl.show()


pl.imshow(np.log10(geodesciDistanceArray),cmap=cm.coolwarm)
pl.show()

#stop

# plotting the basins index DEM ( only for testing purposes)
#defaults.figureNumber = defaults.figureNumber + 1
#plt.figure(defaults.figureNumber)
maskedBasin = np.zeros((nanDemArraybasins.shape))
maskedBasin[nanDemArraybasins==basinIndexList[4]]=1
pl.imshow(maskedBasin,cmap=cm.binary)
#pl.plot(outletPointsList[0],outletPointsList[1],'go')
# For the masked basin get the maximum accumulation are
# location and use that as an outlet for the basin.
maskedBasinFAC = np.zeros((nanDemArraybasins.shape))
maskedBasinFAC[nanDemArraybasins==basinIndexList[4]]=\
    nanDemArrayfac[nanDemArraybasins==basinIndexList[4]]
maskedBasinFAC[maskedBasinFAC==0]=np.nan
pl.imshow(np.log10(maskedBasinFAC),cmap=cm.coolwarm)
pl.xlabel('X[m]')
pl.ylabel('Y[m]')
pl.title('Masked basin with outlets')
pl.show()

# Get the outlet of subbasin
maskedBasinFAC[np.isnan(maskedBasinFAC)]=0
subBasinoutletindices = np.where(maskedBasinFAC==maskedBasinFAC.max())
print subBasinoutletindices
# outlets locations in projection of the input dataset
outletsxx = subBasinoutletindices[0]
outletsyy = subBasinoutletindices[1]

# call the fast marching here
phi = -1*np.ones((nanDemArraycost.shape))
phi[maskedBasinFAC!=0] = nanDemArraycost[maskedBasinFAC!=0]
phi[phi==-1]=np.nan
phi[outletsxx,outletsyy] =-1 # define your start points
speed = phi # this is my speed == cost function
d = skfmm.distance(1/phi, dx=1) # the geodesic distance
t = skfmm.travel_time(d, speed, dx=1) # travel time based on geodesic distance
d[d<0]=np.nan
t[t==0]=np.nan
tlog  = t #np.sqrt(t)
#pl.imshow(tlog,cmap=cm.coolwarm)
pl.plot(outletsyy,outletsxx,'ro')
pl.contour(tlog,100,cmap=cm.coolwarm)
pl.colorbar()
pl.show()
#cm.Pastel2
#cm.gray_r

"""

#-----------------
"""
X, Y = np.meshgrid(np.linspace(-1,1,200), np.linspace(-1,1,200))
print X.shape

phi = 1*np.ones_like(X)
#print phi


phi[X>-0.5] = 1
phi[np.logical_and(np.abs(Y)<0.25, X>-0.75)] = 0

pl.imshow(phi)
pl.show()

pl.contour(X, Y, phi,[0], linewidths=(3), colors='black')
pl.title('Boundary location: the zero contour of phi')
#pl.savefig('2d_phi.png')
pl.show()

d = skfmm.distance(phi, dx=1e-2)
#d = skfmm.distance(nanDemArrayfac, dx=1)
pl.title('Distance from the boundary')
pl.contour(X, Y, phi,[0], linewidths=(3), colors='black')
pl.contour(X, Y, d, 15)
pl.colorbar()
#pl.savefig('2d_phi_distance.png')
pl.show()

speed = np.ones_like(X)
speed[Y>0] = 1.5
t = skfmm.travel_time(phi, speed, dx=1e-2)

pl.title('Travel time from the boundary')
pl.contour(X, Y, phi,[0], linewidths=(3), colors='black')
pl.contour(X, Y, t, 15)
pl.colorbar()
#pl.savefig('2d_phi_travel_time.png')
pl.show()

mask = np.logical_and(abs(X)<0.1, abs(Y)<0.5)
phi  = np.ma.MaskedArray(phi, mask)
t    = skfmm.travel_time(phi, speed, dx=1e-2)
pl.title('Travel time from the boundary with an obstacle')
pl.contour(X, Y, phi, [0], linewidths=(3), colors='black')
pl.contour(X, Y, phi.mask, [0], linewidths=(3), colors='red')
pl.contour(X, Y, t, 15)
pl.colorbar()
#pl.savefig('2d_phi_travel_time_mask.png')
pl.show()
"""
# end here
