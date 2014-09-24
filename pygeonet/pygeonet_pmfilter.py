import numpy as np
import numpy.ma as npma
import matplotlib.pyplot as plt
from matplotlib import cm
from osgeo import gdal
from time import clock
import scipy.signal as conv2
import numpy.ma as npma

def get_boudaryindex(inputArray):
    print 'getting boundary'
    print inputArray.shape
    # we don't want to use the nans
    tmpMask = np.zeros((inputArray.shape))
    tmpMask[np.isnan(inputArray)]= 1
    mx = npma.masked_array(inputArray, mask=tmpMask)
    
    #xL= np.mean(mx[:,0:2],axis=1)

    #print xL
    
    plt.figure(1)
    plt.imshow(mx,cmap=cm.coolwarm)
    plt.title('Mask filter')
    plt.colorbar()
    #plt.show()

    difftemp = np.diff(mx)
    boundaryIndex = np.nonzero(difftemp<0)
    
    print boundaryIndex
    
    plt.figure(2)
    plt.imshow(difftemp,cmap=cm.coolwarm)
    #plt.plot(boundaryIndex[1],boundaryIndex[0],'og')
    plt.title('diff Mask filter')
    plt.colorbar()
    plt.show()
    
    stop
    
    

# pygeonet_nonlinear filtering
def simple_gaussian_smoothing(inputDemArray,kernelWidth,diffusionSigmaSquared):
    """
    smoothing input array with gaussian
    """
    [Ny,Nx]=inputDemArray.shape;
    #print Ny,Nx
    halfKernelWidth=(kernelWidth-1)/2;
    #print "halfKernelWidth",halfKernelWidth
    # Make a ramp array with 5 rows each containing [-2, -1, 0, 1, 2]
    x = np.linspace(-halfKernelWidth, halfKernelWidth, kernelWidth)
    #x= np.ones((kernelWidth,1))* range(-halfKernelWidth,halfKernelWidth+1)
    #y=x.T
    y = x
    xv,yv = np.meshgrid(x,y)
    gaussianFilter = np.exp(-(xv**2+yv**2)/(2*diffusionSigmaSquared))  # 2D Gaussian
    #print "gaussianFilter", gaussianFilter
    gaussianFilter=gaussianFilter/np.sum(gaussianFilter[:]) # Normalize
    #print inputDemArray[:,0:halfKernelWidth]
    xL= np.nanmean(inputDemArray[:,0:halfKernelWidth],axis=1)
    #xL = np.matrix(xL)
    xR= np.nanmean(inputDemArray[:,Nx-halfKernelWidth:Nx],axis =1)
    #xR = np.matrix(xR)
    #print "xR",xR.shape
    part1T = np.vstack((xL,xL))
    part1 = part1T.T
    #part1 = xL * np.matrix(np.ones((1,halfKernelWidth)))
    part2T = np.vstack((xR,xR))
    part2 = part2T.T
    #part2 = xR * np.matrix(np.ones((1,halfKernelWidth)))
    #print part1.shape, part2.shape
    eI = np.hstack((part1,inputDemArray,part2))
    #eI= np.matrix(np.concatenate((part1,inputDemArray,part2),1))
    #print 'eI',eI.shape
    xU= np.nanmean(eI[0:halfKernelWidth,:],axis=0)
    #xU = np.matrix(xU)
    #print "xU",xU.shape
    xD= np.nanmean(eI[Ny-halfKernelWidth:Ny,:],axis =0)
    #xD = np.matrix(xD)
    #print "xD",xD.shape
    part3 = np.vstack((xU,xU))
    part4 = np.vstack((xD,xD))
    #print part3.shape, part4.shape
    #part3 = np.matrix(np.ones((halfKernelWidth,1)))*xU
    #part4 = np.matrix(np.ones((halfKernelWidth,1)))*xD
    # Generate the expanded DTM array, 4 pixels wider in both x,y directions
    eI = np.vstack((part3,eI,part4))
    #eI= np.matrix(np.concatenate((part3, eI, part4),0))
    #print 'eI',eI.shape
    # The 'valid' option forces the 2d convolution to clip 2 pixels off the edges
    # NaNs spread from one pixel to a 5x5 set centered on the NaN
    #smoothedDemArray=conv2.convolve2d(eI,gaussianFilter,'valid'); # original
    fillvalue = np.nanmean(inputDemArray[:])
    smoothedDemArray=conv2.convolve2d(eI,gaussianFilter,'valid')
    return smoothedDemArray


def geonet_diffusion(demArray, diffusionMethod, nFilterIterations,\
    edgeThreshold, diffusionTimeIncrement, diffusionSigmaSquared, pixelSize):
    """
	References:
	   Based on diffusion() by Guy Gilboa
	   Code imported from GeoNet2.1 by Harish Sangireddy June 2014
    """
    print 'Performing Perona Malik'
    print diffusionMethod, nFilterIterations,edgeThreshold, \
    diffusionTimeIncrement, diffusionSigmaSquared, pixelSize
    # DTM dimensions
    [Ny,Nx]=demArray.shape;
    print Ny,Nx
    for i in range(0,nFilterIterations):
        print "iteration",i
        # Gaussian filter the DTM using a 5x5 kernel (Catte et al)
        if diffusionSigmaSquared>0:
            originalDemArray = demArray   # Save original DTM array
            demArrayout = simple_gaussian_smoothing(demArray,5,diffusionSigmaSquared)
        del demArray
        demArray = demArrayout
        #print 'demArray after gaussian smoothing',demArray.shape
        # Now calculate gradient in all directions (N,S,E,W) by simple differencing
        # - with repeat padding in each direction.
        # This step will propagate NaNs one pixel inward in each dirn.
        demArrayMatrix = np.matrix(demArray)
        In = (np.concatenate((demArrayMatrix[0,:],demArrayMatrix[0:Ny-1,:]),0)\
              - demArrayMatrix)/pixelSize
        #print In.shape
        #stop
        Is = (np.concatenate((demArrayMatrix[1:Ny,:],demArrayMatrix[Ny-1,:]),0)\
              - demArrayMatrix)/pixelSize
        #print Isgr.shape, Is.shape
        #print np.nanmean(Isgr - Is)
        #stop
        Ie = (np.concatenate((demArrayMatrix[:,1:Nx],demArrayMatrix[:,Nx-1]),1)\
              - demArrayMatrix)/pixelSize
        #print Ie.shape
        
        Iw = (np.concatenate((demArrayMatrix[:,0],demArrayMatrix[:,0:Nx-1]),1)\
              - demArrayMatrix)/pixelSize
        #print Iw.shape
        In[np.isnan(In)] = 0
        Is[np.isnan(Is)] = 0
        Ie[np.isnan(Ie)] = 0
        Iw[np.isnan(Iw)] = 0
        
        # Calculate diffusion coefficients in all dirns according to diffusionMethod
        if diffusionMethod =='linear':
            Cn=edgeThreshold
            Cs=edgeThreshold
            Ce=edgeThreshold
            Cw=edgeThreshold
        elif diffusionMethod == 'PeronaMalik1':
            Cn=np.exp(-(np.abs(In)/edgeThreshold)**2)
            Cs=np.exp(-(np.abs(Is)/edgeThreshold)**2)
            Ce=np.exp(-(np.abs(Ie)/edgeThreshold)**2)
            Cw=np.exp(-(np.abs(Iw)/edgeThreshold)**2)
        elif diffusionMethod == 'PeronaMalik2':
            Cn=1/(1+ np.array(np.abs(In)/edgeThreshold)**2)
            Cs=1/(1+ np.array(np.abs(Is)/edgeThreshold)**2)
            Ce=1/(1+ np.array(np.abs(Ie)/edgeThreshold)**2)
            Cw=1/(1+ np.array(np.abs(Iw)/edgeThreshold)**2)
        else:
            print 'Unknown smoothing method', diffusionMethod

        if diffusionSigmaSquared>0:
            #Calculate real gradients (not smoothed) - with repeat padding in each
            # direction.  This step will propagate NaNs one pixel inward in each dirn.
            originalDemArrayMatrix = np.matrix(originalDemArray)
            In = (np.concatenate((originalDemArrayMatrix[0,:],originalDemArrayMatrix[0:Ny-1,:]),0)\
                  - originalDemArrayMatrix)/pixelSize
            Is = (np.concatenate((originalDemArrayMatrix[1:Ny,:],originalDemArrayMatrix[Ny-1,:]),0)\
                  - originalDemArrayMatrix)/pixelSize
            Ie = (np.concatenate((originalDemArrayMatrix[:,1:Nx],originalDemArrayMatrix[:,Nx-1]),1)\
                  - originalDemArrayMatrix)/pixelSize
            Iw = (np.concatenate((originalDemArrayMatrix[:,0],originalDemArrayMatrix[:,0:Nx-1]),1)\
                  - originalDemArrayMatrix)/pixelSize
            In[np.isnan(In)] = 0
            Is[np.isnan(Is)] = 0
            Ie[np.isnan(Ie)] = 0
            Iw[np.isnan(Iw)] = 0
            demArray=originalDemArray

        part6 = np.array(Cn)*np.array(In) + np.array(Cs)*np.array(Is) +\
                np.array(Ce)*np.array(In) + np.array(Cw)*np.array(Iw)
        demArrayMatrix = demArrayMatrix + diffusionTimeIncrement*(part6)
        #print demArrayMatrix.shape
        
        demArrayMatrix = demArray
        #print demArrayMatrix.shape
        
    return demArrayMatrix

def compute_slope(inputarray,pixelSize):
    print 'Computing slope of filtered DTM'
    print type(inputarray),len(inputarray.shape),inputarray.shape
    slopeXArray,slopeYArray = np.gradient(inputarray,pixelSize)
    np.savetxt('C:\\Mystuff\\grassgisdatabase\\PMslpikawapyX.txt', slopeXArray, delimiter=',')
    np.savetxt('C:\\Mystuff\\grassgisdatabase\\PMslpikawapyY.txt', slopeYArray, delimiter=',')
    sloeparray = np.sqrt(slopeXArray**2 + slopeYArray**2)
    del slopeXArray,slopeYArray
    return sloeparray
    
    

def main():
    outlettif = 'C:\\Mystuff\\IO_Data\\data\\ikawa_roi1_nutm54_clipped.tif'
    dsout = gdal.Open(outlettif, gdal.GA_ReadOnly)
    aryfdrout = dsout.GetRasterBand(1).ReadAsArray()
    nanDemArray=np.array(aryfdrout)
    np.savetxt('C:\\Mystuff\\grassgisdatabase\\Origikawapy.txt', nanDemArray, delimiter=',')
    nanDemArray[nanDemArray<0] = np.nan
    diffusionMethod='PeronaMalik2'
    nFilterIterations = 50
    edgeThreshold = 1.37034323
    diffusionTimeIncrement = 0.1
    diffusionSigmaSquared = 0.05
    pixelSize = 1
    filteredDem = geonet_diffusion(nanDemArray, diffusionMethod, nFilterIterations,\
    edgeThreshold, diffusionTimeIncrement, diffusionSigmaSquared, pixelSize)
    np.savetxt('C:\\Mystuff\\grassgisdatabase\\PMikawapy.txt', filteredDem, delimiter=',')
    slopedem = compute_slope(filteredDem,pixelSize)
    np.savetxt('C:\\Mystuff\\grassgisdatabase\\PMslpikawapy.txt', slopedem, delimiter=',')


if __name__ == '__main__':
    t0 = clock()
    main()
    t1 = clock()
    print "time taken to complete the script is::",t1-t0," seconds"
    print "script complete"
