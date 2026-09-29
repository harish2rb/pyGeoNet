import numpy as np
from scipy import stats
import statsmodels.api as sm
import matplotlib.pyplot as plt

import prepare_pygeonet_defaults as defaults
import prepare_pygeonet_inputs as parameters
import pygeonet_rasterio as pyg_rio
import pygeonet_plot as pyg_plt


def compute_dem_slope(filteredDemArray, pixelDemScale):
    """
    Compute the slope of the input filtered DEM

    :param filteredDemArray: the input filtered dem array
    :param pixelDemScale: the dem resolution
    :return: the slope array
    """
    slopeXArray,slopeYArray = np.gradient(filteredDemArray, pixelDemScale)
    slopeDemArray = np.sqrt(slopeXArray**2 + slopeYArray**2)
    slopeMagnitudeDemArrayQ = slopeDemArray
    slopeMagnitudeDemArrayQ = np.reshape(slopeMagnitudeDemArrayQ,
                                         np.size(slopeMagnitudeDemArrayQ))
    slopeMagnitudeDemArrayQ = slopeMagnitudeDemArrayQ[~np.isnan(slopeMagnitudeDemArrayQ)]
    # Computation of statistics of slope
    print('slope statistics')
    print('angle min: {}'.format(np.arctan(np.percentile(slopeMagnitudeDemArrayQ,
                                                         0.1))*180/np.pi))
    print('angle max: {}'.format(np.arctan(np.percentile(slopeMagnitudeDemArrayQ,
                                                         99.9))*180/np.pi))
    print('mean slope: {}'.format(np.nanmean(slopeDemArray[:])))
    print('stdev slope: {}'.format(np.nanstd(slopeDemArray[:])))
    return slopeDemArray

    
def compute_dem_curvature(demArray, pixelDemScale, curvatureCalcMethod):
    """
    Compute the dem curvature using a specified cusrvature method

    :param demArray: the input dem array
    :param pixelDemScale: the dem resolution
    :param curvatureCalcMethod: the curvature method to be used for computing curvature
    :return: the curvature dem array, the mean and standard deviation of the curavture
    """
    gradXArray, gradYArray = np.gradient(demArray, pixelDemScale)
    slopeArrayT = np.sqrt(gradXArray**2 + gradYArray**2)
    if curvatureCalcMethod == 'geometric':
        #Geometric curvature
        print(' using geometric curvature')
        gradXArrayT = np.divide(gradXArray,slopeArrayT)
        gradYArrayT = np.divide(gradYArray,slopeArrayT)
    elif curvatureCalcMethod=='laplacian':
        # do nothing..
        print(' using laplacian curvature')
        gradXArrayT = gradXArray
        gradYArrayT = gradYArray
    
    gradGradXArray,tmpy = np.gradient(gradXArrayT,pixelDemScale)
    tmpx,gradGradYArray = np.gradient(gradYArrayT,pixelDemScale)
    curvatureDemArray = gradGradXArray + gradGradYArray
    curvatureDemArray[np.isnan(curvatureDemArray)] = 0
    del tmpy, tmpx
    # Computation of statistics of curvature
    print(' curvature statistics')
    tt = curvatureDemArray[~np.isnan(curvatureDemArray[:])]
    print(' non-nan curvature cell number: {}'.format(tt.shape[0]))
    finiteCurvatureDemList = curvatureDemArray[np.isfinite(curvatureDemArray[:])]
    print(' non-nan finite curvature cell number: {}'.format(finiteCurvatureDemList.shape[0]))
    curvatureDemMean = np.nanmean(finiteCurvatureDemList)
    curvatureDemStdDevn = np.nanstd(finiteCurvatureDemList)
    print(' mean curvature:  {}'.format(curvatureDemMean))
    print(' standard deviation: {}'.format(curvatureDemStdDevn))
    return curvatureDemArray, curvatureDemMean, curvatureDemStdDevn


def compute_quantile_quantile_curve(x):
    """
    Compute the quantile quantile plot

    :param x: the x list
    :return: the quantile quantile plot
    """
    print('getting qqplot estimate')
    if not hasattr(defaults, 'figureNumber'):
        defaults.figureNumber = 0
    defaults.figureNumber = defaults.figureNumber + 1
    plt.figure(defaults.figureNumber)
    res = stats.probplot(x, plot=plt)
    res1 = sm.ProbPlot(x, stats.t, fit=True)
    print("quantile plot {}".format(res1))
    return res


def perform_slope_curvature_computations():
    filteredDemArray = pyg_rio.read_geotif_filteredDEM()
    # Computing slope
    print('computing slope')
    slopeDemArray = compute_dem_slope(filteredDemArray, parameters.demPixelScale)
    slopeDemArray[np.isnan(filteredDemArray)] = np.nan
    # Writing the curvature array
    outfilepath = parameters.geonetResultsDir
    demName = parameters.demFileName.split('.')[0]
    outfilename = demName +'_slope.tif'
    pyg_rio.write_geotif_generic(slopeDemArray, outfilepath, outfilename)
    # Computing curvature
    print('computing curvature')
    curvatureDemArray, curvatureDemMean, \
                       curvatureDemStdDevn = compute_dem_curvature(filteredDemArray,
                                                                   parameters.demPixelScale,
                                                                   defaults.curvatureCalcMethod)
    curvatureDemArray[np.isnan(filteredDemArray)] = np.nan
    # Writing the curvature array
    outfilename = demName +'_curvature.tif'
    pyg_rio.write_geotif_generic(curvatureDemArray, outfilepath, outfilename)
    # plotting the curvature image
    if defaults.doPlot == 1:
        pyg_plt.raster_plot(curvatureDemArray, 'Curvature DEM')

    #*************************************************
    # TODO have to add method to automatically compute the threshold
    # Compute curvature quantile-quantile curve
    # This seems to take a long time ... is commented for now
    #print 'computing curvature quantile-quantile curve'
    #finiteCurvatureDemList = curvatureDemArray[np.isfinite(curvatureDemArray[:])]
    #osm,osr = compute_quantile_quantile_curve(finiteCurvatureDemList)
    #print osm[0]
    #print osr[0]
    thresholdCurvatureQQxx = 1
    # have to add method to automatically compute the thresold
    # .....
    # .....
    #*************************************************

