import warnings
import os
import glob
from time import clock
import numpy as np
from matplotlib import cm
import prepare_pygeonet_defaults as program_defaults
import prepare_pygeonet_inputs as program_inputs
import pygeonet_rasterio as pyg_rio
import pygeonet_plot as pyg_plt
import pygeonet_nonlinear_filter as pyg_nlf
import pygeonet_slope_curvature as pyg_slc
import pygeonet_flow_accumulation as pyg_fac
import pygeonet_skeleton_definition as pyg_skd
import pygeonet_fast_marching as pyg_fm
import pygeonet_channel_head_definition as pyg_chd
import pygeonet_network_delineation as pyg_ntd
import pygeonet_xsbank_extraction as pyg_xsb


def main():
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    print("current working directory : {}".format(os.getcwd()))
    print("Reading input file path : {}".format(parameters.demDataFilePath))
    print("Reading input file : {}".format(parameters.demFileName))
    nanDemArray = pyg_rio.read_dem_from_geotiff(parameters.demFileName,
                                                parameters.demDataFilePath)
    nanDemArray[nanDemArray < defaults.demNanFlag] = np.nan
    defaults.figureNumber = 0
    # plotting the raw DEM
    if defaults.doPlot == 1:
        pyg_plt.raster_plot(nanDemArray, 'Input DEM')

    # Area of analysis
    parameters.yDemSize = np.size(nanDemArray, 0)
    parameters.xDemSize = np.size(nanDemArray, 1)

    # Calculate pixel length scale and assume square
    parameters.maxLowerLeftCoord = np.max([parameters.xDemSize, parameters.yDemSize])
    print("Raw DEM size is {0} rows X {1}  columns".format(parameters.yDemSize,
                                                           parameters.xDemSize))

    # Compute slope magnitude for raw DEM and the threshold lambda used
    # in Perona-Malik nonlinear filtering. The value of lambda (=edgeThresholdValue)
    # is given by the 90th quantile of the absolute value of the gradient.
    edgeThresholdValue = pyg_nlf.lambda_nonlinear_filter(nanDemArray)
    # performing PM filtering using the anisodiff
    print('Performing Perona-Malik nonlinear filtering')
    t0 = clock()
    filteredDemArray = pyg_nlf.anisodiff(nanDemArray, defaults.nFilterIterations,
                                         edgeThresholdValue,
                                         defaults.diffusionTimeIncrement,
                                         (parameters.demPixelScale,
                                          parameters.demPixelScale), 2)
    t1 = clock()
    print("time taken to complete nonlinear filtering: {} seconds".format(t1 - t0))
    # plotting the filtered DEM
    if defaults.doPlot == 1:
        pyg_plt.raster_plot(filteredDemArray, 'Filtered DEM')
    # Writing the filtered DEM as a tif
    pyg_rio.write_geotif_filteredDEM(filteredDemArray, parameters.demDataFilePath,
                                     parameters.demFileName)
    # Computing slope of filtered DEM
    print('Computing slope of filtered DEM')
    slopeDemArray = pyg_slc.compute_dem_slope(filteredDemArray, parameters.demPixelScale)
    # Writing the curvature array
    outfilepath = parameters.geonetResultsDir
    demName = parameters.demFileName.split('.')[0]
    outfilename = demName + '_slope.tif'
    pyg_rio.write_geotif_generic(slopeDemArray, outfilepath, outfilename)
    # Computing curvature
    print('computing curvature of filtered DEM')
    curvatureDemArray, curvatureDemMean, \
    curvatureDemStdDevn = pyg_slc.compute_dem_curvature(filteredDemArray,
                                                        parameters.demPixelScale,
                                                        defaults.curvatureCalcMethod)
    # Writing the curvature array
    outfilename = demName + '_curvature.tif'
    pyg_rio.write_geotif_generic(curvatureDemArray, outfilepath, outfilename)
    # plotting the curvature image
    if defaults.doPlot == 1:
        pyg_plt.raster_plot(curvatureDemArray, 'Curvature DEM')

    # *************************************************
    # Compute curvature quantile-quantile curve
    # This seems to take a long time ... is commented for now
    print('Computing curvature quantile-quantile curve')
    # osm,osr = compute_quantile_quantile_curve(finiteCurvatureDemList)
    # print osm[0]
    # print osr[0]
    thresholdCurvatureQQxx = 1
    # TODO have to add method to automatically compute the thresold
    # *************************************************

    # Computing contributing areas
    print('Computing upstream accumulation areas using MFD from GRASS GIS')
    """
    return {'outlets':outlets, 'fac':nanDemArrayfac ,\
            'fdr':nanDemArrayfdr ,'basins':nanDemArraybasins,\
            'outletsxxProj':outletsxxProj, 'outletsyyProj':outletsyyProj,\
            'bigbasins':allbasins}
    """
    # Call the flow accumulation function
    flowroutingresults = pyg_fac.flowaccumulation(filteredDemArray)

    # Read out the flowroutingresults into appropriate variables
    outletPointsList = flowroutingresults['outlets']
    flowDirectionsArray = flowroutingresults['fdr']
    flowArray = flowroutingresults['fac']
    flowArray[np.isnan(filteredDemArray)] = np.nan
    flowMean = np.mean(flowArray[~np.isnan(flowArray[:])])
    print('Mean upstream flow: {}'.format(flowMean))

    # These are actually not sub basins, if the basin threshold is large, then you might have as
    # nulls, so best practice is to keep the basin threshold close to 1000 default value is
    # 10,000
    basinIndexArray = flowroutingresults['bigbasins']

    # Define a skeleton based on flow alone
    skeletonFromFlowArray = \
        pyg_skd.compute_skeleton_by_single_threshold(flowArray,
                                                     defaults.flowThresholdForSkeleton)
    # Define a skeleton based on curvature alone
    skeletonFromCurvatureArray = \
        pyg_skd.compute_skeleton_by_single_threshold(curvatureDemArray,
                                                     curvatureDemMean +
                                                     thresholdCurvatureQQxx * curvatureDemStdDevn)
    # Writing the skeletonFromCurvatureArray array
    outfilename = demName + '_curvatureskeleton.tif'
    pyg_rio.write_geotif_generic(skeletonFromCurvatureArray,
                                 outfilepath, outfilename)
    # Define a skeleton based on curvature and flow
    skeletonFromFlowAndCurvatureArray = \
        pyg_skd.compute_skeleton_by_dual_threshold(curvatureDemArray,
                                                   flowArray,
                                                   curvatureDemMean +
                                                   thresholdCurvatureQQxx * curvatureDemStdDevn,
                                                   defaults.flowThresholdForSkeleton)
    # plot the skeleton with outlets
    if defaults.doPlot == 1:
        pyg_plt.raster_point_plot(skeletonFromFlowAndCurvatureArray, outletPointsList,
                                  'Skeleton with outlets', cm.binary)
    # Writing the skeletonFromFlowAndCurvatureArray array
    outfilename = demName + '_skeleton.tif'
    pyg_rio.write_geotif_generic(skeletonFromFlowAndCurvatureArray,
                                 outfilepath, outfilename)
    # Making outlets for FMM
    print(type(outletPointsList))
    print(outletPointsList)
    fastMarchingStartPointListFMM = pyg_fm.Fast_Marching_Start_Point_Identification(
        outletPointsList,
        basinIndexArray)
    # Computing the local cost function
    print('Preparing to calculate cost function')
    curvatureDemArray = pyg_fm.Curvature_Preparation(curvatureDemArray)
    # Calculate the local reciprocal cost (weight, or propagation speed in the
    # eikonal equation sense).  If the cost function isn't defined, default to
    # old cost function.
    print('Calculating local costs')
    reciprocalLocalCostArray = pyg_fm.Local_Cost_Computation(flowArray, flowMean,
                                                             skeletonFromFlowAndCurvatureArray,
                                                             curvatureDemArray)
    geodesicDistanceArray = pyg_fm.Fast_Marching(fastMarchingStartPointListFMM,
                                                 basinIndexArray, flowArray,
                                                 reciprocalLocalCostArray)
    if defaults.channelheadPredefined == 0:
        # Locating end points
        print('Defining channel heads')
        xx, yy = pyg_chd.Channel_Head_Definition(skeletonFromFlowAndCurvatureArray,
                                                 geodesicDistanceArray)
    else:
        print('Using given channel heads')
        channelhead_filename = demName + '_channelHeads.tif'
        channelheadArray = pyg_rio.read_geotif_generic(outfilepath, channelhead_filename)
        channelheadArray = np.where(channelheadArray == 1)
        xx = channelheadArray[1]
        yy = channelheadArray[0]
    geodesicPathsCellDic, numberOfEndPoints = \
        pyg_ntd.Channel_Definition(xx, yy,
                                   geodesicDistanceArray,
                                   basinIndexArray,
                                   flowDirectionsArray)
    pyg_xsb.create_xs_bank_wse(filteredDemArray, slopeDemArray, numberOfEndPoints,
                               geodesicPathsCellDic)
    print('Finished pyGeoNet')


if __name__ == '__main__':
    t0 = clock()
    main()
    t1 = clock()
    print("time taken to complete the script is:: {} seconds".format(t1 - t0))
    print("script complete")
