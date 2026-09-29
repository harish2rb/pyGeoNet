import numpy as np
import skfmm
import prepare_pygeonet_defaults as defaults
import prepare_pygeonet_inputs as parameters
import pygeonet_rasterio as pyg_rio
import pygeonet_plot as pyg_plt


def Fast_Marching_Start_Point_Identification(outlet_array, basinIndexArray):
    """
    Identify the starting points for fast marching

    :param outlet_array: The outlet array
    :param basinIndexArray: the basin index array
    :return: The fast marching start point list
    """

    # Computing the percentage drainage areas
    print('Computing percentage drainage area of each indexed basin')
    fastMarchingStartPointList = np.array(outlet_array)
    fastMarchingStartPointListFMMx = []
    fastMarchingStartPointListFMMy = []
    basinsUsedIndexList = np.zeros((len(fastMarchingStartPointList[0]), 1))
    if not hasattr(parameters, 'xDemSize'):
        parameters.xDemSize = np.size(basinIndexArray, 1)
    if not hasattr(parameters, 'yDemSize'):
        parameters.yDemSize = np.size(basinIndexArray, 0)
    nx = parameters.xDemSize
    ny = parameters.yDemSize
    nDempixels = float(nx * ny)
    for label in range(0, len(fastMarchingStartPointList[0])):
        outletbasinIndex = basinIndexArray[fastMarchingStartPointList[0, label],
                                           fastMarchingStartPointList[1, label]]
        numelments = basinIndexArray[basinIndexArray == outletbasinIndex]
        percentBasinArea = float(len(numelments)) * 100 / nDempixels
        print("Basin {0} @ {1}  #Elements {2} and area {3} %".format(
            outletbasinIndex,
            fastMarchingStartPointList[:,
            label],
            len(numelments),
            percentBasinArea))

        if (percentBasinArea > defaults.thresholdPercentAreaForDelineation) and \
                (len(numelments) > parameters.numBasinsElements):
            # Get the watersheds used
            basinsUsedIndexList[label] = label
            # Preparing the outlets used for fast marching in ROI
            fastMarchingStartPointListFMMx.append(fastMarchingStartPointList[1, label])
            fastMarchingStartPointListFMMy.append(fastMarchingStartPointList[0, label])
        # finishing Making outlets for FMM
    # Closing Basin area computation
    fastMarchingStartPointListFMM = np.array([fastMarchingStartPointListFMMy,
                                              fastMarchingStartPointListFMMx])
    return fastMarchingStartPointListFMM


# Normalize Input Array
def normalize(inputArray):
    normalizedArray = inputArray - np.min(inputArray[~np.isnan(inputArray)])
    normalizedArrayR = normalizedArray / np.max(normalizedArray[~np.isnan(normalizedArray)])
    return normalizedArrayR


def Curvature_Preparation(curvatureDemArray):
    # lets normalize the curvature first
    if defaults.doNormalizeCurvature == 1:
        print('normalizing curvature')
        curvatureDemArray = normalize(curvatureDemArray)
        if defaults.doPlot == 1:
            pyg_plt.raster_plot(curvatureDemArray, 'Curvature DEM')
    # set all the nan's to zeros before cost function is computed
    curvatureDemArray[np.isnan(curvatureDemArray)] = 0
    return curvatureDemArray


def Local_Cost_Computation(flowArray, flowMean,
                           skeletonFromFlowAndCurvatureArray,
                           curvatureDemArray):
    if hasattr(defaults, 'reciprocalLocalCostFn'):
        print('Evaluating local cost func.')
        reciprocalLocalCostArray = eval(defaults.reciprocalLocalCostFn)
    else:
        print('Evaluating local cost func. (default)')
        reciprocalLocalCostArray = flowArray + \
                                   (flowMean * skeletonFromFlowAndCurvatureArray) \
                                   + (flowMean * curvatureDemArray)
    if hasattr(defaults, 'reciprocalLocalCostMinimum'):
        if defaults.reciprocalLocalCostMinimum != 'nan':
            reciprocalLocalCostArray[reciprocalLocalCostArray[:] \
                                     < defaults.reciprocalLocalCostMinimum] = 1.0

    print('1/cost min: {}'.format(np.nanmin(reciprocalLocalCostArray[:])))
    print('1/cost max: {}'.format(np.nanmax(reciprocalLocalCostArray[:])))

    # Writing the reciprocal array
    outfilepath = parameters.geonetResultsDir
    outfilename = parameters.demFileName
    outfilename = outfilename.split('.')[0] + '_costfunction.tif'
    pyg_rio.write_geotif_generic(reciprocalLocalCostArray, outfilepath, outfilename)
    return reciprocalLocalCostArray


def Fast_Marching(fastMarchingStartPointListFMM, basinIndexArray, flowArray,
                  reciprocalLocalCostArray):
    # Fast marching
    print('Performing fast marching')
    # Do fast marching for each sub basin
    geodesicDistanceArray = np.zeros((basinIndexArray.shape))
    geodesicDistanceArray[geodesicDistanceArray == 0] = np.Inf
    for i in range(0, len(fastMarchingStartPointListFMM[0])):
        basinIndexList = basinIndexArray[fastMarchingStartPointListFMM[0, i],
                                         fastMarchingStartPointListFMM[1, i]]
        print('basin Index: {}'.format(basinIndexList))
        print('start point : {}'.format(fastMarchingStartPointListFMM[:, i]))
        maskedBasin = np.zeros((basinIndexArray.shape))
        maskedBasin[basinIndexArray == basinIndexList] = 1
        maskedBasinFAC = np.zeros((basinIndexArray.shape))
        maskedBasinFAC[basinIndexArray == basinIndexList] = \
            flowArray[basinIndexArray == basinIndexList]
        # call the fast marching here
        phi = np.nan * np.ones((reciprocalLocalCostArray.shape))
        speed = np.ones((reciprocalLocalCostArray.shape)) * np.nan
        phi[maskedBasinFAC != 0] = 1
        speed[maskedBasinFAC != 0] = reciprocalLocalCostArray[maskedBasinFAC != 0]
        phi[fastMarchingStartPointListFMM[0, i],
            fastMarchingStartPointListFMM[1, i]] = -1
        try:
            travelTimearray = skfmm.travel_time(phi, speed, dx=1)
        except IOError as e:
            print('Error in calculating skfmm travel time')
            print('Error in catchment: {}'.format(basinIndexList))
            print("I/O error({0}): {1}".format(e.errno, e.strerror))
            # setting travel time to empty array
            travelTimearray = np.nan * np.zeros((reciprocalLocalCostArray.shape))
            if defaults.doPlot == 1:
                pyg_plt.raster_point_plot(speed, fastMarchingStartPointListFMM[:, i],
                                          'speed basin Index' + str(basinIndexList))
                # plt.contour(speed,cmap=cm.coolwarm)
                pyg_plt.raster_point_plot(phi, fastMarchingStartPointListFMM[:, i],
                                          'phi basin Index' + str(basinIndexList))
        except ValueError:
            print('Error in calculating skfmm travel time')
            print('Error in catchment: '.format(basinIndexList))
            print("Oops!  That was no valid number.  Try again...")
        geodesicDistanceArray[maskedBasin == 1] = travelTimearray[maskedBasin == 1]
    geodesicDistanceArray[geodesicDistanceArray == np.Inf] = np.nan
    # Plot the geodesic array
    if defaults.doPlot == 1:
        pyg_plt.geodesic_contour_plot(geodesicDistanceArray,
                                      'Geodesic distance array (travel time)')
    # Writing the geodesic distance array
    outfilepath = parameters.geonetResultsDir
    demName = parameters.demFileName.split('.')[0]
    outfilename = demName + '_geodesicDistance.tif'
    pyg_rio.write_geotif_generic(geodesicDistanceArray, outfilepath, outfilename)
    return geodesicDistanceArray
