import numpy as np
import prepare_pygeonet_defaults as defaults
import prepare_pygeonet_inputs as parameters
import pygeonet_plot as pyg_plt
import pygeonet_rasterio as pyg_rio
import pygeonet_vectorio as pyg_vio
import pygeonet_compute_local_extremas as pyg_cle


def create_xs_bank_wse(filteredDemArray, slopeDemArray, nEndPoints, geodesicPathsCellDic):
    filteredDemArray[filteredDemArray < defaults.demNanFlag] = np.nan
    if not hasattr(parameters, 'xDemSize'):
        parameters.xDemSize = np.size(filteredDemArray, 1)
    if not hasattr(parameters, 'yDemSize'):
        parameters.yDemSize = np.size(filteredDemArray, 0)
    tempArray = np.zeros_like(filteredDemArray)
    tempArray[~np.isnan(filteredDemArray)] = 999
    TotalcrossSectionsXYArray = []
    XSIDArray = []
    XSID = 1
    leftBankCellList = []
    rightBankCellList = []
    parameters.skipPixels = 5
    parameters.crosssectiongap = 5
    parameters.resultantVector = 3
    parameters.crosssectionLength = 20
    for i in range(nEndPoints):
        centerlineX = geodesicPathsCellDic[str(i)][1]
        centerlineY = geodesicPathsCellDic[str(i)][0]
        centerline = np.asarray([centerlineY, centerlineX])
        for j in range(centerline.shape[1]):
            tempArray[centerline[0, j], centerline[1, j]] = 1
        lengthOfReach = centerline.shape[1]
        numOfCrossSections = np.around(((lengthOfReach -
                                         parameters.skipPixels -
                                         parameters.skipPixels) / parameters.crosssectiongap))
        if numOfCrossSections > 0:
            crossSectionsXYArray = [[None] * 2 for n in range(
                numOfCrossSections + 1)]
            crossSectionBankXYArray = [[None] * 8 for n in range(
                numOfCrossSections + 1)]
            crossSectionsElevationArray = [[None] for n in range(
                numOfCrossSections + 1)]
            crossSectionsSlopeArray = [[None] for n in range(
                numOfCrossSections + 1)]
            crossSectionCounter = 0
            leftBankCellList.append([[], []])
            rightBankCellList.append([[], []])
            for j in range(parameters.skipPixels, (lengthOfReach - parameters.skipPixels),
                           parameters.crosssectiongap):
                pixelAtCrosssection = centerline[:, j]
                vlx = centerline[1, j + parameters.resultantVector]
                vly = centerline[0, j + parameters.resultantVector]
                vrx = centerline[1, j - parameters.resultantVector]
                vry = centerline[0, j - parameters.resultantVector]
                magnitudeV = np.sqrt((vlx - vrx) ** 2 + (vly - vry) ** 2)
                if vlx > vrx and vly < vry:
                    sintheta = (vry - vly) / magnitudeV
                    costheta = (vrx - vlx) / magnitudeV
                elif vlx < vrx and vly > vry:
                    sintheta = (vly - vry) / magnitudeV
                    costheta = (vlx - vrx) / magnitudeV
                else:
                    sintheta = abs(vry - vly) / magnitudeV
                    costheta = abs(vrx - vlx) / magnitudeV
                orthoResultant = [sintheta, costheta]
                if sintheta < 0:
                    crossSectionX = pixelAtCrosssection[1] + (
                            np.arange(-parameters.crosssectionLength,
                                      parameters.crosssectionLength + 1) * orthoResultant[0])
                    crossSectionY = pixelAtCrosssection[0] + (
                            np.arange(-parameters.crosssectionLength,
                                      parameters.crosssectionLength + 1) * orthoResultant[1])
                else:
                    crossSectionX = pixelAtCrosssection[1] - (
                            np.arange(-parameters.crosssectionLength,
                                      parameters.crosssectionLength + 1) * orthoResultant[0])
                    crossSectionY = pixelAtCrosssection[0] + (
                            np.arange(-parameters.crosssectionLength,
                                      parameters.crosssectionLength + 1) * orthoResultant[1])
                if (crossSectionX[0] - vlx) * (vly - vry) - (crossSectionY[0] - vly) * (
                        vlx - vrx) > 0:
                    crossSectionX = crossSectionX[::-1]
                    crossSectionY = crossSectionY[::-1]
                lengthOfCrosssection = len(crossSectionX)
                crosssectionElevations = np.zeros(lengthOfCrosssection)
                crosssectionSlope = np.zeros(lengthOfCrosssection)
                crosssectionX = np.full(lengthOfCrosssection, np.nan)
                crosssectionY = np.full(lengthOfCrosssection, np.nan)
                for t in range(lengthOfCrosssection):
                    if round(crossSectionX[t]) >= 0 and round(crossSectionY[t]) >= 0 \
                            and round(crossSectionX[t]) <= parameters.xDemSize - 1 \
                            and round(crossSectionY[t]) <= parameters.yDemSize - 1 \
                            and (not np.isnan(filteredDemArray[int(round(crossSectionY[t])), int(
                        round(crossSectionX[t]))])):
                        tempArray[int(round(crossSectionY[t])), int(round(crossSectionX[t]))] = 2
                        crosssectionX[t] = int(round(crossSectionX[t]))
                        crosssectionY[t] = int(round(crossSectionY[t]))
                        crosssectionElevations[t] = filteredDemArray[
                            int(round(crossSectionY[t])), int(round(crossSectionX[t]))]
                        crosssectionSlope[t] = slopeDemArray[
                            int(round(crossSectionY[t])), int(round(crossSectionX[t]))]
                crossSectionsXYArray[crossSectionCounter][1] = crosssectionX
                crossSectionsXYArray[crossSectionCounter][0] = crosssectionY
                crossSectionYArray = crosssectionY[np.logical_and(~np.isnan(crosssectionX),
                                                                  ~np.isnan(crosssectionY))]
                crossSectionXArray = crosssectionX[np.logical_and(~np.isnan(crosssectionX),
                                                                  ~np.isnan(crosssectionY))]
                TotalcrossSectionsXYArray.append([[crossSectionYArray[0], crossSectionYArray[-1]],
                                                  [crossSectionXArray[0], crossSectionXArray[-1]]])
                XSIDArray.append(XSID)
                XSID += 1
                crossSectionsElevationArray[crossSectionCounter][0] = crosssectionElevations
                crossSectionsSlopeArray[crossSectionCounter][0] = crosssectionSlope
                # For Left side of the channel
                [leftSlopeMax, iLeftSlopeMax, _, _] = pyg_cle.find_local_extremas(
                    crosssectionSlope[0:lengthOfCrosssection / 2])
                # Find the Maximum of the left side slope
                maxLeftPeakIndex = np.argmax(leftSlopeMax)
                leftBankLocation = iLeftSlopeMax[maxLeftPeakIndex]
                leftBankElevation = crosssectionElevations[leftBankLocation]
                leftBankX = int(round(crossSectionX[leftBankLocation]))
                leftBankY = int(round(crossSectionY[leftBankLocation]))
                if leftBankX >= np.size(filteredDemArray, 1):
                    leftBankX = np.size(filteredDemArray, 1) - 1
                if leftBankY >= np.size(filteredDemArray, 0):
                    leftBankY = np.size(filteredDemArray, 0) - 1
                # Store left bank locations at position 1 and 2 in the Array
                crossSectionBankXYArray[crossSectionCounter][0] = leftBankX
                crossSectionBankXYArray[crossSectionCounter][1] = leftBankY
                crossSectionBankXYArray[crossSectionCounter][2] = filteredDemArray[
                    leftBankY, leftBankX]
                crossSectionBankXYArray[crossSectionCounter][3] = slopeDemArray[
                    leftBankY, leftBankX]
                tempArray[leftBankY, leftBankX] = 3
                leftBankCellList[-1][0].append(leftBankY)
                leftBankCellList[-1][1].append(leftBankX)
                # For Right side of the channel
                [rightSlopeMax, iRightSlope, _, _] = pyg_cle.find_local_extremas(
                    crosssectionSlope[(lengthOfCrosssection - 1) / 2 + 2:])
                # Find the Maximum of the right side of the slope
                maxRightPeakIndex = np.argmax(rightSlopeMax)
                rightBankLocation = iRightSlope[maxRightPeakIndex]
                rightBankLocation = rightBankLocation + ((lengthOfCrosssection - 1) / 2) + 1
                rightBankElevation = crosssectionElevations[rightBankLocation]
                rightBankX = int(round(crossSectionX[rightBankLocation]))
                rightBankY = int(round(crossSectionY[rightBankLocation]))
                if rightBankX >= np.size(filteredDemArray, 1):
                    rightBankX = np.size(filteredDemArray, 1) - 1
                if rightBankY >= np.size(filteredDemArray, 0):
                    rightBankY = np.size(filteredDemArray, 0) - 1
                # Store right bank locations at position 3 and 4 in the Array
                crossSectionBankXYArray[crossSectionCounter][4] = rightBankX
                crossSectionBankXYArray[crossSectionCounter][5] = rightBankY
                crossSectionBankXYArray[crossSectionCounter][6] = filteredDemArray[
                    rightBankY, rightBankX]
                crossSectionBankXYArray[crossSectionCounter][7] = slopeDemArray[
                    rightBankY, rightBankX]
                tempArray[rightBankY, rightBankX] = 3
                rightBankCellList[-1][0].append(rightBankY)
                rightBankCellList[-1][1].append(rightBankX)
                # The water surface in the channel is assumed to be the minimum
                # of the maximum peaks on both sides
                virtualWaterSurfaceElevation = min(leftBankElevation, rightBankElevation)
                # virtualWaterSurfaceIndex
                crossSectionCounter += 1

    # plotting the cross sections
    if defaults.doPlot == 1:
        pyg_plt.raster_plot(tempArray, 'Reach outline and cross sections')

    # Writing the skeletonFromFlowAndCurvatureArray array
    outfilepath = parameters.geonetResultsDir
    demName = parameters.demFileName.split('.')[0]
    outfilename = demName + '_crosssection.tif'
    pyg_rio.write_geotif_generic(tempArray, outfilepath, outfilename)
    pyg_vio.write_cross_sections(TotalcrossSectionsXYArray, XSIDArray)
    pyg_vio.write_bank_lines(leftBankCellList, rightBankCellList)
