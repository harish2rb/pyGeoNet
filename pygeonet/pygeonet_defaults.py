## Set default parameters for GeoNet

# Setting up Geonet grass location
import shutil
import os


grassGISlocation = 'C:\\Users\\Harish\\Documents\\grassdata\\geonet'
if os.path.exists(grassGISlocation):
    print "Cleaning existing Grass location"
    shutil.rmtree(grassGISlocation)


dirpath = "C:\\Mystuff\\grassgisdatabase\\basinTiffs"
if os.path.exists(dirpath):
   print "Cleaning old basinTiffs"
   shutil.rmtree(dirpath)

print "Making basinTiffs"
os.mkdir(dirpath)

# Reporting, plotting and file handling
doFileOutput=1
doReport=1
doPlot=0
doResetFiguresOnStart=1


# **** Default Parameters for Perona-Malik nonlinear diffusion
# ... could be:  PeronaMalik1, PeronaMalik2, Tukey, rampPreserving
diffusionMethod='PeronaMalik2'
diffusionTimeIncrement=0.1 # this makes the explicit scheme stable
diffusionSigmaSquared=0.05
nFilterIterations=50 # Nonlinear filtering iterations

# Flow routing options and sub basin indexing
yxPixelSizeRatio=1.0
doFloodingFlag=1
doRandomizeFlowDirections=1
flowRandomizeJitterAmplitude=0.15
doSubBasinFastMarchingFlag=0
plotOrderThreshold=7
thresholdAreaSubBasinIndexing=1000

# Define the cost function
# areaArray=D8 accumulation area
# flowArray=Dinf accumulation area
# slopeDemArray=local (scalar, modulus) slope
# curvatureDemArray=geometric curvature
# areaMean=mean D8 accumulation area
# flowMean=mean Dinf accumulation area
# skeletonFromFlowArray=skeleton based on Dinf flow
# skeletonFromCurvatureArray=skeleton based on curvature
# skeletonFromFlowAndCurvatureArray=skeleton based on Dinf flow and curvature
# reciprocalLocalCostFn='(flowArray.^(1/3.0)).*(slopeDemArray.^(-1/3.0))+20'
# reciprocalLocalCostFn = 'flowArray + 30*skeletonFromFlowAndCurvatureArray + exp(curvatureDemArray*3)'
# doNormalizeCurvature=0
reciprocalLocalCostFn='flowArray + flowMean*skeletonFromFlowAndCurvatureArray + flowMean*curvatureDemArray'
doNormalizeCurvature=1
reciprocalLocalCostMinimum='nan'

# What proportion of the DEM should we track drainage?
thresholdPercentAreaForDelineation=0.0
demNanFlag=-9999#3.402823e+038
demErrorFlag=-9999#3.402823e+038

# The demSmoothingQuantile is the quantile of landscape we want to smooth and
# (1-demSmoothingQuantile) is the quantile of landscape we want to enhance.
# A good range of demSmoothingQuantile is 0.5 to 0.9
demSmoothingQuantile=0.9
curvatureCalcMethod='geometric'
thresholdQqCurvature=1
flowThresholdForSkeleton=50

endPointSearchBoxSize = 30
doTrueGradientDescent = 1
