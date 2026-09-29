import unittest
import os
import numpy as np
import pygeonet_rasterio as pyg_rio
import pygeonet_nonlinear_filter as pyg_nlf


currentWorkingDir = os.getcwd()
geoNetHomeDir = currentWorkingDir + "/IO"
demDataFilePath = os.path.join(geoNetHomeDir, "data")
testdemFileName = 'skunk_test1.tif'
nanDemArray = pyg_rio.read_dem_from_geotiff(demFileName=testdemFileName,
                                            demFilePath=demDataFilePath)
nanDemArray[nanDemArray < -9999] = np.nan


class RunPyGeoNet(unittest.TestCase):
    def test_nonlinear_filtering(self):
        import prepare_pygeonet_defaults as program_defaults
        import prepare_pygeonet_inputs as program_inputs

        # Calculate pixel length scale and assume square
        program_inputs.maxLowerLeftCoord = np.max([np.size(nanDemArray, 1),
                                                   np.size(nanDemArray, 0)])

        # Compute slope magnitude for raw DEM and the threshold lambda used
        # in Perona-Malik nonlinear filtering. The value of lambda (=edgeThresholdValue)
        # is given by the 90th quantile of the absolute value of the gradient.
        program_inputs.demFileName = testdemFileName[:-4]
        edgeThresholdValue = pyg_nlf.lambda_nonlinear_filter(nanDemArray,
                                                             defaults=program_defaults,
                                                             parameters=program_inputs)
        self.assertAlmostEquals(edgeThresholdValue, 0.736466879844666)

    def test_run_pygeonet(self):
        import prepare_pygeonet_defaults as program_defaults
        import prepare_pygeonet_inputs as program_inputs
        program_defaults.doBatchProcessing = 1
        import run_pygeonet_processing as run_pyg
        run_pyg


if __name__ == '__main__':
    unittest.main()
