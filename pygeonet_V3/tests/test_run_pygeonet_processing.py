import unittest
import os
import numpy as np
import pygeonet_rasterio as pyg_rio
import pygeonet_nonlinear_filter as pyg_nlf


class RunPyGeoNet(unittest.TestCase):
    def test_nonlinear_filtering(self):
        import prepare_pygeonet_inputs as program_inputs

        program_inputs.currentWorkingDir = os.getcwd()
        program_inputs.geoNetHomeDir = program_inputs.currentWorkingDir + "/IO"
        program_inputs.demDataFilePath = os.path.join(program_inputs.geoNetHomeDir, "data")
        program_inputs.demFileName = 'skunk_test1.tif'
        nanDemArray = pyg_rio.read_dem_from_geotiff(demFileName=program_inputs.demFileName,
                                                    demFilePath=program_inputs.demDataFilePath)
        nanDemArray[nanDemArray < -9999] = np.nan

        # Calculate pixel length scale and assume square
        program_inputs.maxLowerLeftCoord = np.max([np.size(nanDemArray, 1),
                                                   np.size(nanDemArray, 0)])

        # Compute slope magnitude for raw DEM and the threshold lambda used
        # in Perona-Malik nonlinear filtering. The value of lambda (=edgeThresholdValue)
        # is given by the 90th quantile of the absolute value of the gradient.
        program_inputs.demFileName = program_inputs.demFileName[:-4]
        edgeThresholdValue = pyg_nlf.lambda_nonlinear_filter(nanDemArray)
        self.assertAlmostEqual(edgeThresholdValue, 0.736466879844666)


if __name__ == '__main__':
    unittest.main()
