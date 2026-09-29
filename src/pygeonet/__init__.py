"""Modern pyGeoNet numerical and geospatial APIs.

The default pysheds D8 pipeline is not a claim of parity with the legacy GRASS MFD
workflow. See ``SCIENTIFIC_VALIDATION.md`` before comparing scientific outputs.
"""

from .config import PipelineConfig
from .pipeline import PipelineResult, process_dem
from .terrain import curvature, slope

__all__ = ["PipelineConfig", "PipelineResult", "curvature", "process_dem", "slope"]
__version__ = "0.2.0rc2"
