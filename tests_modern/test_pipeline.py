import numpy as np

from pygeonet.config import PipelineConfig
from pygeonet.pipeline import process_dem


def test_end_to_end_synthetic_valley() -> None:
    y, x = np.mgrid[-1:1:65j, -1:1:65j]
    dem = 1000 - 80 * y + 120 * x * x
    config = PipelineConfig(
        filter_iterations=0,
        diffusion_method="none",
        flow_threshold_cells=8,
        require_positive_curvature=False,
        condition_depressions=False,
    )
    result = process_dem(dem, config=config, dx=2, dy=3)
    assert result.accumulation.shape == dem.shape
    assert result.channel_mask.any()
    assert result.network.number_of_edges() > 0
    assert np.isfinite(result.slope).all()


def test_pipeline_deterministic() -> None:
    y, x = np.mgrid[:17, :19]
    dem = 100 - y + 0.1 * (x - 9) ** 2
    config = PipelineConfig(
        filter_iterations=1, flow_threshold_cells=3, require_positive_curvature=False
    )
    first = process_dem(dem, config=config)
    second = process_dem(dem, config=config)
    assert np.array_equal(first.receivers, second.receivers)
    assert np.allclose(first.accumulation, second.accumulation)
