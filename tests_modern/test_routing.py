import numpy as np

from pygeonet.flow import flow_accumulation
from pygeonet.routing import route_dem


def test_pysheds_is_default_and_conserves_d8_area() -> None:
    rows, _ = np.mgrid[:7, :5]
    dem = -rows.astype(float) + np.zeros((7, 5))
    result = route_dem(dem, condition_depressions=False)
    assert result.algorithm == "pysheds-d8"
    assert np.all(result.accumulation[-1] == 7)
    assert float(result.accumulation[-1].sum()) == dem.size


def test_pysheds_conditions_pit_and_resolves_flow() -> None:
    dem = np.array(
        [[9, 8, 7, 6, 5], [9, 4, 4, 4, 3], [9, 4, 0, 4, 2], [9, 4, 4, 4, 1]],
        dtype=float,
    )
    result = route_dem(dem)
    assert result.conditioned_dem[2, 2] >= 4
    assert result.receivers[2, 2] >= 0
    assert np.isfinite(result.accumulation).all()


def test_native_backend_remains_explicitly_available() -> None:
    dem = np.array([[3, 2, 1], [4, 3, 0], [5, 4, -1]], dtype=float)
    result = route_dem(dem, backend="native", condition_depressions=False)
    assert result.algorithm == "native-d8"


def test_pysheds_preserves_nodata_mask() -> None:
    dem = np.array([[5, 4, 3], [6, np.nan, 2], [7, 6, 1]], dtype=float)
    result = route_dem(dem, condition_depressions=False)
    assert result.receivers[1, 1] == -1
    assert np.isnan(result.conditioned_dem[1, 1])
    assert np.isnan(result.accumulation[1, 1])


def test_pysheds_accumulation_matches_converted_receiver_graph() -> None:
    rng = np.random.default_rng(42)
    rows, _ = np.mgrid[:11, :13]
    dem = -10.0 * rows + rng.normal(scale=0.01, size=(11, 13))
    dem[3, 4] = np.nan
    result = route_dem(dem, condition_depressions=False, dx=2, dy=3)
    expected = flow_accumulation(result.receivers, valid_mask=np.isfinite(dem), cell_area=6)
    assert np.allclose(result.accumulation, expected, equal_nan=True)
