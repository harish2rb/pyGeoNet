import numpy as np
import pytest

from pygeonet.flow import basin_labels, d8_receivers, flow_accumulation, priority_flood


def test_priority_flood_fills_isolated_pit() -> None:
    dem = np.array([[5, 5, 5], [5, 0, 5], [5, 5, 5]], dtype=float)
    filled = priority_flood(dem)
    assert filled[1, 1] == 5


def test_d8_routes_planar_slope_and_conserves_area() -> None:
    rows, _ = np.mgrid[:7, :5]
    dem = -rows.astype(float) + np.zeros((7, 5))
    receiver = d8_receivers(dem)
    accumulation = flow_accumulation(receiver, valid_mask=np.ones_like(dem, dtype=bool))
    assert np.all(accumulation[-1] == 7)
    assert float(accumulation[-1].sum()) == dem.size
    labels = basin_labels(receiver, np.ones_like(dem, dtype=bool))
    assert len(np.unique(labels)) == dem.shape[1]


def test_accumulation_rejects_cycle() -> None:
    receiver = np.full((3, 3), -1, dtype=int)
    receiver.flat[0], receiver.flat[1] = 1, 0
    with pytest.raises(ValueError, match="cycle"):
        flow_accumulation(receiver)


def test_nodata_receiver_is_minus_one() -> None:
    dem = np.arange(9, dtype=float).reshape(3, 3)
    dem[1, 1] = np.nan
    assert d8_receivers(dem)[1, 1] == -1


def test_basin_labels_reject_receiver_into_nodata() -> None:
    receiver = np.full((3, 3), -1, dtype=int)
    valid = np.ones((3, 3), dtype=bool)
    valid[1, 1] = False
    receiver[1, 0] = 4
    with pytest.raises(ValueError, match="outside the valid DEM"):
        basin_labels(receiver, valid)


def test_accumulation_rejects_nonfinite_cell_area() -> None:
    with pytest.raises(ValueError, match="cell_area"):
        flow_accumulation(np.full((3, 3), -1), cell_area=float("nan"))
