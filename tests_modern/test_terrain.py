import numpy as np
import pytest

from pygeonet.filtering import perona_malik
from pygeonet.legacy import curvature_square_pixels
from pygeonet.terrain import as_dem, curvature, slope


def test_planar_slope_non_square_pixels() -> None:
    rows, cols = np.mgrid[:31, :41]
    dx, dy = 2.0, 5.0
    dem = 10 + 4 * (cols * dx) - 3 * (rows * dy)
    assert np.allclose(slope(dem, dx=dx, dy=dy), 5.0)


def test_laplacian_of_paraboloid() -> None:
    rows, cols = np.mgrid[-20:21, -25:26]
    dem = (cols * 2.0) ** 2 + (rows * 3.0) ** 2
    measured = curvature(dem, dx=2.0, dy=3.0, method="laplacian")
    assert np.allclose(measured[2:-2, 2:-2], 4.0)


def test_slope_rotation_invariance() -> None:
    rows, cols = np.mgrid[:21, :21]
    dem = 2 * cols + 3 * rows
    assert np.allclose(slope(dem), np.rot90(slope(np.rot90(dem)), -1))


def test_elevation_units_scale_slope_and_laplacian() -> None:
    rows, cols = np.mgrid[-10:11, -10:11]
    dem = cols**2 + rows**2
    assert np.allclose(slope(dem * 100), slope(dem) * 100)
    assert np.allclose(
        curvature(dem * 100, method="laplacian"), curvature(dem, method="laplacian") * 100
    )


def test_nodata_is_preserved_by_filter() -> None:
    dem = np.arange(49, dtype=float).reshape(7, 7)
    dem[3, 3] = np.nan
    filtered = perona_malik(dem, iterations=2, kappa=2, gamma=0.1)
    assert np.isnan(filtered[3, 3])
    assert np.isfinite(filtered[0, 0])


def test_infinite_elevations_are_normalized_to_nodata() -> None:
    dem = np.ones((3, 3), dtype=float)
    dem[1, 1] = np.inf
    assert np.isnan(as_dem(dem)[1, 1])


def test_diffusion_rejects_unstable_spacing_and_timestep() -> None:
    dem = np.arange(25, dtype=float).reshape(5, 5)
    with pytest.raises(ValueError, match="stability limit"):
        perona_malik(dem, iterations=1, kappa=1, gamma=0.15, dx=0.1, dy=0.1)


def test_legacy_curvature_characterization_matches_square_modern() -> None:
    rows, cols = np.mgrid[-10:11, -10:11]
    dem = cols**2 - rows**2 + 0.1 * cols
    assert np.allclose(curvature_square_pixels(dem, 2.0), curvature(dem, dx=2.0, dy=2.0))


def test_rejects_invalid_dem() -> None:
    with pytest.raises(ValueError):
        as_dem(np.ones((2, 2)))


def test_rejects_nonfinite_pixel_spacing() -> None:
    with pytest.raises(ValueError, match="spacings"):
        slope(np.ones((3, 3)), dx=float("nan"))
