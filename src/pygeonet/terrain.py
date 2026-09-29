from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


def as_dem(dem: NDArray[np.floating], nodata: float | None = None) -> FloatArray:
    """Return a two-dimensional float64 DEM with NoData represented by NaN."""
    array = np.asarray(dem, dtype=np.float64)
    if array.ndim != 2 or min(array.shape) < 3:
        raise ValueError("DEM must be a two-dimensional array with at least 3 cells per axis")
    result = array.copy()
    if nodata is not None:
        result[result == nodata] = np.nan
    result[~np.isfinite(result)] = np.nan
    if not np.isfinite(result).any():
        raise ValueError("DEM contains no finite elevation samples")
    return np.asarray(result, dtype=np.float64)


def slope(dem: NDArray[np.floating], *, dx: float = 1.0, dy: float = 1.0) -> FloatArray:
    """Rise/run slope magnitude using row spacing ``dy`` and column spacing ``dx``."""
    if not np.isfinite(dx) or not np.isfinite(dy) or dx <= 0 or dy <= 0:
        raise ValueError("pixel spacings must be positive")
    array = as_dem(dem)
    dz_dy, dz_dx = np.gradient(array, dy, dx)
    result = np.asarray(np.hypot(dz_dx, dz_dy), dtype=np.float64)
    result[~np.isfinite(array)] = np.nan
    return np.asarray(result, dtype=np.float64)


def curvature(
    dem: NDArray[np.floating],
    *,
    dx: float = 1.0,
    dy: float = 1.0,
    method: Literal["geometric", "laplacian"] = "geometric",
) -> FloatArray:
    """Legacy-form curvature with correct independent x/y pixel spacing.

    ``geometric`` is divergence of the normalized elevation gradient.
    ``laplacian`` is d²z/dx² + d²z/dy². The sign follows upstream pyGeoNet.
    """
    if not np.isfinite(dx) or not np.isfinite(dy) or dx <= 0 or dy <= 0:
        raise ValueError("pixel spacings must be positive")
    array = as_dem(dem)
    dz_dy, dz_dx = np.gradient(array, dy, dx)
    if method == "geometric":
        magnitude = np.hypot(dz_dx, dz_dy)
        unit_x = np.divide(dz_dx, magnitude, out=np.zeros_like(dz_dx), where=magnitude > 0)
        unit_y = np.divide(dz_dy, magnitude, out=np.zeros_like(dz_dy), where=magnitude > 0)
        _, dunit_x_dx = np.gradient(unit_x, dy, dx)
        dunit_y_dy, _ = np.gradient(unit_y, dy, dx)
        result = dunit_x_dx + dunit_y_dy
    elif method == "laplacian":
        _, d2z_dx2 = np.gradient(dz_dx, dy, dx)
        d2z_dy2, _ = np.gradient(dz_dy, dy, dx)
        result = d2z_dx2 + d2z_dy2
    else:
        raise ValueError(f"unknown curvature method: {method}")
    result[~np.isfinite(array)] = np.nan
    return np.asarray(result, dtype=np.float64)
