from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.typing import NDArray

from .terrain import FloatArray, as_dem, slope


def diffusion_threshold(
    dem: NDArray[np.floating], quantile: float, *, dx: float, dy: float
) -> float:
    values = slope(dem, dx=dx, dy=dy)
    finite = values[np.isfinite(values)]
    if not 0 < quantile < 1:
        raise ValueError("quantile must lie strictly between zero and one")
    if not finite.size:
        raise ValueError("DEM has no finite slope samples from which to estimate a threshold")
    threshold = float(np.quantile(np.abs(finite), quantile))
    return max(threshold, float(np.finfo(np.float64).eps))


def perona_malik(
    dem: NDArray[np.floating],
    *,
    iterations: int,
    kappa: float,
    gamma: float = 0.15,
    dx: float = 1.0,
    dy: float = 1.0,
    method: Literal["perona-malik-1", "perona-malik-2"] = "perona-malik-2",
) -> FloatArray:
    """Vectorized four-neighbor anisotropic diffusion with fixed NoData mask."""
    if (
        iterations < 0
        or not np.isfinite(kappa)
        or kappa <= 0
        or not np.isfinite(dx)
        or not np.isfinite(dy)
        or dx <= 0
        or dy <= 0
    ):
        raise ValueError("iterations must be non-negative; kappa and spacings must be positive")
    if not np.isfinite(gamma) or not 0 < gamma <= 0.25:
        raise ValueError("gamma must be in (0, 0.25]")
    stability_limit = 0.5 / (dx**-2 + dy**-2)
    if gamma > stability_limit:
        raise ValueError(
            f"gamma={gamma} exceeds the explicit stability limit {stability_limit:g} "
            "for the supplied pixel spacings"
        )
    image = as_dem(dem)
    invalid = ~np.isfinite(image)
    work = image.copy()
    work[invalid] = 0.0
    valid = ~invalid
    for _ in range(iterations):
        south = np.zeros_like(work)
        east = np.zeros_like(work)
        south[:-1] = (work[1:] - work[:-1]) / dy
        east[:, :-1] = (work[:, 1:] - work[:, :-1]) / dx
        south_valid = np.zeros_like(valid)
        east_valid = np.zeros_like(valid)
        south_valid[:-1] = valid[:-1] & valid[1:]
        east_valid[:, :-1] = valid[:, :-1] & valid[:, 1:]
        south[~south_valid] = 0.0
        east[~east_valid] = 0.0
        if method == "perona-malik-1":
            conduct_s = np.exp(-((south / kappa) ** 2))
            conduct_e = np.exp(-((east / kappa) ** 2))
        elif method == "perona-malik-2":
            conduct_s = 1.0 / (1.0 + (south / kappa) ** 2)
            conduct_e = 1.0 / (1.0 + (east / kappa) ** 2)
        else:
            raise ValueError(f"unknown diffusion method: {method}")
        flux_s = conduct_s * south / dy
        flux_e = conduct_e * east / dx
        divergence = flux_s + flux_e
        divergence[1:] -= flux_s[:-1]
        divergence[:, 1:] -= flux_e[:, :-1]
        work[valid] += gamma * divergence[valid]
    work[invalid] = np.nan
    return work
