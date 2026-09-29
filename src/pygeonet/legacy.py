"""Pure numerical characterization of formulas found in upstream V3.

This module is not the default pipeline. It records independently executable pieces of
the legacy science; GRASS MFD/D-infinity and scikit-fmm parity remain unestablished.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from .terrain import FloatArray, as_dem


def curvature_square_pixels(
    dem: NDArray[np.floating], pixel_scale: float, *, geometric: bool = True
) -> FloatArray:
    """Port of V3 ``compute_dem_curvature`` without its NaN-to-zero side effect."""
    array = as_dem(dem)
    grad_row, grad_col = np.gradient(array, pixel_scale)
    if geometric:
        magnitude = np.hypot(grad_row, grad_col)
        grad_row = np.divide(grad_row, magnitude, out=np.zeros_like(grad_row), where=magnitude > 0)
        grad_col = np.divide(grad_col, magnitude, out=np.zeros_like(grad_col), where=magnitude > 0)
    grad_grad_row, _ = np.gradient(grad_row, pixel_scale)
    _, grad_grad_col = np.gradient(grad_col, pixel_scale)
    result = grad_grad_row + grad_grad_col
    result[~np.isfinite(array)] = np.nan
    return np.asarray(result, dtype=np.float64)


def anisodiff(
    dem: NDArray[np.floating],
    *,
    iterations: int,
    kappa: float,
    gamma: float,
    step: tuple[float, float] = (1.0, 1.0),
    option: int = 2,
) -> FloatArray:
    """Faithful vectorized characterization of the V3 ``anisodiff`` function."""
    output = as_dem(dem).copy()
    delta_s = np.zeros_like(output)
    delta_e = np.zeros_like(output)
    ns = np.zeros_like(output)
    ew = np.zeros_like(output)
    for _ in range(iterations):
        delta_s.fill(0)
        delta_e.fill(0)
        delta_s[:-1] = np.diff(output, axis=0)
        delta_e[:, :-1] = np.diff(output, axis=1)
        if option == 1:
            g_s = np.exp(-((delta_s / kappa) ** 2)) / step[0]
            g_e = np.exp(-((delta_e / kappa) ** 2)) / step[1]
        elif option == 2:
            g_s = 1.0 / (1.0 + (delta_s / kappa) ** 2) / step[0]
            g_e = 1.0 / (1.0 + (delta_e / kappa) ** 2) / step[1]
        else:
            raise ValueError("legacy option must be 1 or 2")
        south, east = g_s * delta_s, g_e * delta_e
        ns[:] = south
        ew[:] = east
        ns[1:] -= south[:-1]
        ew[:, 1:] -= east[:, :-1]
        both_nan = np.isnan(ns) & np.isnan(ew)
        ns[np.isnan(ns)] = 0
        ew[np.isnan(ew)] = 0
        update = ns + ew
        update[both_nan] = np.nan
        output += gamma * update
    return output


def local_speed(
    flow: FloatArray, skeleton: NDArray[np.bool_], normalized_curvature: FloatArray
) -> FloatArray:
    """Safe equivalent of the default legacy ``eval`` expression."""
    finite_flow = flow[np.isfinite(flow)]
    flow_mean = float(finite_flow.mean())
    return flow + flow_mean * np.asarray(skeleton) + flow_mean * normalized_curvature
