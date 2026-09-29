from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def channel_skeleton(
    accumulation: NDArray[np.floating],
    curvature: NDArray[np.floating],
    *,
    flow_threshold: float,
    curvature_sigma: float = 1.0,
    require_positive_curvature: bool = True,
) -> NDArray[np.bool_]:
    """Legacy-inspired dual threshold (flow AND positive curvature anomaly)."""
    flow = np.asarray(accumulation, dtype=np.float64)
    curve = np.asarray(curvature, dtype=np.float64)
    if flow.shape != curve.shape:
        raise ValueError("accumulation and curvature shapes must match")
    if (
        not np.isfinite(flow_threshold)
        or not np.isfinite(curvature_sigma)
        or flow_threshold < 0
        or curvature_sigma < 0
    ):
        raise ValueError("flow threshold and curvature sigma must be non-negative")
    mask = np.isfinite(flow) & (flow >= flow_threshold)
    if require_positive_curvature:
        finite = curve[np.isfinite(curve)]
        if not finite.size:
            return np.zeros(flow.shape, dtype=bool)
        threshold = float(np.mean(finite) + curvature_sigma * np.std(finite))
        mask &= np.isfinite(curve) & (curve > threshold)
    return mask


def channel_heads(
    channel_mask: NDArray[np.bool_], receivers: NDArray[np.integer]
) -> list[tuple[int, int]]:
    """Return channel cells with no immediate upstream channel donor."""
    channel = np.asarray(channel_mask, dtype=bool)
    receiver = np.asarray(receivers, dtype=np.int64)
    if channel.shape != receiver.shape:
        raise ValueError("channel mask and receiver shapes must match")
    flat_channel, flat_receiver = channel.ravel(), receiver.ravel()
    channel_receivers = flat_receiver[flat_channel]
    if np.any(channel_receivers < -1) or np.any(channel_receivers >= channel.size):
        raise ValueError("channel receiver index is outside the raster")
    has_channel_donor = np.zeros(channel.size, dtype=bool)
    for source in np.flatnonzero(flat_channel):
        target = flat_receiver[source]
        if target >= 0 and flat_channel[target]:
            has_channel_donor[target] = True
    rows, cols = np.unravel_index(np.flatnonzero(flat_channel & ~has_channel_donor), channel.shape)
    return list(zip(rows.tolist(), cols.tolist(), strict=True))
