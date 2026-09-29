from __future__ import annotations

import math
from heapq import heappop, heappush

import numpy as np
from numpy.typing import NDArray

from .terrain import FloatArray, as_dem

IntArray = NDArray[np.int64]

_OFFSETS = np.asarray(
    [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)],
    dtype=np.int8,
)


def priority_flood(dem: NDArray[np.floating], *, connectivity: int = 8) -> FloatArray:
    """Fill depressions to their lowest spill elevation using Priority-Flood.

    NoData cells are barriers. Valid cells adjacent to the array edge or NoData seed the
    queue, which makes masked internal holes open boundaries by explicit design.
    """
    array = as_dem(dem)
    if connectivity not in (4, 8):
        raise ValueError("connectivity must be 4 or 8")
    offsets = _OFFSETS if connectivity == 8 else _OFFSETS[[1, 3, 4, 6]]
    valid = np.isfinite(array)
    visited = ~valid.copy()
    filled = array.copy()
    queue: list[tuple[float, int, int]] = []
    rows, cols = array.shape
    for row in range(rows):
        for col in range(cols):
            if not valid[row, col]:
                continue
            is_boundary = row in (0, rows - 1) or col in (0, cols - 1)
            if not is_boundary:
                is_boundary = any(not valid[row + int(dr), col + int(dc)] for dr, dc in offsets)
            if is_boundary:
                visited[row, col] = True
                heappush(queue, (float(filled[row, col]), row, col))
    while queue:
        elevation, row, col = heappop(queue)
        for dr, dc in offsets:
            nr, nc = row + int(dr), col + int(dc)
            if 0 <= nr < rows and 0 <= nc < cols and not visited[nr, nc]:
                visited[nr, nc] = True
                filled[nr, nc] = max(float(array[nr, nc]), elevation)
                heappush(queue, (float(filled[nr, nc]), nr, nc))
    return filled


def d8_receivers(dem: NDArray[np.floating], *, dx: float = 1.0, dy: float = 1.0) -> IntArray:
    """Return flattened receiver index per cell, or -1 for outlets/NoData.

    Only strictly downslope moves are routed. Flat resolution is intentionally not hidden;
    depression-filled flats remain outlets in this initial native backend.
    """
    if not np.isfinite(dx) or not np.isfinite(dy) or dx <= 0 or dy <= 0:
        raise ValueError("pixel spacings must be positive")
    array = as_dem(dem)
    rows, cols = array.shape
    receivers = np.full(array.size, -1, dtype=np.int64)
    distances = np.asarray([math.hypot(dc * dx, dr * dy) for dr, dc in _OFFSETS])
    receiver_grid = receivers.reshape(array.shape)
    best_slope = np.zeros(array.shape, dtype=np.float64)
    for offset_index, (raw_dr, raw_dc) in enumerate(_OFFSETS):
        dr, dc = int(raw_dr), int(raw_dc)
        row_start, row_stop = max(0, -dr), min(rows, rows - dr)
        col_start, col_stop = max(0, -dc), min(cols, cols - dc)
        source_slice = (slice(row_start, row_stop), slice(col_start, col_stop))
        target_slice = (slice(row_start + dr, row_stop + dr), slice(col_start + dc, col_stop + dc))
        source = array[source_slice]
        target = array[target_slice]
        candidate = (source - target) / distances[offset_index]
        rr, cc = np.mgrid[row_start + dr : row_stop + dr, col_start + dc : col_stop + dc]
        candidate_index = rr * cols + cc
        current_slope = best_slope[source_slice]
        current_receiver = receiver_grid[source_slice]
        valid = np.isfinite(source) & np.isfinite(target) & (candidate > 0)
        improve = valid & (
            (candidate > current_slope)
            | (
                (candidate == current_slope)
                & ((current_receiver < 0) | (candidate_index < current_receiver))
            )
        )
        current_slope[improve] = candidate[improve]
        current_receiver[improve] = candidate_index[improve]
    return receiver_grid


def flow_accumulation(
    receivers: IntArray,
    *,
    valid_mask: NDArray[np.bool_] | None = None,
    cell_area: float = 1.0,
) -> FloatArray:
    """Accumulate D8 contributing area in O(n) time via a topological queue."""
    receiver = np.asarray(receivers, dtype=np.int64)
    if receiver.ndim != 2 or not np.isfinite(cell_area) or cell_area <= 0:
        raise ValueError("receivers must be 2-D and cell_area must be positive")
    valid = (
        np.ones(receiver.shape, dtype=bool) if valid_mask is None else np.asarray(valid_mask, bool)
    )
    if valid.shape != receiver.shape:
        raise ValueError("valid_mask shape must match receivers")
    flat = receiver.ravel()
    valid_flat = valid.ravel()
    if np.any(flat[valid_flat] < -1):
        raise ValueError("valid cells may use only -1 as the outlet receiver sentinel")
    indegree = np.zeros(flat.size, dtype=np.int64)
    for source_index in np.flatnonzero(valid_flat):
        target = flat[source_index]
        if target >= 0:
            if target >= flat.size or not valid_flat[target]:
                raise ValueError("receiver points outside the valid DEM")
            indegree[target] += 1
    queue: list[int] = [int(index) for index in np.flatnonzero(valid_flat & (indegree == 0))]
    accumulation = np.where(valid_flat, cell_area, np.nan).astype(np.float64)
    head = 0
    processed = 0
    while head < len(queue):
        queue_source = queue[head]
        head += 1
        processed += 1
        target = flat[queue_source]
        if target >= 0:
            accumulation[target] += accumulation[queue_source]
            indegree[target] -= 1
            if indegree[target] == 0:
                queue.append(int(target))
    if processed != int(valid_flat.sum()):
        raise ValueError("receiver graph contains a cycle")
    return accumulation.reshape(receiver.shape)


def basin_labels(receivers: IntArray, valid_mask: NDArray[np.bool_]) -> IntArray:
    """Label each valid cell by its terminal outlet using path compression."""
    receiver_grid = np.asarray(receivers, dtype=np.int64)
    valid_grid = np.asarray(valid_mask, dtype=bool)
    if receiver_grid.ndim != 2 or valid_grid.shape != receiver_grid.shape:
        raise ValueError("receivers must be 2-D and valid_mask must have the same shape")
    receiver = receiver_grid.ravel()
    valid = valid_grid.ravel()
    valid_receivers = receiver[valid]
    if np.any(valid_receivers < -1):
        raise ValueError("valid cells may use only -1 as the outlet receiver sentinel")
    downstream = valid_receivers[valid_receivers >= 0]
    if np.any(downstream >= receiver.size) or np.any(~valid[downstream]):
        raise ValueError("receiver points outside the valid DEM")
    labels = np.full(receiver.size, -1, dtype=np.int64)
    next_label = 0
    for start in np.flatnonzero(valid):
        if labels[start] >= 0:
            continue
        path: list[int] = []
        node = int(start)
        seen: set[int] = set()
        while node >= 0 and labels[node] < 0 and node not in seen:
            seen.add(node)
            path.append(node)
            node = int(receiver[node])
        if node in seen:
            raise ValueError("receiver graph contains a cycle")
        if node >= 0:
            label = int(labels[node])
        else:
            label = next_label
            next_label += 1
        for cell in path:
            labels[cell] = label
    return labels.reshape(receiver_grid.shape)
