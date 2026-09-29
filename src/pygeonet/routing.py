from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from .flow import basin_labels, d8_receivers, flow_accumulation, priority_flood
from .terrain import FloatArray, as_dem

IntArray = NDArray[np.int64]
FlowBackend = Literal["pysheds", "native"]

_PYSHEDS_D8 = {
    64: (-1, 0),
    128: (-1, 1),
    1: (0, 1),
    2: (1, 1),
    4: (1, 0),
    8: (1, -1),
    16: (0, -1),
    32: (-1, -1),
}


@dataclass(frozen=True, slots=True)
class RoutingResult:
    conditioned_dem: FloatArray
    receivers: IntArray
    accumulation: FloatArray
    basins: IntArray
    algorithm: str


def _directions_to_receivers(
    directions: NDArray[np.integer], valid_mask: NDArray[np.bool_]
) -> IntArray:
    """Convert pysheds D8 codes into flattened downstream cell indices."""
    fdir = np.asarray(directions)
    valid = np.asarray(valid_mask, dtype=bool)
    rows, cols = fdir.shape
    receivers = np.full(fdir.shape, -1, dtype=np.int64)
    known = np.isin(fdir, (*_PYSHEDS_D8, -2, -1, 0))
    if np.any(valid & ~known):
        unknown = np.unique(fdir[valid & ~known]).tolist()
        raise ValueError(f"unsupported pysheds D8 direction code(s): {unknown}")
    for code, (dr, dc) in _PYSHEDS_D8.items():
        sources = np.argwhere((fdir == code) & valid)
        if not sources.size:
            continue
        targets = sources + np.asarray((dr, dc))
        inside = (
            (targets[:, 0] >= 0)
            & (targets[:, 0] < rows)
            & (targets[:, 1] >= 0)
            & (targets[:, 1] < cols)
        )
        if not np.all(inside):
            raise ValueError("pysheds direction points outside the DEM")
        target_valid = valid[targets[:, 0], targets[:, 1]]
        if not np.all(target_valid):
            raise ValueError("pysheds direction points into a NoData cell")
        receivers[sources[:, 0], sources[:, 1]] = targets[:, 0] * cols + targets[:, 1]
    return receivers


def _route_pysheds(
    dem: FloatArray, *, condition_depressions: bool, dx: float, dy: float
) -> RoutingResult:
    try:
        from affine import Affine
        from pyproj import CRS
        from pysheds.grid import Grid
        from pysheds.sview import Raster, ViewFinder
    except ImportError as exc:  # pragma: no cover - package dependency guard
        raise ImportError("the default routing backend requires pysheds>=0.5") from exc

    valid = np.isfinite(dem)
    # Routing needs cell-axis lengths, not world orientation. The original affine remains on
    # the pipeline result and network; this normalized internal affine also handles rotated DEMs.
    view = ViewFinder(
        affine=Affine(dx, 0.0, 0.0, 0.0, -dy, 0.0),
        shape=dem.shape,
        nodata=np.nan,
        mask=valid,
        crs=CRS.from_epsg(3857),
    )
    grid = Grid(viewfinder=view)
    prepared = Raster(dem.copy(), viewfinder=view)
    if condition_depressions:
        prepared = grid.fill_pits(prepared, nodata_out=np.nan)
        prepared = grid.fill_depressions(prepared, nodata_out=np.nan)
        prepared = grid.resolve_flats(prepared, nodata_out=np.nan)
    directions = grid.flowdir(prepared, routing="d8", nodata_out=0)
    receivers = _directions_to_receivers(np.asarray(directions), valid)
    accumulation = np.asarray(
        grid.accumulation(directions, routing="d8", nodata_out=np.nan), dtype=np.float64
    )
    accumulation[~valid] = np.nan
    accumulation *= dx * dy
    conditioned = np.asarray(prepared, dtype=np.float64).copy()
    conditioned[~valid] = np.nan
    return RoutingResult(
        conditioned_dem=conditioned,
        receivers=receivers,
        accumulation=accumulation,
        basins=basin_labels(receivers, valid),
        algorithm="pysheds-d8",
    )


def _route_native(
    dem: FloatArray,
    *,
    condition_depressions: bool,
    connectivity: Literal[4, 8],
    dx: float,
    dy: float,
) -> RoutingResult:
    conditioned = (
        priority_flood(dem, connectivity=connectivity) if condition_depressions else dem.copy()
    )
    valid = np.isfinite(conditioned)
    receivers = d8_receivers(conditioned, dx=dx, dy=dy)
    return RoutingResult(
        conditioned_dem=conditioned,
        receivers=receivers,
        accumulation=flow_accumulation(receivers, valid_mask=valid, cell_area=dx * dy),
        basins=basin_labels(receivers, valid),
        algorithm="native-d8",
    )


def route_dem(
    dem: NDArray[np.floating],
    *,
    backend: FlowBackend = "pysheds",
    condition_depressions: bool = True,
    connectivity: Literal[4, 8] = 8,
    dx: float = 1.0,
    dy: float = 1.0,
) -> RoutingResult:
    """Condition and route a DEM with the selected explicit D8 backend."""
    if not np.isfinite(dx) or not np.isfinite(dy) or dx <= 0 or dy <= 0:
        raise ValueError("pixel spacings must be positive")
    source = as_dem(dem)
    if backend == "pysheds":
        if connectivity != 8:
            raise ValueError("pysheds D8 routing requires connectivity=8")
        return _route_pysheds(source, condition_depressions=condition_depressions, dx=dx, dy=dy)
    if backend == "native":
        return _route_native(
            source,
            condition_depressions=condition_depressions,
            connectivity=connectivity,
            dx=dx,
            dy=dy,
        )
    raise ValueError(f"unsupported flow backend: {backend}")
