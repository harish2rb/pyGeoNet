from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from uuid import uuid4

import numpy as np
from numpy.typing import NDArray

from .terrain import FloatArray

Transform = tuple[float, float, float, float, float, float]


@dataclass(frozen=True, slots=True)
class Raster:
    """A dependency-light, single-band raster and its spatial metadata."""

    data: FloatArray
    transform: Transform
    crs: str | None
    nodata: float | None
    tags: tuple[tuple[str, str], ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if self.data.ndim != 2:
            raise ValueError("raster data must be two-dimensional")
        if len(self.transform) != 6 or not np.all(np.isfinite(self.transform)):
            raise ValueError("transform must contain six finite GDAL affine coefficients")
        if (
            self.dx == 0
            or self.dy == 0
            or self.transform[1] * self.transform[5] == self.transform[2] * self.transform[4]
        ):
            raise ValueError("transform must define two non-degenerate pixel axes")

    @property
    def dx(self) -> float:
        return float(np.hypot(self.transform[1], self.transform[4]))

    @property
    def dy(self) -> float:
        return float(np.hypot(self.transform[2], self.transform[5]))


def read_raster(path: str | Path, *, band: int = 1) -> Raster:
    """Read one GeoTIFF band while retaining CRS, affine, NoData, mask, and tags."""
    try:
        import rasterio
    except ImportError as exc:  # pragma: no cover - package dependency guard
        raise ImportError(
            "GeoTIFF support requires Rasterio; reinstall pygeonet-modernized"
        ) from exc
    if band < 1:
        raise ValueError("band must be a positive one-based index")
    with rasterio.open(path) as dataset:
        if band > dataset.count:
            raise ValueError(f"band {band} does not exist in {dataset.count}-band raster")
        data = dataset.read(band, masked=True).astype(np.float64).filled(np.nan)
        return Raster(
            data=data,
            transform=tuple(dataset.transform.to_gdal()),
            crs=dataset.crs.to_wkt() if dataset.crs else None,
            nodata=dataset.nodata,
            tags=tuple(sorted(dataset.tags(band).items())),
        )


def _output_nodata(array: NDArray[np.generic], dtype: np.dtype[Any], reference: Raster) -> float:
    finite_values = array[np.isfinite(array)]
    if np.issubdtype(dtype, np.floating):
        candidate = reference.nodata
        float_limits = np.finfo(dtype)
        if (
            candidate is not None
            and (np.isnan(candidate) or float_limits.min <= candidate <= float_limits.max)
            and (np.isnan(candidate) or not np.any(finite_values == candidate))
        ):
            return float(candidate)
        return float("nan")
    if not np.issubdtype(dtype, np.integer):
        raise ValueError(f"unsupported raster dtype: {dtype}")
    int_limits = np.iinfo(dtype)
    for candidate in (reference.nodata, int_limits.min, int_limits.max):
        if (
            candidate is not None
            and float(candidate).is_integer()
            and int_limits.min <= candidate <= int_limits.max
            and not np.any(finite_values == candidate)
        ):
            return float(candidate)
    raise ValueError(f"no non-colliding NoData value is available for {dtype}")


def write_raster(
    path: str | Path,
    data: NDArray[np.generic],
    reference: Raster,
    *,
    dtype: str = "float32",
    driver: str = "GTiff",
    compress: str = "deflate",
) -> None:
    """Atomically write a masked, compressed GeoTIFF or Cloud Optimized GeoTIFF.

    The source raster's invalid-data footprint, affine transform, CRS, and band tags are
    retained. ``driver`` must be ``"GTiff"`` or ``"COG"``.
    """
    try:
        import rasterio
        from affine import Affine
    except ImportError as exc:  # pragma: no cover - package dependency guard
        raise ImportError(
            "GeoTIFF support requires Rasterio; reinstall pygeonet-modernized"
        ) from exc
    array = np.asarray(data)
    if array.ndim != 2 or array.shape != reference.data.shape:
        raise ValueError("output data must be two-dimensional and match the reference shape")
    if driver not in {"GTiff", "COG"}:
        raise ValueError("driver must be 'GTiff' or 'COG'")
    target_dtype = np.dtype(dtype)
    valid = np.isfinite(reference.data) & np.isfinite(array)
    nodata = _output_nodata(array, target_dtype, reference)
    if np.issubdtype(target_dtype, np.integer) and np.any(valid):
        limits = np.iinfo(target_dtype)
        valid_values = array[valid]
        if np.any(valid_values < limits.min) or np.any(valid_values > limits.max):
            raise OverflowError(f"values cannot be represented as {target_dtype}")
        if np.any(valid_values != np.floor(valid_values)):
            raise ValueError(f"fractional values cannot be represented as {target_dtype}")
    elif np.issubdtype(target_dtype, np.floating) and np.any(valid):
        float_limits = np.finfo(target_dtype)
        if np.any(np.abs(array[valid]) > float_limits.max):
            raise OverflowError(f"values cannot be represented as {target_dtype}")
    written = np.where(valid, array, nodata).astype(target_dtype)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.name}.{uuid4().hex}.tmp.tif")
    creation_options: dict[str, object] = {
        "driver": driver,
        "height": array.shape[0],
        "width": array.shape[1],
        "count": 1,
        "dtype": target_dtype.name,
        "crs": reference.crs,
        "transform": Affine.from_gdal(*reference.transform),
        "nodata": nodata,
        "compress": compress,
    }
    if driver == "GTiff":
        creation_options.update(tiled=True, blockxsize=256, blockysize=256)
    else:
        creation_options.update(blocksize=512, overview_resampling="nearest")
    try:
        with (
            rasterio.Env(GDAL_TIFF_INTERNAL_MASK=True),
            rasterio.open(temporary, "w", **creation_options) as dataset,
        ):
            dataset.write(written, 1)
            dataset.write_mask(valid.astype(np.uint8) * 255)
            if reference.tags:
                dataset.update_tags(1, **dict(reference.tags))
            dataset.update_tags(software="pygeonet-modernized")
        temporary.replace(target)
    finally:
        temporary.unlink(missing_ok=True)
