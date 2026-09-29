from pathlib import Path

import numpy as np
import pytest

from pygeonet.raster import Raster, read_raster, write_raster


def test_geotiff_roundtrip_preserves_spatial_metadata(tmp_path: Path) -> None:
    rasterio = pytest.importorskip("rasterio")
    data = np.arange(20, dtype=float).reshape(4, 5)
    data[1, 2] = np.nan
    reference = Raster(
        data=data,
        transform=(500_000.0, 2.0, 0.25, 4_500_000.0, 0.1, -3.0),
        crs="EPSG:32610",
        nodata=-9999.0,
        tags=(("source", "analytical-test"),),
    )
    path = tmp_path / "roundtrip.tif"
    write_raster(path, data, reference)
    restored = read_raster(path)
    assert np.allclose(restored.data, data, equal_nan=True)
    assert np.allclose(restored.transform, reference.transform)
    assert restored.crs is not None and "UTM zone 10N" in restored.crs
    assert ("source", "analytical-test") in restored.tags
    with rasterio.open(path) as dataset:
        assert dataset.profile["tiled"] is True
        assert dataset.compression.name.lower() == "deflate"


def test_integer_output_uses_noncolliding_nodata_and_source_mask(tmp_path: Path) -> None:
    rasterio = pytest.importorskip("rasterio")
    source = np.ones((3, 4), dtype=float)
    source[0, 0] = np.nan
    reference = Raster(source, (100, 1, 0, 200, 0, -1), None, None)
    receivers = np.array([[-1, -1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11]], dtype=np.int64)
    path = tmp_path / "receivers.tif"
    write_raster(path, receivers, reference, dtype="int32")
    with rasterio.open(path) as dataset:
        restored = dataset.read(1, masked=True)
        assert dataset.nodata == np.iinfo(np.int32).min
        assert bool(restored.mask[0, 0])
        assert restored[0, 1] == -1


def test_cog_output_is_cloud_optimized(tmp_path: Path) -> None:
    rasterio = pytest.importorskip("rasterio")
    data = np.arange(4096, dtype=float).reshape(64, 64)
    reference = Raster(data, (10, 2, 0, 20, 0, -2), "EPSG:3857", None)
    path = tmp_path / "terrain-cog.tif"
    write_raster(path, data, reference, driver="COG")
    with rasterio.open(path) as dataset:
        assert dataset.driver == "GTiff"
        assert dataset.tags(ns="IMAGE_STRUCTURE").get("LAYOUT") == "COG"
        assert np.array_equal(dataset.read(1), data.astype(np.float32))


def test_integer_output_rejects_fractional_values(tmp_path: Path) -> None:
    pytest.importorskip("rasterio")
    data = np.ones((3, 3), dtype=float)
    data[1, 1] = 1.5
    reference = Raster(data, (100, 1, 0, 200, 0, -1), None, None)
    with pytest.raises(ValueError, match="fractional"):
        write_raster(tmp_path / "invalid.tif", data, reference, dtype="int32")
