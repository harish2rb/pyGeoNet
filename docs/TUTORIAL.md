# Analytical valley tutorial

Create and process a 129×129 concave valley:

```bash
pygeonet synthetic --size 129 --output synthetic-valley.npy
pygeonet process synthetic-valley.npy --config configs/pysheds-d8.json --output tutorial-output
```

The surface is `z = 1000 - 80y + 120x²`: it has a regional downstream gradient and a
known valley axis at `x=0`. It is generated from an equation and carries no external dataset
license. Inspect `run.json`, load stages with `numpy.load`, and read `network.geojson` with a
GIS or JSON tool. NumPy input uses unit pixel coordinates; it has no CRS.

For a real GeoTIFF, install the `geo` extra and run:

```bash
pygeonet process dem.tif --config configs/pysheds-d8.json --output result --geotiff
```

The default outputs are tiled, Deflate-compressed GeoTIFFs. For HTTP range-friendly Cloud
Optimized GeoTIFF output, add `--raster-format cog`. Both layouts preserve the source mask,
affine transform, CRS, and band tags. Integer stages use a NoData value that does not collide
with the valid `-1` outlet value.

Confirm elevation and horizontal units are compatible before interpreting slope, ensure the
CRS is projected for metric analysis, inspect NoData boundaries, and perform a threshold
sensitivity analysis. Never compare pysheds D8 output directly with legacy GRASS MFD output as
if the routing models were equivalent. Use `configs/native-d8.json` for an explicit comparison
with the retained native implementation.
