# pyGeoNet modernization release candidate

This branch modernizes the original GPL-3.0 pyGeoNet research code into an installable,
tested Python package for terrain derivatives and drainage-network extraction from DEMs.
It preserves the complete upstream source under `legacy/` and provides a separate modern
implementation under `src/pygeonet`.

> **Scientific status:** the default backend uses pysheds 0.5 depression/flat handling and
> D8 routing. It is analytically tested, but it is **not equivalent** to the original GRASS
> MFD + fast-marching/geodesic workflow. Do not interpret a pysheds-D8 result as a
> reproduction of a published GeoNet result. See
> [SCIENTIFIC_VALIDATION.md](SCIENTIFIC_VALIDATION.md).

## What it does

The working modern pipeline:

1. reads a NumPy DEM or GeoTIFF;
2. preserves masks, NoData, CRS, affine transforms, band tags, and separate x/y pixel sizes;
3. optionally applies Perona–Malik nonlinear diffusion;
4. computes slope and geometric or Laplacian curvature;
5. conditions depressions and resolves flats with pysheds;
6. computes pysheds D8 receivers and contributing area, then explicit basin labels;
7. thresholds a candidate channel skeleton;
8. exports a topology-preserving NetworkX graph and GeoJSON.

## Installation

Python 3.11 was clean-install tested locally. CI is configured for 3.11, 3.12, and 3.13;
3.12/3.13 should not be treated as verified until that CI runs successfully. On the local
Intel macOS host, dependency markers select NumPy below 2.4, Numba 0.62.1, LLVMlite 0.45.1,
and Rasterio below 1.4.4 so every compiled dependency resolves from a wheel. Rasterio 1.4.3
uses bundled GDAL 3.9.3; no system GDAL or LLVM installation is required.

```bash
python3.11 -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[dev]'
pygeonet validate-installation
```

Pysheds and Rasterio are installed by default. The empty `geo` extra remains only as an
installation-compatibility alias for prerelease users.

Conda is useful where Rasterio/GDAL wheels are unsuitable:

```bash
conda create -n pygeonet-modern python=3.12 numpy networkx rasterio numba scipy
conda activate pygeonet-modern
python -m pip install pysheds==0.5
python -m pip install -e . --no-deps
```

## Reproducible tutorial

The tutorial creates an analytical DEM rather than downloading data of uncertain provenance.

```bash
pygeonet synthetic --size 129 --output synthetic-valley.npy
pygeonet process synthetic-valley.npy --config configs/pysheds-d8.json --output tutorial-output
```

Outputs are `.npy` stage arrays, `network.geojson`, and `run.json`. Because `.npy` has no
georeferencing, its coordinates are pixel-space centers. For a GeoTIFF input, pass
`--geotiff` to write tiled, Deflate-compressed GeoTIFF stages too, or add
`--raster-format cog` for Cloud Optimized GeoTIFFs. Full instructions are in
[docs/TUTORIAL.md](docs/TUTORIAL.md).

Programmatic use:

```python
import numpy as np
from pygeonet import PipelineConfig, process_dem

y, x = np.mgrid[-1:1:129j, -1:1:129j]
dem = 1000 - 80 * y + 120 * x**2
config = PipelineConfig(flow_threshold_cells=30, require_positive_curvature=False)
result = process_dem(dem, config=config, dx=2.0, dy=3.0)
print(result.channel_heads, result.network.number_of_edges())
```

## Commands

```text
pygeonet process DEM [--config JSON] [--output DIR] [--geotiff] [--raster-format gtiff|cog]
pygeonet synthetic [--size N] [--output FILE.npy]
pygeonet benchmark [--sizes 128 256 512] [--output benchmark.json]
pygeonet validate-installation
```

## GRASS and legacy compatibility

The original code requires old GRASS 7 modules, GDAL bindings, scikit-fmm, mutable module
globals, and platform-specific configuration. It remains in `legacy/` for provenance and
characterization, but is not imported or installed. GRASS 8.5 is restricted to a
scheduled/containerized legacy-validation workflow under `validation/grass`; it is not an
application dependency or selectable production backend. The old source describes its GRASS
result as D-infinity, but the invoked `r.watershed` mode is default MFD. A true historical-parity
claim still requires reviewed reference outputs.

## Routing backends

`PipelineConfig()` defaults to `flow_backend="pysheds"`. This release uses pysheds D8 because
the downstream graph has exactly one receiver per cell. Pysheds MFD and D-infinity require a
weighted multi-receiver result model and are not silently collapsed into a dominant D8 path.
Set `flow_backend="native"` or use `configs/native-d8.json` for the retained native reference
implementation.

Pysheds eagerly compiles Numba kernels on first import, so the first routing process can have
noticeable startup latency. Unit tests disable JIT for speed; a separate CI job runs the compiled
default backend.

## Development and verification

```bash
ruff check .
ruff format --check .
mypy
pytest
python -m build
pygeonet benchmark --output benchmark-results.json
```

Architecture, migration, audit, performance, and release details are in:

- [MODERNIZATION_AUDIT.md](MODERNIZATION_AUDIT.md)
- [ARCHITECTURE.md](ARCHITECTURE.md)
- [SCIENTIFIC_VALIDATION.md](SCIENTIFIC_VALIDATION.md)
- [PERFORMANCE_REPORT.md](PERFORMANCE_REPORT.md)
- [MIGRATION_GUIDE.md](MIGRATION_GUIDE.md)
- [FINAL_REPORT.md](FINAL_REPORT.md)

## Citation, provenance, and license

Please cite Sangireddy et al. (2016), *GeoNet: An open source software for the automatic
and objective extraction of channel heads, channel network, and channel morphology from
high resolution topography data*, DOI `10.1016/j.envsoft.2016.04.026`, and the foundational
method paper listed in [CITATION.cff](CITATION.cff).

The upstream `License.txt` is retained verbatim. This modernization does not relicense the
project. Bundled upstream DEMs have unclear per-file provenance and are excluded from the
modern tests and tutorial.
