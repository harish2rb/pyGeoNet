# Final modernization report

## Implemented

Upstream commit `d29805f` was cloned and work was isolated on
`modernization/py312-core`. Legacy source was preserved under `legacy/`. A new PEP 517 package
provides explicit configuration, structured logging, safe numerical APIs, NoData-aware terrain
derivatives/filtering, default pysheds D8 routing, retained native Priority-Flood/D8 routing,
basin labels,
candidate channels, channel heads, topology-preserving NetworkX construction, full-affine
coordinates, Rasterio I/O, GeoJSON, a CLI, analytical tutorial, tests, benchmarks,
CI, dependency auditing, and release documentation.

The modern path has no `eval`, shell subprocess, global GIS environment, or recursive workspace
deletion. The upstream GPL v3 text and authorship are preserved; no public release was created.

## Verification performed locally

- Python 3.11.15: all 43 tests passed, including pysheds and Rasterio integration tests.
- The pysheds 0.5 default passed in non-JIT unit mode and in a separate locally JIT-compiled
  process; CI has a dedicated JIT-enabled smoke job. Intel macOS uses the wheel-supported
  Numba 0.62.1/LLVMlite 0.45.1 pair, including across the declared Python 3.11-3.13 range.
- Analytical plane, paraboloid, rotation, units, NoData, drainage conservation, topology, and
  deterministic end-to-end checks passed.
- Measured vectorized slope benchmark completed at 128², 256², and 512² with zero interior
  difference from the loop reference; details are in `PERFORMANCE_REPORT.md`.
- Python 3.12/3.13 and non-macOS checks are configured in CI but were not available locally.
- Rasterio 1.4.3 installed from its Intel macOS wheel with bundled GDAL 3.9.3. Round-trip
  metadata/mask retention, non-colliding integer NoData, tiled compression, and COG layout
  were verified locally; the Linux CI job independently exercises geospatial I/O.
- The final wheel installed and passed both `validate-installation` and a COG round-trip in an
  isolated temporary environment. Ruff, formatting, strict mypy, and package build checks passed. `pip-audit`
  found no known vulnerabilities after the environment's pip/setuptools were updated; the
  unpublished local package itself was necessarily skipped by the registry-backed audit.
- Release artifacts: wheel SHA-256
  `4aca94af092dfa093aa6fec03ae96cb59f7e51f888dae2756c258df003e81294` and source archive
  SHA-256 `0182aa2cd2a0a52d815f8500d43c9858582d4e11f180ad458e2ac4ce05c3f512`.

## Known limitations

Legacy GRASS MFD, scikit-fmm, channel-head, cross-section, and bank parity is not established.
Pysheds D8 is not a drop-in scientific replacement for GRASS MFD. Real DEM
accuracy, resolution convergence, reprojection sensitivity, full pipeline RSS profiling, and
out-of-core processing remain incomplete. Bundled upstream DEMs
were not used because their per-file provenance is not documented.

## Prioritized scientific roadmap

1. Pin the configured GRASS 8.5 image by digest and publish hashed MFD/D8 reference arrays for a
   legally documented small DEM.
2. Pin the `r.stream.basins` add-on revision in the validation image and compare basins, FMM,
   and channel paths without exposing GRASS as a production backend.
3. Generalize results to weighted multi-receiver flow before enabling pysheds MFD/D-infinity;
   then compare D8/MFD/D-infinity sensitivity across rotations and resolutions.
4. Add geodesic cost and joint channel-head extraction as experimental strategies with explicit
   uncertainty and no default promotion until field/reference validation.
5. Validate cross sections and bank detection independently; do not port Python-2 indices
   mechanically.
6. Add tiled derivatives/conditioning with documented halos, external-memory routing, and full
   stage RSS/temporary-disk benchmarks on representative large DEMs.
7. Establish curated real-terrain benchmarks with survey/reference networks, documented licenses,
   accuracy metrics, and parameter sensitivity rather than visual-only assessment.
