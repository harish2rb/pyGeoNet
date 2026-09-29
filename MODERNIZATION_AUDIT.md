# Modernization audit

Audit basis: upstream commit `d29805f`, inspected on 2026-09-27. Every tracked Python file
was inventoried; the V2.1/V3 pairs were hash-compared and the implementation-bearing V3
modules, monoliths, entry points, tests, configuration, I/O, and license were inspected.
Identical or near-identical version copies are grouped below rather than repeating findings.

## Critical findings

| Priority | Affected files | Finding | Implemented response | Residual risk / success criterion |
|---|---|---|---|---|
| P0 security | V2/V3 `pygeonet_fast_marching.py` | `eval(defaults.reciprocalLocalCostFn)` executes arbitrary configured Python. | Modern code exposes a typed `local_speed` function; no expression evaluation. | Any future expression language must use a fixed parser/AST allowlist. `rg 'eval\(' src` must stay empty. |
| P0 data loss | V2/V3 `pygeonet_flow_accumulation.py` | Recursively deletes a fixed `~/grassdata/geonet` location before GRASS setup. | Legacy code is non-importable archival material; native backend uses no external workspace. | A GRASS adapter must use a caller-selected temporary workspace and never recursively delete an unresolved path. |
| P0 reproducibility | all legacy entry points/config files | Import-time globals, current-working-directory discovery, hard-coded Windows/macOS paths, and environment mutation control results. | Frozen `PipelineConfig`, explicit arguments, structured run manifest, no import-time I/O. | Configuration migration cannot promise old-result parity without saved original configuration and environment. |
| P0 scientific | V3 flow/FMM/network modules | Published-style workflow is split across GRASS MFD routing, `r.stream.basins`, scikit-fmm, and custom greedy path descent. None can run in the available modern environment as-is. | Original is retained; component formulas are characterized. Pysheds D8 is the explicit non-equivalent default; GRASS 8.5 is isolated to validation. | Establish reference outputs in the containerized legacy-validation environment before making any parity claim. |
| P1 correctness | V3 `pygeonet_nonlinear_filter.py` | `geonet_diffusion` uses east coefficient with north gradient (`Ce * In`), has an invalid `originalDemArray[:0]` slice, deletes variables conditionally, and Python-3 float slice widths. | New vectorized, stable, NoData-aware four-neighbor Perona–Malik implementation; separate faithful `anisodiff` characterization. | Do not claim equality to the defective `geonet_diffusion`; compare against the actually called V3 `anisodiff`. |
| P1 correctness | raster/vector I/O | Assumes square north-up pixels, drops rotation terms, uses pixel corners inconsistently, mislabels projected x/y as latitude/longitude, and does not consistently write NoData. | Rasterio adapter retains full affine/CRS/masks/NoData/tags; pixel centers apply all six affine coefficients. Atomic tiled GeoTIFF and COG output are supported. | Real-terrain reprojection and large-raster behavior remain to be validated. |
| P1 compatibility | monoliths, tests, V2 modules | Python 2 `print`, `xrange`, integer-division assumptions, `time.clock`, `np.asscalar`, `np.Inf`, `ndimage.measurements`, and old pip internals. | Modern package targets Python 3.11 syntax/APIs; legacy excluded from lint/build/import. | CI must actually pass before claiming 3.12/3.13 support. |
| P1 topology | network delineation/vector output | Greedy raster walks, broad exception handling, implicit transposes, and line export do not expose a formal graph model. | Directed NetworkX output with embedded geometry; only 1-in/1-out chains are simplified. | Confluence geometry and GRASS geodesic equivalence need real-data validation. |

## File-by-file inventory

- `pygeonet/pygeonet_processing*.py` and `pygeonet_v1.0/`: Python 2 monoliths mixing
  GDAL, GRASS, algorithms, plots, mutable globals, and output. They duplicate later modules,
  contain large procedural functions, and cannot be safely imported on Python 3.
- `pygeonet_defaults.py` / `prepare_pygeonet_defaults.py`: module-level switches and a
  string cost function; thresholds have changed between versions without a migration record.
- `pygeonet_prepare.py` / `prepare_pygeonet_inputs.py`: absolute local paths, assumed file
  names, duplicated assignments, and output paths created from CWD.
- `pygeonet_nonlinear_filter.py`: Gaussian boundary construction assumes width five;
  legacy anisotropic diffusion is vectorized but float32 and contains unexplained edge rules.
- `pygeonet_slope_curvature.py`: assumes square pixels; converts geometric singularities to
  zero; imports plotting/statsmodels in the numerical module; QQ threshold automation is TODO.
- `pygeonet_flow_accumulation.py`: GRASS executable paths are fixed to obsolete releases;
  `shell=True`, `sys.exit`, global environment, external add-on dependency, destructive
  cleanup, redundant GeoTIFF serialization, and output rereads prevent library use.
- `pygeonet_fast_marching.py`: unsafe evaluation, writes from computational functions,
  uses removed NumPy names, catches incomplete exception sets, and conflates speed/cost terms.
- `pygeonet_channel_head_definition.py`: nested full-raster Python loops, data-derived
  histogram indexing that can fail for few components, plotting and writes in the algorithm.
- `pygeonet_network_delineation.py`: deeply stateful greedy walk, repeated allocations,
  fragile boundary masks, bare `except`, GRASS side route, and pandas used as an intermediate.
- `pygeonet_xsbank_extraction.py` / `compute_local_extremas.py`: Python-2 division used as
  indices/ranges, hidden parameters assigned globally, empty-peak cases, and pixel-only
  geometry. This capability is preserved only in legacy pending validation.
- `pygeonet_rasterio.py` / `pygeonet_vectorio.py`: direct old GDAL/OGR API, no context
  managers, square/north-up transform assumptions, silent Float32 conversion, no robust
  creation checks, and ambiguous coordinate field naming.
- `pygeonet_plot.py`: computation imports plotting; global figure numbers and interactive
  calls make headless and parallel use fragile.
- `run_pygeonet_processing*.py`: `time.clock`, warnings globally suppressed, all stages
  hard-wired, output after every stage, and no error boundary or machine-readable manifest.
- old tests: scripts rather than isolated assertions, absolute paths, large bundled DEMs,
  print-based inspection, obsolete APIs, and no deterministic expected outputs.
- `checkmod.py`: imports a removed private pip API. Package metadata now declares dependencies.
- tracked `.pyc`, `.DS_Store`, `Thumbs.db`, duplicate TIFFs: supply-chain noise and unclear
  provenance. Retained in `legacy/` or upstream data for historical fidelity, excluded from
  distributions by the `src` build and from modern validation.

## Dependency and compliance findings

The repository supplies a GPL v3 license text but few file headers and no package metadata.
No license was invented or changed. The original DEM files have no per-file provenance or
license record, so modern tests generate analytical arrays and do not redistribute those DEMs
through wheels. The default runtime now includes pysheds 0.5 and its numerical/geospatial
dependencies; Intel macOS markers constrain NumPy, Numba, LLVMlite, and Rasterio to available
wheels. Scikit-fmm and GRASS are not application dependencies. GRASS remains validation-only
and FMM remains documented legacy-only.

## Measurable release criteria

1. Build wheel/sdist with a PEP 517 frontend and install the wheel into a clean environment.
2. Ruff, strict mypy, unit/scientific tests, and dependency audit pass in CI.
3. Plane slope error < 1e-12; interior paraboloid Laplacian equals 4 within floating tolerance.
4. NoData remains masked; non-square pixels and full affine pixel-center coordinates pass tests.
5. D8 accumulation conserves source-cell area and detects receiver cycles.
6. Graph simplification preserves confluences and never invents an edge across disconnected data.
7. Benchmarks record real runtime, Python-level peak memory, numerical error, hardware, and
   cold/steady phase. No unmeasured speed claim enters release notes.
8. Legacy parity remains explicitly “not established” until reference artifacts exist.
