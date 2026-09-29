# Performance report

Measured 2026-09-28 on the local x86_64 macOS 13.7 host (Darwin 22.6.0) with Python 3.11.15
and NumPy 2.2.6. The raw,
machine-readable records are in `benchmark-results.json`. Timing uses `perf_counter` and peak
memory uses `tracemalloc`, which measures Python-managed allocations and is not total RSS.

| DEM | Loop reference | Vector cold | Vector steady | Cold speedup | Steady speedup | Vector peak |
|---:|---:|---:|---:|---:|---:|---:|
| 128² | 0.2896 s | 0.00120 s | 0.00058 s | 242x | 498x | 0.63 MiB |
| 256² | 1.2356 s | 0.00196 s | 0.00104 s | 629x | 1190x | 2.13 MiB |
| 512² | 7.9391 s | 0.00747 s | 0.00509 s | 1063x | 1559x | 8.25 MiB |

Steady-state values are medians of five runs; cold and loop values are single runs. Interior
numerical maximum absolute error is zero. Whole-array checksums differ because the
transparent loop reference intentionally leaves its one-cell boundary undefined, whereas
`numpy.gradient` uses one-sided boundary differences. This benchmark demonstrates one scoped
optimization only; it is not an upstream end-to-end speedup claim.

The selectable native D8 receiver implementation is also vectorized one neighbor direction at a time to bound
temporary memory instead of allocating an eight-layer cube. Accumulation uses an O(n)
topological queue. Priority-Flood remains `O(n log n)` and Python heap overhead is a likely
large-raster bottleneck. Full pipeline stage/RSS benchmarks, tiled halo semantics, and out-of-
core routing remain release-roadmap items. Pysheds is now the default, but no speed comparison
is claimed: its eager Numba compilation adds process-start latency that must be reported
separately from steady-state routing before a meaningful benchmark is published.
