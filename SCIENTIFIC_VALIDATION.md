# Scientific validation

## What upstream actually implements

The latest upstream V3 entry point calls: a 90th-percentile slope threshold; Perona–Malik
option-2 diffusion; slope magnitude; either divergence-of-unit-gradient (“geometric”) or
Laplacian curvature; GRASS watershed accumulation/directions and stream basins; a dual
flow/curvature threshold; reciprocal local speed

`A + mean(A) * S + mean(A) * K_normalized`;

scikit-fmm travel time from basin outlets; channel heads as local travel-time maxima on
connected skeleton components; and greedy eight-neighbor descent to outlets. Cross sections
and banks are subsequently estimated from centerline-normal profiles. Several alternatives
appear in comments or dead functions; they are not treated as implemented defaults.

## Modern equations

For elevation `z`, pixel spacings `dx, dy`, slope is
`sqrt((dz/dx)^2 + (dz/dy)^2)`. Laplacian curvature is
`d²z/dx² + d²z/dy²`. Legacy geometric curvature is
`div(grad(z) / |grad(z)|)`, with zero contribution where gradient magnitude is zero.

Perona–Malik option 2 uses conductance `g(s)=1/(1+(s/kappa)^2)` and an explicit four-neighbor
flux update. `kappa` defaults to the configured slope quantile. NoData flux is zero and the
original NoData mask is restored.

The default hydrology backend uses pysheds 0.5 to fill pits/depressions, resolve flats, assign
D8 directions, and accumulate contributing area. Direction codes are converted to explicit
receiver indices before basin and network construction. The retained native backend uses
Priority-Flood and steepest-downslope D8; unlike pysheds, its filled flats remain outlets.

The legacy code calls its GRASS result D-infinity, but it invokes `r.watershed` without the
single-flow flag and therefore uses GRASS MFD accumulation. Neither modern D8 backend is
equivalent to that MFD behavior on divergent hillslopes, flats, or cell-scale orientations.
Pysheds MFD/D-infinity are not enabled because the current result and graph model represents
exactly one receiver per cell; silently selecting a dominant branch would change the method.

## Automated analytical evidence

Tests cover:

- planar slope on non-square pixels;
- the known Laplacian (`4`) of `x²+y²`;
- slope rotation behavior;
- linear response to elevation-unit scaling;
- NoData preservation and boundaries;
- parity between square-pixel modern and characterized V3 curvature formulas;
- isolated depression filling;
- D8 area conservation, basin labeling, cycle rejection, and NoData;
- pysheds default selection, pit/flat handling, determinism, NoData, and area conservation;
- full-affine pixel centers and confluence-preserving simplification;
- deterministic end-to-end extraction on a synthetic concave valley.
- GeoTIFF mask, affine, CRS, NoData, tag, compression, integer-sentinel, and COG behavior.

All 43 tests pass on Python 3.11.15. Rasterio 1.4.3 with bundled GDAL 3.9.3 was verified
locally on Intel macOS; geospatial I/O is also isolated in a dedicated Linux CI job.

## Unestablished claims

Original-output parity was not established: the legacy GRASS 7 + add-on + Python 2 environment
is unavailable, upstream tests have no numerical assertions, and bundled DEM provenance is
unclear. Fast marching, GRASS MFD/D∞, legacy channel-head histograms, cross sections, banks,
resolution convergence, reprojection sensitivity, and real-landscape accuracy remain outside
the validated modern core. The GRASS 8.5 container currently verifies executable MFD/D8 modes,
not numerical parity. Pysheds D8 is therefore the supported default, not a claim of scientific
superiority or legacy equivalence.
