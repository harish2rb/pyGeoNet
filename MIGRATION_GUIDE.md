# Migration guide

| Legacy concept | Modern equivalent |
|---|---|
| edit `prepare_pygeonet_inputs.py` | pass CLI paths and `PipelineConfig` JSON |
| mutate `prepare_pygeonet_defaults` | construct frozen `PipelineConfig` |
| `compute_dem_slope(a, pixelScale)` | `slope(a, dx=..., dy=...)` |
| `compute_dem_curvature` | `curvature(..., method="geometric"|"laplacian")` |
| `anisodiff` | `perona_malik`; exact formula characterization remains in `legacy.anisodiff` |
| GRASS flow accumulation | default pysheds D8 adapter; native D8 remains selectable |
| shapefile side effects | returned NetworkX graph + explicit GeoJSON export |
| implicit transposes | always `(row, column)` arrays and `(x, y)` geometry |

The pysheds D8 routing substitution changes science. Users needing published-workflow comparability
must retain their legacy environment and compare reference arrays; this release cannot migrate
those outputs automatically. Thresholds expressed in cells remain resolution-sensitive.

Legacy sources moved, without content rewrites, beneath `legacy/` so they no longer shadow the
installable package. Imports must change from sibling script imports to public `pygeonet` APIs.
