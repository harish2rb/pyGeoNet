# Architecture

```text
Raster I/O (optional) -> Raster(data, affine, CRS, NoData)
                              |
                         PipelineConfig
                              |
 DEM -> filtering -> terrain derivatives -> routing backend (pysheds D8 default)
                                                        |-> accumulation/basins
                                                        |-> channel candidates
                                                        |-> embedded DiGraph
                                                               |-> GeoJSON
```

`terrain`, `filtering`, `flow`, and `features` are array-only numerical modules. They do not
read files, plot, mutate global configuration, or write intermediate rasters. `raster` is the
optional Rasterio boundary and owns validation, metadata/mask retention, atomic tiled GeoTIFF,
and COG serialization. `network` owns pixel-to-world conversion and topology-preserving
graph simplification. `pipeline` orchestrates stages and returns all arrays for inspection.
`legacy` captures runnable numerical formulas but is not an end-to-end compatibility backend.

`routing` is the only hydrology backend boundary. It defaults to pysheds D8, converts pysheds
direction codes to the package's explicit flattened-receiver representation, and retains the
native Priority-Flood/D8 implementation as a selectable reference. GRASS is outside this graph:
its pinned 8.5 container is used only by the legacy-validation workflow.

Coordinates use NumPy `(row, column)` inside arrays and `(x, y)` in graph geometry. `dx` is
column spacing and `dy` is row spacing. Pixel geometry refers to centers. A GDAL-order affine
tuple `(x0, a, b, y0, d, e)` is retained without assuming north-up or square pixels.

Future algorithms should implement explicit stage protocols and identify themselves in
metadata. Experimental uncertainty-aware reconstruction belongs behind a non-default interface;
it must never add connections during ordinary graph simplification.
