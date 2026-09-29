# Public API

- `PipelineConfig`: validated immutable pipeline parameters.
- `process_dem`: in-memory orchestration returning `PipelineResult`.
- `slope`, `curvature`: coordinate-aware derivatives.
- `perona_malik`, `diffusion_threshold`: nonlinear filtering.
- `route_dem`: routing boundary; defaults to pysheds D8 and also exposes native D8.
- `priority_flood`, `d8_receivers`, `flow_accumulation`, `basin_labels`: retained native
  routing primitives.
- `channel_skeleton`, `channel_heads`: feature extraction.
- `build_network`, `simplify_network`, `network_to_geojson`: graph construction/export.
- `read_raster`, `write_raster`: optional Rasterio adapter with mask/metadata retention,
  atomic writes, safe typed NoData, tiled compression, and COG output.
- `legacy`: component characterization only, not a compatibility pipeline.

Every numerical function accepts arrays and explicit parameters. File I/O is confined to the
raster and CLI boundaries. Exceptions are raised for invalid dimensions, units, routing
cycles, and missing optional dependencies.
