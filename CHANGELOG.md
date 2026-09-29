# Changelog

## 0.2.0rc2 - 2026-09-28

- Make pysheds 0.5 D8 the default conditioned-routing backend while retaining native D8 as an
  explicit reference backend.
- Add a pinned GRASS 8.5 container workflow for legacy MFD/D8 validation only.
- Promote validated Rasterio support into the default dependency set.

## 0.2.0rc1 - 2026-09-27

- Preserve all upstream implementations under `legacy/`.
- Add PEP 517 packaging and a typed `src/pygeonet` library.
- Add array-only terrain, filtering, conditioning, native D8, feature, and graph stages.
- Add optional Rasterio I/O with atomic tiled GeoTIFF/COG writes, retained spatial metadata
  and masks, and non-colliding typed NoData values.
- Add GeoJSON export and four CLI commands.
- Add analytical scientific tests, benchmark records, lint/type/build configuration, and CI.
- Remove unsafe expression evaluation and implicit workspace deletion from the modern path.
- Explicitly mark original GRASS/FMM parity and real-landscape validation as incomplete.
