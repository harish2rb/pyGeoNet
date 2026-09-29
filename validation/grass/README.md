# GRASS 8.5 legacy validation

GRASS is an isolated validation oracle, not an application dependency or selectable production
backend. The workflow runs the legacy-equivalent default MFD mode (`r.watershed -a`) and an
explicit D8 comparison (`r.watershed -sa`) in an automatically deleted temporary project.

The workflow pins the official GRASS 8.5.0 Ubuntu image by immutable registry digest. The smoke
test verifies that environment and both routing modes. Scientific
reference-array comparison remains gated on a legally documented DEM and reviewed expected
outputs. The legacy `r.stream.basins` add-on is intentionally excluded until its source revision
is pinned in a derived image.
