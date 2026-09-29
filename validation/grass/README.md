# GRASS 8.5 legacy validation

GRASS is an isolated validation oracle, not an application dependency or selectable production
backend. The workflow runs the legacy-equivalent default MFD mode (`r.watershed -a`) and an
explicit D8 comparison (`r.watershed -sa`) in an automatically deleted temporary project.

The workflow checks out source on the GitHub-hosted runner, then mounts it read-only into the
official GRASS 8.5.0 Ubuntu image pinned by immutable registry digest. Container networking is
disabled because validation is fully local; this also avoids depending on the image's CA bundle.
The smoke test verifies that environment and both routing modes. Scientific
reference-array comparison remains gated on a legally documented DEM and reviewed expected
outputs. The legacy `r.stream.basins` add-on is intentionally excluded until its source revision
is pinned in a derived image.
