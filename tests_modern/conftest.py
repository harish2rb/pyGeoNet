import os

# pysheds eagerly compiles its Numba kernels at import. Unit tests exercise the same functions
# without JIT so the matrix stays fast; a separate CI smoke job exercises JIT-enabled execution.
os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
