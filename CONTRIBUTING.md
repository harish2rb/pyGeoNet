# Contributing

Create focused branches and keep numerical changes separate from formatting. Every algorithm
change needs an analytical or reference-output test, a statement of coordinate/NoData behavior,
and benchmark evidence when performance is claimed. Do not update legacy source content except
for an explicitly reviewed archival correction.

Run `ruff check .`, `ruff format --check .`, `mypy`, `pytest`, and `python -m build`. New routing
or network methods must remain opt-in until analytical validation and an explicit maintainer
decision promote them, must preserve topology metadata, and must describe theoretical differences
from existing methods. Ordinary tests disable Numba JIT; changes to the pysheds adapter must also
pass the JIT-enabled CI smoke job.
