# Release process

1. Confirm CI has passed on every declared Python version and supported platform.
2. Run scientific tests and benchmarks; update reports only from saved machine-readable data.
3. Review dependency audit, license metadata, wheel contents, and `twine check dist/*`.
4. Install the built wheel into a fresh environment and run `pygeonet validate-installation`.
5. Update semantic version and changelog, then tag a reviewed commit.
6. Publishing to PyPI or creating a public GitHub release requires explicit owner authorization.
