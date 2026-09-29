"""Smoke-test the isolated GRASS 8.5 legacy hydrology reference environment."""

from __future__ import annotations

import json

import grass.script as gs


def raster_stats(name: str) -> dict[str, float]:
    values = gs.parse_command("r.univar", flags="g", map=name)
    return {key: float(value) for key, value in values.items()}


def main() -> None:
    gs.run_command("g.region", n=9, s=0, e=9, w=0, rows=9, cols=9)
    gs.mapcalc("dem = 1000 - row() * 10 + abs(col() - 5) * 2", overwrite=True)
    records: dict[str, dict[str, float]] = {}
    for method, flags in (("legacy-mfd", "a"), ("comparison-d8", "sa")):
        accumulation = f"acc_{method.replace('-', '_')}"
        drainage = f"dir_{method.replace('-', '_')}"
        gs.run_command(
            "r.watershed",
            flags=flags,
            elevation="dem",
            accumulation=accumulation,
            drainage=drainage,
            threshold=3,
            convergence=5,
            overwrite=True,
        )
        accumulation_stats = raster_stats(accumulation)
        direction_stats = raster_stats(drainage)
        if accumulation_stats["max"] < 9 or direction_stats["n"] != 81:
            raise RuntimeError(f"unexpected GRASS output for {method}")
        records[method] = {
            "accumulation_max": accumulation_stats["max"],
            "accumulation_sum": accumulation_stats["sum"],
            "direction_cells": direction_stats["n"],
        }
    version = gs.parse_command("g.version", flags="g")
    if not version["version"].startswith("8.5"):
        raise RuntimeError(f"expected GRASS 8.5, found {version['version']}")
    print(json.dumps({"grass_version": version["version"], "results": records}, sort_keys=True))


if __name__ == "__main__":
    main()
