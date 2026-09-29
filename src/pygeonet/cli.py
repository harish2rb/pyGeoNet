from __future__ import annotations

import argparse
import hashlib
import json
import logging
import platform
from importlib.metadata import version
from pathlib import Path
from typing import Any

import numpy as np

from . import __version__
from .benchmark import run_benchmarks
from .config import PipelineConfig
from .network import network_to_geojson
from .pipeline import process_dem
from .raster import Raster, read_raster, write_raster
from .terrain import slope


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load(path: Path) -> Raster:
    if path.suffix.lower() == ".npy":
        return Raster(np.load(path, allow_pickle=False), (0, 1, 0, 0, 0, -1), None, None)
    return read_raster(path)


def _process(args: argparse.Namespace) -> None:
    input_path = Path(args.dem)
    source = _load(input_path)
    config = PipelineConfig.from_json(args.config) if args.config else PipelineConfig()
    result = process_dem(
        source.data,
        config=config,
        dx=source.dx,
        dy=source.dy,
        transform=source.transform,
        crs=source.crs,
    )
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    arrays = {
        "filtered_dem": result.filtered_dem,
        "conditioned_dem": result.conditioned_dem,
        "slope": result.slope,
        "curvature": result.curvature,
        "receivers": result.receivers,
        "accumulation": result.accumulation,
        "basins": result.basins,
        "channel_mask": result.channel_mask,
    }
    for name, array in arrays.items():
        np.save(output / f"{name}.npy", array)
        if args.geotiff and Path(args.dem).suffix.lower() != ".npy":
            dtype = "int32" if name in {"receivers", "basins"} else "float32"
            driver = "COG" if args.raster_format == "cog" else "GTiff"
            write_raster(output / f"{name}.tif", array, source, dtype=dtype, driver=driver)
    _write_json(output / "network.geojson", network_to_geojson(result.network))
    _write_json(
        output / "run.json",
        {
            "version": __version__,
            "input": str(input_path.resolve()),
            "input_sha256": _sha256(input_path),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pysheds": version("pysheds"),
            "numba": version("numba"),
            "config": config.to_dict(),
            "shape": list(source.data.shape),
            "dx": source.dx,
            "dy": source.dy,
            "crs_present": source.crs is not None,
            "raster_format": args.raster_format if args.geotiff else None,
            "flow_backend": result.flow_backend,
            "channel_heads": result.channel_heads,
            "network_nodes": result.network.number_of_nodes(),
            "network_edges": result.network.number_of_edges(),
        },
    )
    print(
        json.dumps(
            {"output": str(output), "network_edges": result.network.number_of_edges()}, indent=2
        )
    )


def _synthetic(args: argparse.Namespace) -> None:
    size = args.size
    y, x = np.mgrid[-1 : 1 : complex(size), -1 : 1 : complex(size)]
    # Southward regional slope plus a concave valley centered on x=0.
    dem = 1000.0 - 80.0 * y + 120.0 * x * x
    np.save(args.output, dem)
    print(args.output)


def _validate() -> None:
    y, x = np.mgrid[:25, :31]
    plane = 7.0 + 2.0 * x + 3.0 * y
    measured = slope(plane)
    error = float(np.max(np.abs(measured - np.hypot(2, 3))))
    report = {"version": __version__, "plane_slope_max_error": error, "ok": error < 1e-12}
    print(json.dumps(report, indent=2))
    if not report["ok"]:
        raise SystemExit(1)


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="pygeonet", description="Modern pyGeoNet terrain-network tools"
    )
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("--verbose", action="store_true")
    commands = parser.add_subparsers(dest="command", required=True)

    process = commands.add_parser("process", help="run a DEM-to-network pipeline")
    process.add_argument("dem", help="input .npy or GeoTIFF")
    process.add_argument("--config", help="JSON PipelineConfig overrides")
    process.add_argument("--output", default="pygeonet-output")
    process.add_argument(
        "--geotiff", action="store_true", help="also write GeoTIFF arrays with Rasterio"
    )
    process.add_argument(
        "--raster-format",
        choices=("gtiff", "cog"),
        default="gtiff",
        help="GeoTIFF layout: tiled GTiff or cloud-optimized COG (default: gtiff)",
    )

    synthetic = commands.add_parser(
        "synthetic", help="create a legally redistributable analytical valley DEM"
    )
    synthetic.add_argument("--size", type=int, default=129)
    synthetic.add_argument("--output", default="synthetic-valley.npy")

    benchmark = commands.add_parser(
        "benchmark", help="measure slope runtime and peak Python memory"
    )
    benchmark.add_argument("--sizes", nargs="+", type=int, default=[128, 256, 512])
    benchmark.add_argument("--output", default="benchmark.json")
    commands.add_parser("validate-installation", help="run a dependency-free analytical smoke test")

    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s %(name)s: %(message)s",
    )
    if args.command == "process":
        _process(args)
    elif args.command == "synthetic":
        _synthetic(args)
    elif args.command == "benchmark":
        records = run_benchmarks(tuple(args.sizes))
        _write_json(Path(args.output), records)
        print(json.dumps(records, indent=2))
    else:
        _validate()


if __name__ == "__main__":
    main()
