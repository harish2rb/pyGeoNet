from __future__ import annotations

import platform
import tracemalloc
from dataclasses import asdict, dataclass
from statistics import median
from time import perf_counter

import numpy as np
from numpy.typing import NDArray

from .terrain import slope


@dataclass(frozen=True, slots=True)
class BenchmarkRecord:
    size: int
    implementation: str
    seconds: float
    peak_megabytes: float
    checksum: float
    python: str
    numpy: str
    platform: str
    machine: str
    max_abs_error_interior: float


def _loop_slope(dem: NDArray[np.float64]) -> NDArray[np.float64]:
    result = np.full_like(dem, np.nan, dtype=float)
    for row in range(1, dem.shape[0] - 1):
        for col in range(1, dem.shape[1] - 1):
            dz_dy = (dem[row + 1, col] - dem[row - 1, col]) / 2
            dz_dx = (dem[row, col + 1] - dem[row, col - 1]) / 2
            result[row, col] = np.hypot(dz_dx, dz_dy)
    return result


def run_benchmarks(sizes: tuple[int, ...] = (128, 256, 512)) -> list[dict[str, object]]:
    """Measure vectorized and transparent loop-reference slope implementations."""
    records = []
    for size in sizes:
        y, x = np.mgrid[:size, :size]
        dem = 0.001 * x * x + 0.002 * y * y + np.sin(x / 13)
        reference: NDArray[np.float64] | None = None
        for name, function in (
            ("loop-reference", _loop_slope),
            ("vectorized-cold", slope),
            ("vectorized-steady", slope),
        ):
            durations: list[float] = []
            peaks: list[int] = []
            repeats = 5 if name == "vectorized-steady" else 1
            result = np.empty_like(dem)
            for _ in range(repeats):
                tracemalloc.start()
                started = perf_counter()
                result = function(dem)
                durations.append(perf_counter() - started)
                _, peak = tracemalloc.get_traced_memory()
                peaks.append(peak)
                tracemalloc.stop()
            seconds = median(durations)
            peak = max(peaks)
            if reference is None:
                reference = result
            difference = float(np.nanmax(np.abs(result[1:-1, 1:-1] - reference[1:-1, 1:-1])))
            records.append(
                asdict(
                    BenchmarkRecord(
                        size,
                        name,
                        seconds,
                        peak / 1024**2,
                        float(np.nansum(result)),
                        platform.python_version(),
                        np.__version__,
                        platform.platform(),
                        platform.machine(),
                        difference,
                    )
                )
            )
    return records
