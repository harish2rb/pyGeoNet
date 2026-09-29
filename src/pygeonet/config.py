from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from numbers import Integral, Real
from pathlib import Path
from typing import Any, Literal


@dataclass(frozen=True, slots=True)
class PipelineConfig:
    """Explicit configuration for the modern pipeline."""

    filter_iterations: int = 10
    filter_quantile: float = 0.90
    diffusion_gamma: float = 0.15
    diffusion_method: Literal["perona-malik-1", "perona-malik-2", "none"] = "perona-malik-2"
    curvature_method: Literal["geometric", "laplacian"] = "geometric"
    flow_threshold_cells: float = 50.0
    curvature_sigma: float = 1.0
    require_positive_curvature: bool = True
    condition_depressions: bool = True
    connectivity: Literal[4, 8] = 8
    flow_backend: Literal["pysheds", "native"] = "pysheds"

    def __post_init__(self) -> None:
        if (
            not isinstance(self.filter_iterations, Integral)
            or isinstance(self.filter_iterations, bool)
            or self.filter_iterations < 0
        ):
            raise ValueError("filter_iterations must be non-negative")
        if not _finite_real(self.filter_quantile) or not 0 < self.filter_quantile < 1:
            raise ValueError("filter_quantile must lie strictly between 0 and 1")
        if not _finite_real(self.diffusion_gamma) or not 0 < self.diffusion_gamma <= 0.25:
            raise ValueError("diffusion_gamma must be in (0, 0.25] for explicit stability")
        if not _finite_real(self.flow_threshold_cells) or self.flow_threshold_cells < 1:
            raise ValueError("flow_threshold_cells must be at least one cell")
        if not _finite_real(self.curvature_sigma) or self.curvature_sigma < 0:
            raise ValueError("curvature_sigma must be non-negative")
        if (
            not isinstance(self.connectivity, Integral)
            or isinstance(self.connectivity, bool)
            or self.connectivity not in (4, 8)
        ):
            raise ValueError("connectivity must be 4 or 8")
        if self.diffusion_method not in {"perona-malik-1", "perona-malik-2", "none"}:
            raise ValueError(f"unsupported diffusion method: {self.diffusion_method}")
        if self.curvature_method not in {"geometric", "laplacian"}:
            raise ValueError(f"unsupported curvature method: {self.curvature_method}")
        if self.flow_backend not in {"pysheds", "native"}:
            raise ValueError(f"unsupported flow backend: {self.flow_backend}")
        if not isinstance(self.require_positive_curvature, bool):
            raise ValueError("require_positive_curvature must be boolean")
        if not isinstance(self.condition_depressions, bool):
            raise ValueError("condition_depressions must be boolean")
        if self.flow_backend == "pysheds" and self.connectivity != 8:
            raise ValueError("pysheds D8 routing requires connectivity=8")

    @classmethod
    def from_json(cls, path: str | Path) -> PipelineConfig:
        loaded: Any = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(loaded, dict):
            raise ValueError("pipeline configuration must be a JSON object")
        data: dict[str, Any] = loaded
        return cls(**data)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _finite_real(value: object) -> bool:
    return isinstance(value, Real) and not isinstance(value, bool) and math.isfinite(float(value))
