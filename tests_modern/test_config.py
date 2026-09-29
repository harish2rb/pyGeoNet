import json
from pathlib import Path

import pytest

from pygeonet.config import PipelineConfig


def test_default_backend_is_pysheds() -> None:
    assert PipelineConfig().flow_backend == "pysheds"


@pytest.mark.parametrize(
    ("field", "value"),
    (("curvature_sigma", -1), ("connectivity", 6)),
)
def test_config_rejects_invalid_routing_parameters(field: str, value: float) -> None:
    with pytest.raises(ValueError):
        PipelineConfig(**{field: value})  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("field", "value"),
    (
        ("flow_threshold_cells", float("nan")),
        ("curvature_sigma", float("inf")),
        ("filter_iterations", 1.5),
        ("diffusion_method", "unknown"),
        ("curvature_method", "unknown"),
        ("flow_backend", "unknown"),
        ("condition_depressions", "yes"),
    ),
)
def test_config_rejects_nonfinite_and_wrong_runtime_types(field: str, value: object) -> None:
    with pytest.raises(ValueError):
        PipelineConfig(**{field: value})  # type: ignore[arg-type]


def test_config_json_must_be_an_object(tmp_path: Path) -> None:
    path = tmp_path / "config.json"
    path.write_text(json.dumps(["not", "an", "object"]), encoding="utf-8")
    with pytest.raises(ValueError, match="JSON object"):
        PipelineConfig.from_json(path)
