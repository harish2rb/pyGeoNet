import numpy as np
import pytest

from pygeonet.features import channel_heads, channel_skeleton


def test_positive_curvature_requirement_handles_all_nodata() -> None:
    accumulation = np.ones((3, 3), dtype=float)
    curvature = np.full((3, 3), np.nan)
    assert not channel_skeleton(
        accumulation, curvature, flow_threshold=1, require_positive_curvature=True
    ).any()


def test_channel_heads_rejects_invalid_receiver() -> None:
    channel = np.ones((3, 3), dtype=bool)
    receivers = np.full((3, 3), -1, dtype=int)
    receivers[0, 0] = 100
    with pytest.raises(ValueError, match="outside the raster"):
        channel_heads(channel, receivers)
