import json

import numpy as np
import pytest

from pygeonet.network import build_network, network_to_geojson, pixel_center, simplify_network


def test_pixel_center_uses_affine_without_axis_swap() -> None:
    transform = (100, 2, 0, 200, 0, -3)
    assert pixel_center(1, 2, transform) == (105.0, 195.5)


def test_simplification_preserves_branch_and_geometry() -> None:
    mask = np.zeros((5, 5), dtype=bool)
    mask[0:4, 2] = True
    mask[2, 1:4] = True
    receiver = np.full((5, 5), -1, dtype=int)
    for source, target in [
        ((0, 2), (1, 2)),
        ((1, 2), (2, 2)),
        ((2, 1), (2, 2)),
        ((2, 3), (2, 2)),
        ((2, 2), (3, 2)),
    ]:
        receiver[source] = np.ravel_multi_index(target, mask.shape)
    graph = simplify_network(build_network(mask, receiver))
    junction = np.ravel_multi_index((2, 2), mask.shape)
    assert junction in graph
    assert graph.in_degree(junction) == 3
    payload = network_to_geojson(graph)
    json.dumps(payload)
    assert len(payload["features"]) == 4


def test_network_rejects_out_of_bounds_receiver() -> None:
    mask = np.ones((3, 3), dtype=bool)
    receiver = np.full((3, 3), -1, dtype=int)
    receiver[0, 0] = 99
    with pytest.raises(ValueError, match="outside the raster"):
        build_network(mask, receiver)
