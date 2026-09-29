from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import networkx as nx
import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True, slots=True)
class NetworkMetadata:
    crs: str | None
    transform: tuple[float, float, float, float, float, float]
    algorithm: str
    source_shape: tuple[int, int]


def pixel_center(
    row: int, col: int, transform: tuple[float, float, float, float, float, float]
) -> tuple[float, float]:
    """Apply GDAL-style affine coefficients to the center of a pixel."""
    x0, a, b, y0, d, e = transform
    x = x0 + a * (col + 0.5) + b * (row + 0.5)
    y = y0 + d * (col + 0.5) + e * (row + 0.5)
    return float(x), float(y)


def build_network(
    channel_mask: NDArray[np.bool_],
    receivers: NDArray[np.integer],
    *,
    transform: tuple[float, float, float, float, float, float] = (0, 1, 0, 0, 0, -1),
    crs: str | None = None,
    algorithm: str = "native-d8",
) -> nx.DiGraph[int]:
    """Create a directed embedded channel graph without adding inferred connections."""
    channel = np.asarray(channel_mask, dtype=bool)
    receiver = np.asarray(receivers, dtype=np.int64)
    if channel.ndim != 2 or receiver.shape != channel.shape:
        raise ValueError("channel mask must be 2-D and match the receiver shape")
    if len(transform) != 6 or not np.all(np.isfinite(transform)):
        raise ValueError("transform must contain six finite affine coefficients")
    if transform[1] * transform[5] == transform[2] * transform[4]:
        raise ValueError("transform must define two non-degenerate pixel axes")
    channel_receivers = receiver.ravel()[channel.ravel()]
    if np.any(channel_receivers < -1) or np.any(channel_receivers >= channel.size):
        raise ValueError("channel receiver index is outside the raster")
    graph: nx.DiGraph[int] = nx.DiGraph()
    source_shape = (int(channel.shape[0]), int(channel.shape[1]))
    graph.graph.update(metadata=NetworkMetadata(crs, transform, algorithm, source_shape))
    _, cols = channel.shape
    for flat_index in np.flatnonzero(channel.ravel()):
        row, col = divmod(int(flat_index), cols)
        graph.add_node(int(flat_index), row=row, col=col, pos=pixel_center(row, col, transform))
    for source in list(graph.nodes):
        target = int(receiver.ravel()[source])
        if target in graph:
            a, b = graph.nodes[source]["pos"], graph.nodes[target]["pos"]
            graph.add_edge(
                source, target, geometry=[a, b], length=float(np.hypot(b[0] - a[0], b[1] - a[1]))
            )
    return graph


def simplify_network(graph: nx.DiGraph[int]) -> nx.DiGraph[int]:
    """Collapse only unambiguous in-degree=out-degree=1 chains, preserving geometry."""
    result = graph.copy()
    changed = True
    while changed:
        changed = False
        for node in list(result.nodes):
            if result.in_degree(node) == 1 and result.out_degree(node) == 1:
                predecessor = next(result.predecessors(node))
                successor = next(result.successors(node))
                if predecessor == successor or result.has_edge(predecessor, successor):
                    continue
                first = result.edges[predecessor, node]
                second = result.edges[node, successor]
                geometry = [*first["geometry"], *second["geometry"][1:]]
                result.add_edge(
                    predecessor,
                    successor,
                    geometry=geometry,
                    length=first["length"] + second["length"],
                )
                result.remove_node(node)
                changed = True
                break
    return result


def network_to_geojson(graph: nx.DiGraph[int]) -> dict[str, Any]:
    metadata: NetworkMetadata = graph.graph["metadata"]
    features = []
    for u, v, attrs in graph.edges(data=True):
        features.append(
            {
                "type": "Feature",
                "properties": {
                    "source": u,
                    "target": v,
                    "length": attrs["length"],
                    "algorithm": metadata.algorithm,
                },
                "geometry": {"type": "LineString", "coordinates": attrs["geometry"]},
            }
        )
    return {
        "type": "FeatureCollection",
        "name": "pyGeoNet channel network",
        "crs_wkt": metadata.crs,
        "features": features,
    }
