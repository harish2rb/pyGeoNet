from __future__ import annotations

import logging
from dataclasses import dataclass

import networkx as nx
import numpy as np
from numpy.typing import NDArray

from .config import PipelineConfig
from .features import channel_heads, channel_skeleton
from .filtering import diffusion_threshold, perona_malik
from .network import build_network, simplify_network
from .routing import route_dem
from .terrain import FloatArray, as_dem, curvature, slope

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class PipelineResult:
    filtered_dem: FloatArray
    conditioned_dem: FloatArray
    slope: FloatArray
    curvature: FloatArray
    receivers: NDArray[np.int64]
    accumulation: FloatArray
    basins: NDArray[np.int64]
    channel_mask: NDArray[np.bool_]
    channel_heads: tuple[tuple[int, int], ...]
    network: nx.DiGraph[int]
    flow_backend: str


def process_dem(
    dem: NDArray[np.floating],
    *,
    config: PipelineConfig | None = None,
    dx: float = 1.0,
    dy: float = 1.0,
    transform: tuple[float, float, float, float, float, float] | None = None,
    crs: str | None = None,
) -> PipelineResult:
    """Execute the modern workflow entirely in memory."""
    cfg = config or PipelineConfig()
    source = as_dem(dem)
    if transform is None:
        transform = (0.0, dx, 0.0, 0.0, 0.0, -dy)
    if cfg.diffusion_method == "none" or cfg.filter_iterations == 0:
        filtered = source.copy()
    else:
        kappa = diffusion_threshold(source, cfg.filter_quantile, dx=dx, dy=dy)
        LOGGER.info("Perona-Malik threshold kappa=%g", kappa)
        filtered = perona_malik(
            source,
            iterations=cfg.filter_iterations,
            kappa=kappa,
            gamma=cfg.diffusion_gamma,
            dx=dx,
            dy=dy,
            method=cfg.diffusion_method,
        )
    routing = route_dem(
        filtered,
        backend=cfg.flow_backend,
        condition_depressions=cfg.condition_depressions,
        connectivity=cfg.connectivity,
        dx=dx,
        dy=dy,
    )
    slope_array = slope(filtered, dx=dx, dy=dy)
    curvature_array = curvature(filtered, dx=dx, dy=dy, method=cfg.curvature_method)
    channel = channel_skeleton(
        routing.accumulation / (dx * dy),
        curvature_array,
        flow_threshold=cfg.flow_threshold_cells,
        curvature_sigma=cfg.curvature_sigma,
        require_positive_curvature=cfg.require_positive_curvature,
    )
    heads = tuple(channel_heads(channel, routing.receivers))
    graph = simplify_network(
        build_network(
            channel,
            routing.receivers,
            transform=transform,
            crs=crs,
            algorithm=routing.algorithm,
        )
    )
    LOGGER.info(
        "extracted %d channel cells, %d heads, %d simplified edges",
        int(channel.sum()),
        len(heads),
        graph.number_of_edges(),
    )
    return PipelineResult(
        filtered_dem=filtered,
        conditioned_dem=routing.conditioned_dem,
        slope=slope_array,
        curvature=curvature_array,
        receivers=routing.receivers,
        accumulation=routing.accumulation,
        basins=routing.basins,
        channel_mask=channel,
        channel_heads=heads,
        network=graph,
        flow_backend=routing.algorithm,
    )
