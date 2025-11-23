from __future__ import annotations

from .clustering import clustering, run_clustering
from .core_analysis import compute_coremap_embedding, core_analysis
from .deg import marker_genes, run_differential_expression
from .preprocess import (
    filtering,
    hvg,
    neighbors,
    normalization,
    pca,
    preprocess,
    qc,
    scaling,
)
from .visualization import coremap, render_visualizations, umap

__all__ = [
    # Preprocess steps
    "qc",
    "filtering",
    "normalization",
    "hvg",
    "scaling",
    "pca",
    "neighbors",
    "preprocess",
    # Visualization
    "umap",
    "coremap",
    "render_visualizations",
    # Clustering
    "clustering",
    "run_clustering",
    # DEG
    "marker_genes",
    "run_differential_expression",
    # CoreAnalysis
    "core_analysis",
    "compute_coremap_embedding",  # Backward compatibility alias
]
