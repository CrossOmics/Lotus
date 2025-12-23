'''
User APIs
'''
from ._preprocess_core import (
    # Core Pipeline
    run_preprocessing,

    # Quality Control & Filtering
    calculate_qc_metrics,
    filter_cells,
    filter_genes,
    downsample_counts,
    sample,

    # Doublet Detection
    scrublet,
    scrublet_simulate_doublets,

    # Normalization & Transformation
    normalization,
    log1p,
    scale,

    # Feature Selection & Correction
    highly_variable_genes,
    regress_out,
    combat,

    # Dimensionality Reduction & Neighborhood
    pca,
    neighbors,

    # Recipes (Pre-defined workflows)
    recipe_zheng17,
    recipe_weinreb17,
    recipe_seurat,
)

__all__ = [
    # Core Pipeline
    "run_preprocessing",

    # Quality Control & Filtering
    "calculate_qc_metrics",
    "filter_cells",
    "filter_genes",
    "downsample_counts",
    "sample",

    # Doublet Detection
    "scrublet",
    "scrublet_simulate_doublets",

    # Normalization & Transformation
    "normalization",
    "log1p",
    "scale",

    # Feature Selection & Correction
    "highly_variable_genes",
    "regress_out",
    "combat",

    # Dimensionality Reduction & Neighborhood
    "pca",
    "neighbors",

    # Recipes
    "recipe_zheng17",
    "recipe_weinreb17",
    "recipe_seurat",
]
