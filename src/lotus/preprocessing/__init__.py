from importlib import import_module
from typing import Any

__all__ = [
    "calculate_qc_metrics",
    "filter_cells",
    "filter_genes",
    "downsample_counts",
    "sample",
    "scrublet",
    "scrublet_simulate_doublets",
    "normalize_total",
    "log1p",
    "scale",
    "highly_variable_genes",
    "regress_out",
    "combat",
    "pca",
    "neighbors",
    "recipe_zheng17",
    "recipe_weinreb17",
    "recipe_seurat",
]

_CORE_EXPORTS = set(__all__)


def __getattr__(name: str) -> Any:
    if name in _CORE_EXPORTS:
        module = import_module("lotus.preprocessing._preprocess_core")
        value = getattr(module, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module 'lotus.preprocessing' has no attribute {name!r}")
