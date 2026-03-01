from importlib import import_module
from typing import Any

__all__ = [
    "leiden",
    "louvain",
    "dendrogram",
    "umap",
    "tsne",
]

_CORE_EXPORTS = set(__all__)


def __getattr__(name: str) -> Any:
    if name in _CORE_EXPORTS:
        module = import_module("lotus.clustering.clustering_core")
        value = getattr(module, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module 'lotus.clustering' has no attribute {name!r}")
