from __future__ import annotations

import importlib
import sys
from types import ModuleType
from typing import Any

from ...methods.deg import find_markers as _scRNA_seq_find_markers
from ...methods.deg.find_markers import DEAnalyzer, DEOptions, de_from_model, de_from_adata

# Import wrapped functions from external module
from .external.cplearn import corespect, coremap_embedding

# Create a namespace object for backward compatibility
class _CplearnNamespace:
    """Namespace object for cplearn functions and classes."""
    def __init__(self, module_getattr):
        self._module_getattr = module_getattr
        self.corespect = corespect
        self.coremap_embedding = coremap_embedding
        # Add DEG analysis classes and functions
        self.DEOptions = DEOptions
        self.DEAnalyzer = DEAnalyzer
        self.de_from_model = de_from_model
        self.de_from_adata = de_from_adata
        self.find_markers = find_markers
    
    def __getattr__(self, name: str) -> Any:
        # Delegate to module-level __getattr__ for dynamic imports
        # This allows access to CorespectModel, Coremap, etc.
        return self._module_getattr(name)

__all__ = [
    "cplearn",
    "coremap",
    "corespect",
    "coremap_embedding",
    "utils",
    "anchored_map",
    "Coremap",
    "CorespectModel",
    "FlowRankConfig",
    "StableCoreConfig",
    "FineGrainedConfig",
    "ClusterConfig",
    "PropagationConfig",
    "find_anchors",
    "cluster_core",
    "stable_core",
    "fine_grained_core",
    "propagate_from_core",
    "FlowRank",
    "find_markers",
    "DEAnalyzer",
    "DEOptions",
    "de_from_model",
    "de_from_adata",
]

_IMPORT_ERROR: ModuleNotFoundError | None = None
_MODULE_CACHE: dict[str, ModuleType] = {}

_ATTR_SOURCES = {
    "anchored_map": "cplearn.coremap.coremap",
    "Coremap": "cplearn.coremap.coremap",
    "CorespectModel": "cplearn.corespect.corespect",
    "FlowRankConfig": "cplearn.corespect.config",
    "StableCoreConfig": "cplearn.corespect.config",
    "FineGrainedConfig": "cplearn.corespect.config",
    "ClusterConfig": "cplearn.corespect.config",
    "PropagationConfig": "cplearn.corespect.config",
    "find_anchors": "cplearn.coremap.find_anchors",
    "cluster_core": "cplearn.corespect.corespect",
    "stable_core": "cplearn.corespect.corespect",
    "fine_grained_core": "cplearn.corespect.corespect",
    "propagate_from_core": "cplearn.corespect.corespect",
    "FlowRank": "cplearn.corespect.corespect",
}

_SUBPACKAGES = ("coremap", "corespect", "utils")

find_markers = _scRNA_seq_find_markers


def _import_module(name: str) -> ModuleType:
    global _IMPORT_ERROR
    if name in _MODULE_CACHE:
        return _MODULE_CACHE[name]
    try:
        module = importlib.import_module(name)
    except ModuleNotFoundError as exc:  # pragma: no cover - optional dependency
        _IMPORT_ERROR = exc
        raise ModuleNotFoundError(
            "cplearn is not available. Install it in the current environment to "
            "use lotus.cplearn."
        ) from exc
    _MODULE_CACHE[name] = module
    return module


def _ensure_subpackage(name: str) -> ModuleType:
    full_name = f"{__name__}.{name}"
    if full_name in sys.modules:
        return sys.modules[full_name]
    module = _import_module(f"cplearn.{name}")
    sys.modules[full_name] = module
    return module


def _module_getattr(name: str) -> Any:
    """Internal function for module-level attribute access."""
    if name in _SUBPACKAGES:
        return _ensure_subpackage(name)

    # Direct imports (already available)
    if name == "corespect":
        return corespect
    if name == "coremap_embedding":
        return coremap_embedding
    
    # DEG analysis classes and functions
    if name == "DEOptions":
        return DEOptions
    if name == "DEAnalyzer":
        return DEAnalyzer
    if name == "de_from_model":
        return de_from_model
    if name == "de_from_adata":
        return de_from_adata
    if name == "find_markers":
        return find_markers

    module_name = _ATTR_SOURCES.get(name)
    if module_name is None:
        raise AttributeError(f"cplearn has no attribute {name!r}") from None
    module = _import_module(module_name)
    return getattr(module, name)


def __getattr__(name: str) -> Any:
    if name == "cplearn":
        return _cplearn_namespace
    
    return _module_getattr(name)


# Create the namespace object
_cplearn_namespace = _CplearnNamespace(_module_getattr)


def __dir__() -> list[str]:
    return sorted(set(__all__) | set(_SUBPACKAGES))
