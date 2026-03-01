from importlib import import_module
from typing import Any

__all__ = [
    "diffmap",
    "draw_graph",
    "dpt",
    "paga",
    "embedding_density",
    "sim",
]

_CORE_EXPORTS = set(__all__)


def __getattr__(name: str) -> Any:
    if name in _CORE_EXPORTS:
        module = import_module("lotus.trajectory._trajectory_core")
        value = getattr(module, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module 'lotus.trajectory' has no attribute {name!r}")
