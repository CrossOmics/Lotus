from importlib import import_module
from typing import Any

__all__ = [
    "annotation_analysis",
    "clustering",
    "core_analysis",
    "io",
    "lineagetracker",
    "preprocessing",
    "trajectory",
    "visualization",
]


def __getattr__(name: str) -> Any:
    if name in __all__:
        module = import_module(f"lotus.{name}")
        globals()[name] = module
        return module
    raise AttributeError(f"module 'lotus' has no attribute {name!r}")
