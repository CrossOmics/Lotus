from importlib import import_module
from typing import Any

__all__ = ["core_selection"]


def __getattr__(name: str) -> Any:
    if name == "core_selection":
        module = import_module("lotus.core_analysis._cplearn_core")
        value = getattr(module, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module 'lotus.core_analysis' has no attribute {name!r}")
