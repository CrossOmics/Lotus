from importlib import import_module
from typing import Any

__all__ = [
    "read",
    "write",
    "read_10x_h5",
    "read_10x_mtx",
    "read_visium",
    "read_h5ad",
    "read_csv",
    "read_excel",
    "read_hdf",
    "read_loom",
    "read_mtx",
    "read_text",
    "read_umi_tools",
]

_CORE_EXPORTS = set(__all__)


def __getattr__(name: str) -> Any:
    if name in _CORE_EXPORTS:
        module = import_module("lotus.io._io_core")
        value = getattr(module, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module 'lotus.io' has no attribute {name!r}")
