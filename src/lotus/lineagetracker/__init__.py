"""Lineage tracking for AnnData objects.

Importing this package automatically:

1. Monkey-patches ``AnnData.__getitem__``, ``AnnData.copy``, and
   ``anndata.concat`` so that slice / copy / concat operations are
   tracked transparently.
2. Registers an ``atexit`` hook that persists the final DAG to JSON
   and renders it as a PNG in ``<cwd>/lineage-tracer_data/``.
"""

import atexit

from .metadata import LoggingMeta, logged, track_module, tracked
from .models import LineageNode, OperationRecord
from .patch import install_patches
from .tracker import LineageTracker
from .visualizer import visualize

# Apply AnnData monkey-patches on import
install_patches()


def _on_exit():
    """Persist the DAG and render the graph when the process exits."""
    if LineageTracker._instance is None:
        return  # Tracker was never used — nothing to save
    tracker = LineageTracker._instance
    tracker.save()
    try:
        visualize(tracker.lineage_file, tracker.graph_file)
    except Exception:
        pass  # Never crash the exit sequence


atexit.register(_on_exit)

__all__ = [
    "LineageNode",
    "OperationRecord",
    "LineageTracker",
    "LoggingMeta",
    "logged",
    "tracked",
    "track_module",
    "visualize",
]
