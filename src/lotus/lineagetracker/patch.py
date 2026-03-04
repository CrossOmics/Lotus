"""Monkey-patches for AnnData to transparently track derivation operations.

Patches are installed once at ``import lotus`` time (via ``__init__.py``).
Each patch wraps the original method, calls it normally, then registers
the resulting AnnData as a child node in the lineage DAG.
"""

import functools

import anndata
from anndata import AnnData

from .tracker import LineageTracker


def _parent_name(tracker: LineageTracker, lid: str) -> str:
    """Resolve a parent lid to its display_name, falling back to the lid prefix."""
    node = tracker._nodes.get(lid)
    if node and node.display_name:
        return node.display_name
    return lid[:8]


def install_patches():
    """Apply all AnnData monkey-patches.  Safe to call multiple times."""
    _patch_getitem()
    _patch_copy()
    _patch_concat()


def _patch_getitem():
    """Intercept ``adata[index]`` (slicing) to create a child node."""
    original = AnnData.__getitem__

    @functools.wraps(original)
    def patched(self, index):
        result = original(self, index)

        tracker = LineageTracker.instance()
        if tracker.is_nested_operation:
            return result
        parent_lid = tracker.get_lid(self)
        if parent_lid:
            tracker.register(
                result,
                parents=[parent_lid],
                creation_op="slice",
                description=f"slice from {_parent_name(tracker, parent_lid)}",
            )
        return result

    AnnData.__getitem__ = patched


def _patch_copy():
    """Intercept ``adata.copy()`` to create a child node."""
    original = AnnData.copy

    @functools.wraps(original)
    def patched(self):
        result = original(self)

        tracker = LineageTracker.instance()
        if tracker.is_nested_operation:
            return result
        parent_lid = tracker.get_lid(self)
        if parent_lid:
            tracker.register(
                result,
                parents=[parent_lid],
                creation_op="copy",
                description=f"copy of {_parent_name(tracker, parent_lid)}",
            )
        return result

    AnnData.copy = patched


def _patch_concat():
    """Intercept ``anndata.concat(...)`` to create a multi-parent node."""
    original_concat = anndata.concat

    @functools.wraps(original_concat)
    def patched(*args, **kwargs):
        result = original_concat(*args, **kwargs)

        tracker = LineageTracker.instance()
        if tracker.is_nested_operation:
            return result

        # The first positional arg (or kwarg "adatas") is the list of inputs
        adatas = args[0] if args else kwargs.get("adatas", [])
        if isinstance(adatas, dict):
            adatas = list(adatas.values())

        # Collect parent lids from all tracked inputs
        parent_lids = []
        for a in adatas:
            lid = tracker.get_lid(a)
            if lid:
                parent_lids.append(lid)

        if parent_lids:
            names = [_parent_name(tracker, l) for l in parent_lids]
            tracker.register(
                result,
                parents=parent_lids,
                creation_op="concat",
                description=f"concat of [{', '.join(names)}]",
            )
        return result

    anndata.concat = patched
