"""Decorators and metaclass that intercept function / method calls to
record lineage operations and snapshot ``adata.X`` into layers."""

import functools
import inspect

import anndata
from anndata import AnnData
from loguru import logger

from .tracker import LineageTracker


def _extract_adata(*args, **kwargs) -> AnnData | None:
    """Return the first AnnData found in positional or keyword arguments."""
    for arg in args:
        if isinstance(arg, anndata.AnnData):
            return arg
    for value in kwargs.values():
        if isinstance(value, anndata.AnnData):
            return value
    return None


def _build_args_dict(func, args, kwargs) -> dict:
    """Bind positional + keyword arguments to their parameter names.

    Falls back to just ``kwargs`` when introspection fails (e.g. C
    extensions, or functions with unusual signatures).
    """
    try:
        sig = inspect.signature(func)
        bound = sig.bind(*args, **kwargs)
        bound.apply_defaults()
        result = dict(bound.arguments)
        result.pop("self", None)
        return result
    except (ValueError, TypeError):
        return dict(kwargs)


class LoggingMeta(type):
    """Metaclass that wraps every public method with:

    1. ``adata.X`` snapshotting into ``adata.layers``
    2. ``tracker.record_op()`` for lineage tracking
    3. ``loguru`` entry / exit logging
    4. Auto-registration of any new AnnData returned by the method
    """

    def __new__(mcs, name, bases, namespace):
        for attr, value in namespace.items():
            if callable(value) and not attr.startswith("_"):
                namespace[attr] = mcs._wrap(attr, value)
        return super().__new__(mcs, name, bases, namespace)

    @staticmethod
    def _wrap(method_name, method):
        @functools.wraps(method)
        def wrapper(self, *args, **kwargs):
            cls_name = type(self).__name__
            tracker = LineageTracker.instance()

            adata = _extract_adata(*args, **kwargs)
            if adata is not None:
                # Snapshot current X so the user can compare before / after
                adata.layers[f"before_{cls_name}.{method_name}_copy"] = (
                    adata.X.copy()
                )
                # Record the operation with serialised arguments
                all_args = _build_args_dict(method, (self, *args), kwargs)
                tracker.record_op(
                    adata, f"{cls_name}.{method_name}", all_args
                )

            logger.info(f"[BEFORE] {cls_name}.{method_name}")
            result = method(self, *args, **kwargs)
            logger.info(f"[AFTER]  {cls_name}.{method_name}")

            # If the method returned a *new* AnnData (e.g. copy=True),
            # register it as a child of the input adata.
            if isinstance(result, AnnData) and not tracker.get_lid(result):
                parent_lid = tracker.get_lid(adata) if adata else None
                tracker.register(
                    result,
                    parents=[parent_lid] if parent_lid else [],
                    creation_op=f"{cls_name}.{method_name}",
                    description=f"{cls_name}.{method_name}",
                )

            return result

        return wrapper


def logged(func):
    """Decorator for standalone functions — same behaviour as
    :class:`LoggingMeta` but for functions that live outside a class.

    When applied, every call will:

    1. Snapshot ``adata.X`` into ``adata.layers['before_<name>_copy']``
    2. Call ``tracker.record_op()`` with the function name and arguments
    3. Log entry / exit via loguru
    4. Auto-register any new AnnData returned by the function
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        tracker = LineageTracker.instance()

        adata = _extract_adata(*args, **kwargs)
        if adata is not None:
            # Snapshot current X before the function mutates it
            adata.layers[f"before_{func.__name__}_copy"] = adata.X.copy()
            # Record this operation on the node that owns *adata*
            all_args = _build_args_dict(func, args, kwargs)
            tracker.record_op(adata, func.__name__, all_args)

        logger.info(f"[BEFORE] {func.__name__}")
        result = func(*args, **kwargs)
        logger.info(f"[AFTER]  {func.__name__}")

        # If the function returned a brand-new AnnData (e.g. an IO reader
        # or a function called with copy=True), register it in the DAG.
        if isinstance(result, AnnData) and not tracker.get_lid(result):
            parent_lid = tracker.get_lid(adata) if adata else None
            parents = [parent_lid] if parent_lid else []
            # For root nodes (no parent), include the first arg in the
            # description — it is usually the file path for IO functions.
            desc = func.__name__
            if not parents and args:
                desc = f"{func.__name__}('{args[0]}')"
            tracker.register(
                result,
                parents=parents,
                creation_op=func.__name__,
                description=desc,
            )

        return result

    return wrapper


# ── lightweight tracking (no X snapshot, no loguru) ─────────────────

def tracked(func):
    """Minimal decorator that only records lineage — no X snapshot,
    no loguru logging.  Used by :func:`track_module` to auto-wrap
    every public function in a Lotus submodule."""

    @functools.wraps(func)
    def _tracked_wrapper(*args, **kwargs):
        tracker = LineageTracker.instance()

        # Record the operation if an AnnData is among the arguments
        adata = _extract_adata(*args, **kwargs)
        if adata is not None:
            all_args = _build_args_dict(func, args, kwargs)
            tracker.record_op(adata, func.__name__, all_args)

        result = func(*args, **kwargs)

        # Auto-register any brand-new AnnData returned by the function
        if isinstance(result, AnnData) and not tracker.get_lid(result):
            parent_lid = tracker.get_lid(adata) if adata else None
            parents = [parent_lid] if parent_lid else []
            desc = func.__name__
            if not parents and args:
                desc = f"{func.__name__}('{args[0]}')"
            tracker.register(
                result,
                parents=parents,
                creation_op=func.__name__,
                description=desc,
            )

        return result

    # Sentinel attribute so track_module() can skip already-wrapped functions
    _tracked_wrapper._lotus_tracked = True
    return _tracked_wrapper


def track_module(module):
    """Auto-wrap every public function in *module* with :func:`tracked`.

    Iterates over ``module.__all__`` (falling back to ``dir(module)``)
    and replaces each public, callable, non-class attribute with a
    ``tracked``-wrapped version.  Already-wrapped functions are skipped.
    """
    names = getattr(module, "__all__", dir(module))
    for name in names:
        if name.startswith("_"):
            continue
        obj = getattr(module, name, None)
        if obj is None or not callable(obj) or isinstance(obj, type):
            continue
        # Skip functions that are already wrapped
        if getattr(obj, "_lotus_tracked", False):
            continue
        setattr(module, name, tracked(obj))
