"""Decorators and metaclass that intercept function / method calls to
record lineage operations and snapshot ``adata.X`` into layers."""

import functools
import inspect

import anndata
from anndata import AnnData
from loguru import logger

from .tracker import LineageTracker


# ---- Helper functions
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


def _run_with_tracking(
    *,
    func,
    call_name: str,
    adata: AnnData | None,
    args_for_record: tuple,
    kwargs_for_record: dict,
    invoke,
    description: str,
    before_log: str | None = None,
):
    """Shared tracking flow for both method and function wrappers."""
    tracker = LineageTracker.instance()
    is_nested = tracker.is_nested_operation

    if adata is not None and not is_nested:
        all_args = _build_args_dict(func, args_for_record, kwargs_for_record)
        tracker.record_op(adata, call_name, all_args)

    if before_log is not None:
        logger.info(before_log)
    with tracker.operation_scope():
        result = invoke()

    if (
        not is_nested
        and isinstance(result, AnnData)
        and not tracker.get_lid(result)
    ):
        parent_lid = tracker.get_lid(adata) if adata else None
        tracker.register(
            result,
            parents=[parent_lid] if parent_lid else [],
            creation_op=call_name,
            description=description,
        )
    return result


# ---- Core classes and methods
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
            adata = _extract_adata(*args, **kwargs)
            call_name = f"{cls_name}.{method_name}"
            return _run_with_tracking(
                func=method,
                call_name=call_name,
                adata=adata,
                args_for_record=(self, *args),
                kwargs_for_record=kwargs,
                invoke=lambda: method(self, *args, **kwargs),
                description=call_name,
                before_log=f"[Now in] {call_name}",
            )

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
        adata = _extract_adata(*args, **kwargs)
        description = func.__name__
        if adata is None and args:
            description = f"{func.__name__}('{args[0]}')"
        return _run_with_tracking(
            func=func,
            call_name=func.__name__,
            adata=adata,
            args_for_record=args,
            kwargs_for_record=kwargs,
            invoke=lambda: func(*args, **kwargs),
            description=description,
            before_log=f"[Now in] {func.__name__}"
        )

    return wrapper
