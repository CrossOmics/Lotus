"""Central lineage tracker — singleton that owns the DAG and persists it to JSON."""

from __future__ import annotations

import ast
import inspect
import json
import textwrap
import uuid
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
from anndata import AnnData

from .models import LineageNode, OperationRecord

_INTERNAL_DIR = Path(__file__).resolve().parent
_WORKSPACE_DIR = Path.cwd().resolve()
_GENERIC_LOCAL_NAMES = {"result"}


def _is_internal_frame(filename: str) -> bool:
    try:
        return Path(filename).resolve().is_relative_to(_INTERNAL_DIR)
    except (RuntimeError, OSError, ValueError):
        return False


def _is_workspace_frame(filename: str) -> bool:
    try:
        return Path(filename).resolve().is_relative_to(_WORKSPACE_DIR)
    except (RuntimeError, OSError, ValueError):
        return False


def _is_valid_inferred_name(name: str) -> bool:
    if not name or name.startswith("_"):
        return False
    if name in _GENERIC_LOCAL_NAMES:
        return False
    return True


# ── argument serialisation helpers ───────────────────────────────


def _serialize_value(v: Any, tracker: LineageTracker) -> Any:
    """Best-effort conversion of a single value to a JSON-safe type.

    * AnnData → ``"<AnnData:lid_prefix>"``
    * numpy arrays / scalars → shape description or plain Python scalar
    * Path → string
    * Containers → recursed
    * Anything else → ``repr()``
    """
    if isinstance(v, AnnData):
        lid = tracker.get_lid(v)
        if lid and lid in tracker._nodes:
            name = tracker._nodes[lid].display_name or lid[:8]
        else:
            name = "untracked"
        return f"<AnnData:{name}>"
    if isinstance(v, np.ndarray):
        return f"<ndarray shape={v.shape} dtype={v.dtype}>"
    if isinstance(v, np.generic):
        return v.item()
    if isinstance(v, Path):
        return str(v)
    if isinstance(v, (str, int, float, bool, type(None))):
        return v
    if isinstance(v, (list, tuple)):
        return [_serialize_value(i, tracker) for i in v]
    if isinstance(v, dict):
        return {str(k): _serialize_value(val, tracker) for k, val in v.items()}
    return repr(v)


def serialize_args(args: dict, tracker: LineageTracker) -> dict:
    """Serialise a full keyword-argument dict for JSON storage."""
    return {k: _serialize_value(v, tracker) for k, v in args.items()}


# ── tracker singleton ────────────────────────────────────────────


class LineageTracker:
    """Session-scoped singleton that maintains a DAG of AnnData lineage.

    * Each AnnData is assigned a UUID (``lid``) on first contact.
    * In-place operations are appended to the node's ``operations`` list.
    * Derivations (slice / copy / concat) create new child nodes with
      parent edges.
    * The full DAG is persisted to ``<cwd>/lineage-tracer_data/lineage.json``
      after every mutation so it survives crashes.
    """

    _instance: LineageTracker | None = None

    def __init__(self):
        self._root_dir = Path.cwd() / "lineage-tracer_data"
        self._root_dir.mkdir(exist_ok=True)
        self._session_ts = (
            datetime.now().strftime("%Y%m%d_%H%M%S") + "_" + uuid.uuid4().hex
        )
        self._data_dir = self._root_dir / self._session_ts
        self._data_dir.mkdir(exist_ok=True)
        self._lineage_file = self._data_dir / "lineage.json"
        self._graph_file = self._data_dir / "lineage_graph.png"

        # lid → node: the authoritative DAG (includes nodes from prior sessions)
        self._nodes: dict[str, LineageNode] = {}
        # id(adata) → lid: session-local mapping that also covers AnnData views
        # (views share their parent's .uns, so we can't rely on uns alone)
        self._lid_map: dict[int, str] = {}
        # Depth counter for nested tracked calls. Depth>0 means we are executing
        # inside a top-level tracked operation and should suppress nested records.
        self._operation_depth: int = 0

        self._load()

    # ── singleton access ─────────────────────────────────────────

    @classmethod
    def instance(cls) -> LineageTracker:
        """Return the global tracker, creating it on first call."""
        if cls._instance is None:
            cls._instance = cls()
        return cls._instance

    @classmethod
    def reset(cls):
        """Discard the singleton (useful in tests)."""
        cls._instance = None

    # ── persistence ──────────────────────────────────────────────

    def _load(self):
        """Each session starts with a fresh DAG (its own timestamped folder)."""
        pass

    def save(self):
        """Write the full DAG to disk (overwrite)."""
        data = {
            lid: node.model_dump(mode="json")
            for lid, node in self._nodes.items()
        }
        self._lineage_file.write_text(
            json.dumps(data, indent=2, ensure_ascii=False, default=str),
            encoding="utf-8",
        )

    @property
    def root_dir(self) -> Path:
        return self._root_dir

    @property
    def lineage_file(self) -> Path:
        return self._lineage_file

    @property
    def graph_file(self) -> Path:
        return self._graph_file

    @property
    def is_nested_operation(self) -> bool:
        return self._operation_depth > 0

    @contextmanager
    def operation_scope(self):
        self._operation_depth += 1
        try:
            yield
        finally:
            self._operation_depth -= 1

    # ── lid helpers ──────────────────────────────────────────────

    def get_lid(self, adata: AnnData) -> str | None:
        """Look up the lineage ID for *adata*.

        Checks the session-local ``id()`` map first (works for views),
        then falls back to ``adata.uns['_lotus_lid']`` (persisted in h5ad).
        """
        lid = self._lid_map.get(id(adata))
        if lid:
            return lid
        return adata.uns.get("_lotus_lid")

    def ensure_registered(self, adata: AnnData) -> str:
        """Return *adata*'s lid, auto-registering it as an external root
        node if it has never been seen before."""
        lid = self.get_lid(adata)
        if lid is None:
            lid = self.register(
                adata, parents=[], creation_op=None, description="external"
            )
        return lid

    # ── core operations ──────────────────────────────────────────

    @staticmethod
    def _normalize_name(value: Any) -> str | None:
        if value is None:
            return None
        name = str(value).strip()
        return name or None

    @classmethod
    def _resolve_display_name(
        cls,
        adata: AnnData,
        description: str,
        variable_name: str | None = None,
    ) -> str:
        """Name priority: variable name > ``adata.uns['name']`` > description."""
        resolved_variable_name = cls._normalize_name(variable_name)
        if resolved_variable_name:
            return resolved_variable_name

        uns_name = cls._normalize_name(adata.uns.get("name"))
        if uns_name:
            return uns_name

        return description

    @staticmethod
    def _iter_target_names(target: ast.AST):
        if isinstance(target, ast.Name):
            yield target.id
            return
        if isinstance(target, (ast.Tuple, ast.List)):
            for elt in target.elts:
                yield from LineageTracker._iter_target_names(elt)

    @staticmethod
    def _parse_frame_tree(frameinfo: inspect.FrameInfo) -> ast.Module | None:
        if not frameinfo.code_context:
            return None
        source = textwrap.dedent("".join(frameinfo.code_context)).strip()
        if not source:
            return None
        try:
            parsed = ast.parse(source)
        except SyntaxError:
            return None
        if isinstance(parsed, ast.Module):
            return parsed
        return None

    @classmethod
    def _extract_name_from_call_args(
        cls,
        frameinfo: inspect.FrameInfo,
        adata: AnnData,
    ) -> str | None:
        tree = cls._parse_frame_tree(frameinfo)
        if tree is None:
            return None

        local_scope = frameinfo.frame.f_locals
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            for arg in node.args:
                if (
                    isinstance(arg, ast.Name)
                    and _is_valid_inferred_name(arg.id)
                    and local_scope.get(arg.id) is adata
                ):
                    return arg.id
            for kw in node.keywords:
                if (
                    kw.arg
                    and isinstance(kw.value, ast.Name)
                    and _is_valid_inferred_name(kw.value.id)
                    and local_scope.get(kw.value.id) is adata
                ):
                    return kw.value.id
        return None

    @staticmethod
    def _is_register_call(call: ast.Call) -> bool:
        if isinstance(call.func, ast.Attribute):
            return call.func.attr == "register"
        if isinstance(call.func, ast.Name):
            return call.func.id == "register"
        return False

    @classmethod
    def _extract_name_from_assignment(
        cls,
        frameinfo: inspect.FrameInfo,
    ) -> str | None:
        tree = cls._parse_frame_tree(frameinfo)
        if tree is None:
            return None

        for stmt in tree.body:
            if isinstance(stmt, ast.Assign):
                if isinstance(stmt.value, ast.Call) and cls._is_register_call(stmt.value):
                    continue
                for target in stmt.targets:
                    for name in cls._iter_target_names(target):
                        if not _is_valid_inferred_name(name):
                            continue
                        return name
            if isinstance(stmt, ast.AnnAssign):
                if isinstance(stmt.value, ast.Call) and cls._is_register_call(stmt.value):
                    continue
                for name in cls._iter_target_names(stmt.target):
                    if not _is_valid_inferred_name(name):
                        continue
                    return name
        return None

    @staticmethod
    def _extract_name_from_locals(
        frameinfo: inspect.FrameInfo,
        adata: AnnData,
    ) -> str | None:
        for name, value in frameinfo.frame.f_locals.items():
            if not _is_valid_inferred_name(name):
                continue
            if value is adata:
                return name
        return None

    def _detect_variable_name(self, adata: AnnData) -> str | None:
        """Best-effort variable-name detection from non-lineagetracker frames."""
        stack = inspect.stack(context=1)
        try:
            for frameinfo in stack[2:]:
                if _is_internal_frame(frameinfo.filename):
                    continue
                if not _is_workspace_frame(frameinfo.filename):
                    continue

                from_call_args = self._extract_name_from_call_args(frameinfo, adata)
                if from_call_args:
                    return from_call_args

                from_assignment = self._extract_name_from_assignment(frameinfo)
                if from_assignment:
                    return from_assignment

                from_locals = self._extract_name_from_locals(frameinfo, adata)
                if from_locals:
                    return from_locals
        finally:
            del stack
        return None

    def register(
        self,
        adata: AnnData,
        *,
        parents: list[str],
        creation_op: str | None,
        description: str = "",
        variable_name: str | None = None,
    ) -> str:
        """Create a new DAG node for *adata* and persist the change.

        For non-view objects the lid is also written into
        ``adata.uns['_lotus_lid']`` so it survives h5ad round-trips.
        """
        lid = str(uuid.uuid4())
        resolved_variable_name = (
            self._normalize_name(variable_name) or self._detect_variable_name(adata)
        )
        node = LineageNode(
            lid=lid,
            parents=parents,
            creation_op=creation_op,
            description=description,
            variable_name=resolved_variable_name,
            display_name=self._resolve_display_name(
                adata,
                description,
                resolved_variable_name,
            ),
            shape=(adata.n_obs, adata.n_vars),
        )
        self._nodes[lid] = node
        self._lid_map[id(adata)] = lid
        # Views share .uns with their parent — writing here would
        # overwrite the parent's lid, so skip views.
        if not adata.is_view:
            adata.uns["_lotus_lid"] = lid
        self.save()
        return lid

    def bind_variable_name(self, adata: AnnData, variable_name: str) -> str:
        """Explicitly bind or replace the variable name associated with *adata*."""
        normalized = self._normalize_name(variable_name)
        if normalized is None:
            raise ValueError("variable_name must be a non-empty string.")

        lid = self.ensure_registered(adata)
        node = self._nodes[lid]
        node.variable_name = normalized
        node.display_name = self._resolve_display_name(
            adata, node.description, node.variable_name
        )
        self.save()
        return lid

    def record_op(
        self,
        adata: AnnData,
        method: str,
        args: dict,
        *,
        variable_name: str | None = None,
    ):
        """Append an in-place operation record to *adata*'s node."""
        lid = self.ensure_registered(adata)
        node = self._nodes[lid]
        observed_variable_name = (
            self._normalize_name(variable_name) or self._detect_variable_name(adata)
        )
        if observed_variable_name:
            node.variable_name = observed_variable_name
        node.display_name = self._resolve_display_name(
            adata, node.description, node.variable_name
        )
        record = OperationRecord(
            method=method,
            args=serialize_args(args, self),
        )
        node.operations.append(record)
        self.save()


def bind_variable_name(adata: AnnData, variable_name: str) -> str:
    """Explicit helper API for binding an AnnData variable name to its lineage ID."""
    return LineageTracker.instance().bind_variable_name(adata, variable_name)
