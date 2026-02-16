"""Central lineage tracker — singleton that owns the DAG and persists it to JSON."""

from __future__ import annotations

import json
import uuid
from pathlib import Path
from typing import Any

import numpy as np
from anndata import AnnData

from .models import LineageNode, OperationRecord


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
        return f"<AnnData:{lid[:8] if lid else 'untracked'}>"
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
        self._data_dir = Path.cwd() / "lineage-tracer_data"
        self._data_dir.mkdir(exist_ok=True)
        self._lineage_file = self._data_dir / "lineage.json"
        self._graph_file = self._data_dir / "lineage_graph.png"

        # lid → node: the authoritative DAG (includes nodes from prior sessions)
        self._nodes: dict[str, LineageNode] = {}
        # id(adata) → lid: session-local mapping that also covers AnnData views
        # (views share their parent's .uns, so we can't rely on uns alone)
        self._lid_map: dict[int, str] = {}

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
        """Restore the DAG from the JSON file if it exists."""
        if self._lineage_file.exists():
            try:
                raw = json.loads(self._lineage_file.read_text(encoding="utf-8"))
                for lid, node_data in raw.items():
                    self._nodes[lid] = LineageNode.model_validate(node_data)
            except (json.JSONDecodeError, Exception):
                pass  # Start fresh if the file is corrupted

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
    def lineage_file(self) -> Path:
        return self._lineage_file

    @property
    def graph_file(self) -> Path:
        return self._graph_file

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

    def register(
        self,
        adata: AnnData,
        *,
        parents: list[str],
        creation_op: str | None,
        description: str = "",
    ) -> str:
        """Create a new DAG node for *adata* and persist the change.

        For non-view objects the lid is also written into
        ``adata.uns['_lotus_lid']`` so it survives h5ad round-trips.
        """
        lid = str(uuid.uuid4())
        node = LineageNode(
            lid=lid,
            parents=parents,
            creation_op=creation_op,
            description=description,
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

    def record_op(self, adata: AnnData, method: str, args: dict):
        """Append an in-place operation record to *adata*'s node."""
        lid = self.ensure_registered(adata)
        record = OperationRecord(
            method=method,
            args=serialize_args(args, self),
        )
        self._nodes[lid].operations.append(record)
        self.save()
