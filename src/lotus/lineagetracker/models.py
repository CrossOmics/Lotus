"""Pydantic models for lineage tracking records."""

from __future__ import annotations

import uuid
from datetime import datetime

from pydantic import BaseModel, Field


class OperationRecord(BaseModel):
    """A single in-place operation performed on an AnnData object.

    Recorded every time a tracked function mutates an existing AnnData
    (e.g. normalize_total, log1p, scale).
    """

    method: str  # Fully-qualified method name, e.g. "preprocessing.normalize_total"
    args: dict = {}  # Serialised arguments (AnnData refs replaced with lid stubs)
    timestamp: datetime = Field(default_factory=datetime.now)


class LineageNode(BaseModel):
    """One node in the lineage DAG, representing a single AnnData object.

    Root nodes (loaded from disk or created externally) have an empty
    ``parents`` list.  Derived nodes carry one parent (slice / copy) or
    multiple parents (concat).
    """

    lid: str = Field(default_factory=lambda: str(uuid.uuid4()))
    parents: list[str] = []  # Lineage IDs of parent AnnData objects
    created_at: datetime = Field(default_factory=datetime.now)
    creation_op: str | None = None  # "slice", "copy", "concat", or func name
    description: str = ""  # Human-readable origin summary
    variable_name: str | None = None  # Latest observed Python variable name
    display_name: str | None = None  # Preferred label for graph rendering
    shape: tuple[int, int] = (0, 0)  # (n_obs, n_vars) at registration time
    operations: list[OperationRecord] = []  # Chronological list of in-place ops
