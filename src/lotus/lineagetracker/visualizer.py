"""Build a lineage DAG from the lineage JSON and render it as a flowchart PNG
using the ``diagrams`` library (backed by Graphviz)."""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Any
from loguru import logger

from diagrams import Diagram, Edge, Node

from .models import LineageNode

# ── colour constants ─────────────────────────────────────────────
_ADATA_COLOR = "#BAE6FD"   # Light blue  — all AnnData data nodes
_OP_COLOR    = "#FDE68A"   # Light amber — all operation nodes


def _format_arg_value(value: Any, max_len: int = 28) -> str:
    """Render a compact, bounded string for an operation argument value."""
    if isinstance(value, str):
        rendered = value
    else:
        rendered = repr(value)
    if len(rendered) > max_len:
        return rendered[: max_len - 3] + "..."
    return rendered


def _data_node_label(node: LineageNode, lid: str) -> str:
    """Build a compact label for an AnnData lineage node."""
    display_name = node.display_name or node.description or lid[:8]
    if not node.parents:
        op_label = "root"
    else:
        op_label = (node.creation_op or "process").lower()
    return "\n".join(
        [
            display_name,
            f"{node.shape[0]} x {node.shape[1]}",
            f"op: {op_label}",
        ]
    )


def _operation_node_label(method: str, args: dict[str, Any], index: int) -> str:
    """Build a compact label for a single operation node."""
    short_name = method.split(".")[-1]
    lines = [f"{index}. {short_name}"]
    if args:
        items = list(args.items())
        for key, value in items[:3]:
            lines.append(f"{key}={_format_arg_value(value)}")
    else:
        lines.append("no args")
    return "\n".join(lines)


def _parent_anchor_index(parent: LineageNode, child: LineageNode) -> int:
    """Return the index of the latest parent operation that happened before
    child creation time. Returns -1 when the edge should start from parent data."""

    def find_less_equal_bound(nums: list[datetime], target: datetime) -> int:
        """Find the index of the last value that is <= target."""
        l, r = 0, len(nums) - 1
        while l <= r:
            mid = (l + r) // 2
            if nums[mid] <= target:
                l = mid + 1
            else:
                r = mid - 1
        return r

    if not parent.operations:
        return -1
    return find_less_equal_bound([op.timestamp for op in parent.operations], child.created_at)


def visualize(lineage_file: Path, output_file: Path):
    """Read the lineage JSON file, build a DAG, and save as a flowchart PNG.

    All AnnData nodes are rendered as filled boxes in one color.
    All operation nodes are rendered as filled boxes in another color.
    Directed arrows point from parent to child.
    """
    if not lineage_file.exists():
        return

    raw = json.loads(lineage_file.read_text(encoding="utf-8"))
    if not raw:
        return

    nodes = {
        lid: LineageNode.model_validate(data) for lid, data in raw.items()
    }

    if not nodes:
        return


    # Strip the extension — Diagram appends it automatically.
    out_stem = str(output_file.with_suffix(""))

    graph_attr = {
        "labelloc": "t",
        "labeljust": "c",
        "fontsize": "14",
        "bgcolor": "white",
        "pad": "0.5",
        "nodesep": "1.0",
        "ranksep": "1.2",
    }
    node_attr = {
        "fixedsize": "false",
        "labelloc": "c",
        "fontname": "Sans-Serif",
        "fontsize": "11",
    }

    with Diagram(
        filename=out_stem,
        outformat="png",
        show=False,
        direction="TB",
        graph_attr=graph_attr,
        node_attr=node_attr,
    ):
        diagram_nodes: dict[str, object] = {}
        operation_nodes: dict[str, list[object]] = {}

        for lid, node in nodes.items():
            label = _data_node_label(node, lid)
            diagram_nodes[lid] = Node(
                label,
                shape="box",
                style="rounded,filled",
                fillcolor=_ADATA_COLOR,
            )

            op_chain: list[object] = []
            for idx, op in enumerate(node.operations, start=1):
                op_label = _operation_node_label(op.method, op.args, idx)
                op_chain.append(
                    Node(
                        op_label,
                        shape="box",
                        style="rounded,filled",
                        fillcolor=_OP_COLOR,
                    )
                )
            operation_nodes[lid] = op_chain

        for lid, node in nodes.items():
            for parent_lid in node.parents:
                parent_node = nodes.get(parent_lid)
                if parent_node is None:
                    continue
                anchor_index = _parent_anchor_index(parent_node, node)
                parent_chain = operation_nodes.get(parent_lid, [])
                if 0 <= anchor_index < len(parent_chain):
                    anchor = parent_chain[anchor_index]
                else:
                    anchor = diagram_nodes.get(parent_lid)
                if anchor is not None:
                    anchor >> Edge(color="#666666") >> diagram_nodes[lid]

            # Render operation history as chained nodes from each data node.
            chain = operation_nodes.get(lid, [])
            if chain:
                diagram_nodes[lid] >> Edge(color="#4B5563", style="dashed") >> chain[0]
                for prev, curr in zip(chain, chain[1:]):
                    prev >> Edge(color="#4B5563", style="dashed") >> curr


def render_missing_pngs(root_dir: Path) -> None:
    """Scan all subfolders of *root_dir* and render a PNG for any that have a
    ``lineage.json`` but no ``lineage_graph.png``.

    Rules:
    - Subfolder has ``lineage.json`` but **no** ``lineage_graph.png`` → render
    - Subfolder is empty (neither file) → skip
    - Subfolder already has both files → skip
    """
    if not root_dir.is_dir():
        return
    for subfolder in sorted(root_dir.iterdir()):
        if not subfolder.is_dir():
            continue
        json_file = subfolder / "lineage.json"
        png_file  = subfolder / "lineage_graph.png"
        if json_file.exists() and not png_file.exists():
            try:
                visualize(json_file, png_file)
            except Exception:
                logger.exception("error visualizing png for {}", json_file)
