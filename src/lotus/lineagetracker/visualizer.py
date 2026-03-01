"""Build a lineage DAG from the lineage JSON and render it as a flowchart PNG
using the ``diagrams`` library (backed by Graphviz)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from diagrams import Diagram, Edge, Node

from .models import LineageNode

_PREPROCESSING_OPS = {
    "calculate_qc_metrics",
    "filter_cells",
    "filter_genes",
    "normalize_total",
    "log1p",
    "scale",
    "highly_variable_genes",
    "regress_out",
    "combat",
    "neighbors",
    "pca",
}
_IO_OPS = {
    "read",
    "write",
    "standardize_load",
    "read_h5ad",
    "read_10x_h5",
    "read_10x_mtx",
    "read_visium",
    "read_csv",
    "read_text",
    "read_excel",
    "read_hdf",
    "read_mtx",
    "read_loom",
    "read_umi_tools",
}
_CLUSTERING_OPS = {
    "leiden",
    "louvain",
    "kmeans",
    "cluster_core",
}
_CORE_ANALYSIS_OPS = {
    "flowrank",
    "stable_core",
    "fine_grained_core",
    "propagation_from_core",
    "core_selection",
}
_OPERATION_COLOR_MAP = {
    "preprocessing": "#93C5FD",
    "io": "#86EFAC",
    "clustering": "#FCD34D",
    "visualization": "#C4B5FD",
    "core_analysis": "#FCA5A5",
    "other": "#D1D5DB",
}
_OPERATION_LEGEND_ITEMS = [
    ("preprocessing", _OPERATION_COLOR_MAP["preprocessing"]),
    ("io", _OPERATION_COLOR_MAP["io"]),
    ("clustering", _OPERATION_COLOR_MAP["clustering"]),
    ("visualization", _OPERATION_COLOR_MAP["visualization"]),
    ("core_analysis", _OPERATION_COLOR_MAP["core_analysis"]),
    ("other", _OPERATION_COLOR_MAP["other"]),
]


def _format_arg_value(value: Any, max_len: int = 28) -> str:
    """Render a compact, bounded string for an operation argument value."""
    if isinstance(value, str):
        rendered = value
    else:
        rendered = repr(value)
    if len(rendered) > max_len:
        return rendered[: max_len - 3] + "..."
    return rendered


def _format_operation_lines(
    method: str,
    args: dict[str, Any],
    max_args: int = 3,
) -> list[str]:
    """Format one operation as a compact multi-line call-like block."""
    short_name = method.split(".")[-1]
    if not args:
        return [f"- {short_name}()"]

    items = list(args.items())
    lines = [f"- {short_name}("]
    for key, value in items[:max_args]:
        lines.append(f"  {key}={_format_arg_value(value)}")
    if len(items) > max_args:
        lines.append(f"  ...+{len(items) - max_args} args")
    lines.append(")")
    return lines


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
        # if len(items) > 3:
        #     lines.append(f"...+{len(items) - 3} args")
    else:
        lines.append("no args")
    return "\n".join(lines)


def _operation_category(method: str) -> str:
    """Return operation category key for a method string."""
    short_name = method.split(".")[-1]

    if method.startswith("preprocessing.") or short_name in _PREPROCESSING_OPS:
        return "preprocessing"
    if method.startswith("io.") or short_name in _IO_OPS:
        return "io"
    if method.startswith("clustering.") or short_name in _CLUSTERING_OPS:
        return "clustering"
    if method.startswith("visualization."):
        return "visualization"
    if (
        method.startswith("core_analysis.")
        or "CorespectModel." in method
        or short_name in _CORE_ANALYSIS_OPS
    ):
        return "core_analysis"
    return "other"


def _operation_fillcolor(method: str) -> str:
    """Return operation node fill color by operation category."""
    return _OPERATION_COLOR_MAP[_operation_category(method)]


def _make_node(node: LineageNode, label: str):
    """Pick a Graphviz shape based on how the node was created."""
    op = (node.creation_op or "").lower()
    shape = "box"
    if not node.parents:
        shape = "cylinder"
    elif op == "slice":
        shape = "note"
    elif op == "copy":
        shape = "box3d"
    elif op == "concat":
        shape = "diamond"
    return Node(label, shape=shape)


def _make_operation_node(label: str, method: str):
    """Create a visual node for one operation entry."""
    return Node(
        label,
        shape="box",
        style="rounded,filled",
        fillcolor=_operation_fillcolor(method),
    )


def _add_operation_legend(anchor_node: object | None = None):
    """Add a compact legend explaining operation-node colors."""
    legend_nodes = [
        Node(
            f"legend: {label}",
            shape="box",
            style="rounded,filled",
            fillcolor=color,
        )
        for label, color in _OPERATION_LEGEND_ITEMS
    ]
    if anchor_node is not None and legend_nodes:
        anchor_node >> Edge(style="invis") >> legend_nodes[0]
    for prev, curr in zip(legend_nodes, legend_nodes[1:]):
        prev >> Edge(style="invis") >> curr


def visualize(lineage_file: Path, output_file: Path):
    """Read the lineage JSON file, build a DAG, and save as a flowchart PNG.

    Each AnnData node shows: display name and shape.
    Operation history is rendered as dedicated operation nodes.
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
        # "label": "Lotus Lineage DAG",
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
        # "Lotus Lineage DAG",
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
            diagram_nodes[lid] = _make_node(node, label)

            op_chain: list[object] = []
            for idx, op in enumerate(node.operations, start=1):
                op_label = _operation_node_label(op.method, op.args, idx)
                op_chain.append(_make_operation_node(op_label, op.method))
            operation_nodes[lid] = op_chain

        # legend_anchor = next(iter(diagram_nodes.values()), None)
        # _add_operation_legend(legend_anchor)

        for lid, node in nodes.items():
            for parent_lid in node.parents:
                if parent_lid in diagram_nodes:
                    diagram_nodes[parent_lid] >> Edge(color="#666666") >> diagram_nodes[lid]

            # Render operation history as chained nodes from each data node.
            chain = operation_nodes.get(lid, [])
            if chain:
                diagram_nodes[lid] >> Edge(color="#4B5563", style="dashed") >> chain[0]
                for prev, curr in zip(chain, chain[1:]):
                    prev >> Edge(color="#4B5563", style="dashed") >> curr
