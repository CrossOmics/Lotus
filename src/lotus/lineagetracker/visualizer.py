"""Build a lineage DAG from the lineage JSON and render it as a flowchart PNG
using the ``diagrams`` library (backed by Graphviz)."""

from __future__ import annotations

from html import escape
import json
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any, Literal
from loguru import logger

from diagrams import Diagram, Edge, Node

from .models import LineageNode

# ── colour constants ─────────────────────────────────────────────
_ADATA_COLOR = "#BAE6FD"   # Light blue  — all AnnData data nodes
_OP_COLOR    = "#FDE68A"   # Light amber — all operation nodes
_OUTPUT_EDGE_COLOR = "#111111"
_PARENT_EDGE_COLOR = "#666666"
_OP_EDGE_COLOR = "#4B5563"
_DATA_TEXT_COLOR = "#0F172A"
_DATA_META_COLOR = "#475569"
_INPUT_COLOR = "#DCFCE7"
_INPUT_BORDER_COLOR = "#166534"

NodeRole = Literal["input", "intermediate", "output"]


def _graphviz_available() -> bool:
    """Return ``True`` when the Graphviz ``dot`` executable is available."""
    return shutil.which("dot") is not None


def _log_missing_graphviz(output_file: Path) -> None:
    """Emit a clear warning when Graphviz is unavailable."""
    logger.warning(
        "Skipping lineage PNG render for {} because Graphviz 'dot' was not found on PATH. "
        "Install Graphviz with your system package manager "
        "(for example: 'brew install graphviz' on macOS, 'sudo apt-get install graphviz' on Ubuntu/Debian, "
        "or 'winget install Graphviz.Graphviz' on Windows).",
        output_file,
    )


def _format_arg_value(value: Any, max_len: int = 28) -> str:
    """Render a compact, bounded string for an operation argument value."""
    if isinstance(value, str):
        rendered = value
    else:
        rendered = repr(value)
    if len(rendered) > max_len:
        return rendered[: max_len - 3] + "..."
    return rendered


def _compact_text(value: str, max_len: int = 24) -> str:
    """Render a bounded label-friendly string."""
    if len(value) > max_len:
        return value[: max_len - 3] + "..."
    return value


def _node_role(node: LineageNode, has_children: bool) -> NodeRole:
    """Classify a lineage node for workflow-oriented styling."""
    if not node.parents:
        return "input"
    if not has_children:
        return "output"
    return "intermediate"


def _data_op_label(node: LineageNode) -> str:
    """Return the operation context line for an AnnData node."""
    if not node.parents:
        return "root"
    return (node.creation_op or "process").lower()


def _data_node_label(node: LineageNode, lid: str, role: NodeRole) -> str:
    """Build a workflow-style HTML label for an AnnData lineage node."""
    display_name = escape(_compact_text(node.display_name or node.description or lid[:8]))
    role_title = {
        "input": "input",
        "intermediate": "anndata",
        "output": "output",
    }[role]
    op_label = escape(_data_op_label(node))
    return """<
<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="2">
  <TR><TD ALIGN="CENTER"><FONT POINT-SIZE="9" COLOR="#1E3A5F">{role_title}</FONT></TD></TR>
  <TR><TD ALIGN="CENTER"><B>{display_name}</B></TD></TR>
  <TR><TD ALIGN="CENTER"><FONT POINT-SIZE="10" COLOR="{text_color}">AnnData  {n_obs} x {n_vars}</FONT></TD></TR>
  <TR><TD ALIGN="CENTER"><FONT POINT-SIZE="8" COLOR="{meta_color}">op: {op_label}</FONT></TD></TR>
</TABLE>
>""".format(
        role_title=role_title,
        display_name=display_name,
        n_obs=node.shape[0],
        n_vars=node.shape[1],
        op_label=op_label,
        text_color=_DATA_TEXT_COLOR,
        meta_color=_DATA_META_COLOR,
    )


def _data_node_attrs(node: LineageNode, role: NodeRole) -> dict[str, str]:
    """Return Graphviz attributes for data nodes."""
    if role == "input":
        return {
            "shape": "box",
            "style": "rounded,filled,bold",
            "fillcolor": _INPUT_COLOR,
            "color": _INPUT_BORDER_COLOR,
            "penwidth": "2.4",
            "margin": "0.12,0.09",
        }
    if role == "output":
        return {
            "shape": "box",
            "style": "rounded,filled",
            "fillcolor": _ADATA_COLOR,
            "color": "#0F172A",
            "penwidth": "2.0",
            "margin": "0.12,0.08",
        }
    return {
        "shape": "box",
        "style": "rounded,filled",
        "fillcolor": _ADATA_COLOR,
        "color": "#2563EB",
        "penwidth": "1.5",
        "margin": "0.12,0.08",
    }


def _operation_node_label(method: str, args: dict[str, Any], index: int) -> str:
    """Build an HTML label for a single operation node."""
    short_name = method.split(".")[-1]
    arg_rows: list[str] = []
    if args:
        items = list(args.items())
        for key, value in items[:2]:
            arg_rows.append(
                "<TR><TD ALIGN=\"LEFT\"><FONT POINT-SIZE=\"9\" COLOR=\"#6B7280\">{}</FONT></TD></TR>".format(
                    escape(f"{key}={_format_arg_value(value, max_len=20)}")
                )
            )
        if len(items) > 2:
            arg_rows.append(
                "<TR><TD ALIGN=\"LEFT\"><FONT POINT-SIZE=\"8\" COLOR=\"#94A3B8\">+{} more</FONT></TD></TR>".format(
                    len(items) - 2
                )
            )
    else:
        arg_rows.append(
            "<TR><TD ALIGN=\"LEFT\"><FONT POINT-SIZE=\"9\" COLOR=\"#6B7280\">no args</FONT></TD></TR>"
        )
    return """<
<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="2">
    <TR><TD ALIGN="CENTER"><B>{method}</B></TD></TR>
  {arg_rows}
</TABLE>
>""".format(method=escape(short_name), arg_rows="".join(arg_rows))


def _parent_edge_attrs(role: NodeRole) -> dict[str, str]:
    """Return Graphviz edge styling for parent-child lineage edges."""
    if role == "output":
        return {"color": _OUTPUT_EDGE_COLOR, "style": "solid", "penwidth": "2.4"}
    return {"color": _PARENT_EDGE_COLOR, "style": "solid", "penwidth": "1.4"}


def _child_lookup(nodes: dict[str, LineageNode]) -> dict[str, set[str]]:
    """Build a reverse index of lineage node children."""
    children = {lid: set() for lid in nodes}
    for child_lid, node in nodes.items():
        for parent_lid in node.parents:
            children.setdefault(parent_lid, set()).add(child_lid)
    return children


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

    if not _graphviz_available():
        _log_missing_graphviz(output_file)
        return

    raw = json.loads(lineage_file.read_text(encoding="utf-8"))
    if not raw:
        return

    nodes = {
        lid: LineageNode.model_validate(data) for lid, data in raw.items()
    }

    if not nodes:
        return

    child_lookup = _child_lookup(nodes)
    node_roles = {
        lid: _node_role(node, has_children=bool(child_lookup.get(lid)))
        for lid, node in nodes.items()
    }


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
        "margin": "0.12,0.08",
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
            role = node_roles[lid]
            label = _data_node_label(node, lid, role)
            diagram_nodes[lid] = Node(
                label,
                **_data_node_attrs(node, role),
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
                        color="#7C5E10",
                        penwidth="1.4",
                        margin="0.18,0.12",
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
                    anchor >> Edge(**_parent_edge_attrs(node_roles[lid])) >> diagram_nodes[lid]

            # Render operation history as chained nodes from each data node.
            chain = operation_nodes.get(lid, [])
            if chain:
                diagram_nodes[lid] >> Edge(color=_OP_EDGE_COLOR, style="dashed") >> chain[0]
                for prev, curr in zip(chain, chain[1:]):
                    prev >> Edge(color=_OP_EDGE_COLOR, style="dashed") >> curr


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
    if not _graphviz_available():
        _log_missing_graphviz(root_dir / "lineage_graph.png")
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
