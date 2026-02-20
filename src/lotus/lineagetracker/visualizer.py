"""Build a networkx DAG from the lineage JSON and render it as a PNG."""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # Non-interactive backend — safe for headless servers
import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402

from .models import LineageNode  # noqa: E402


def _topological_layout(G: nx.DiGraph) -> dict:
    """Assign (x, y) positions using a top-down layered layout.

    Nodes are placed in rows by their topological depth (root = top).
    Siblings within the same row are centred horizontally.
    """
    # Compute the depth (level) of every node
    levels: dict[str, int] = {}
    for node in nx.topological_sort(G):
        parents = list(G.predecessors(node))
        if not parents:
            levels[node] = 0
        else:
            levels[node] = max(levels[p] for p in parents) + 1

    # Group nodes by level
    buckets: dict[int, list[str]] = defaultdict(list)
    for node, lvl in levels.items():
        buckets[lvl].append(node)

    # Assign coordinates — centre each row around x=0
    pos = {}
    for lvl, nodes in buckets.items():
        for i, node in enumerate(nodes):
            x = (i - (len(nodes) - 1) / 2) * 4
            y = -lvl * 2.5
            pos[node] = (x, y)
    return pos


def visualize(lineage_file: Path, output_file: Path):
    """Read the lineage JSON file, build a DAG, and save as PNG.

    Each node box shows: display name, shape, and operation
    history.  Directed arrows point from parent to child.
    """
    if not lineage_file.exists():
        return

    raw = json.loads(lineage_file.read_text(encoding="utf-8"))
    if not raw:
        return

    nodes = {
        lid: LineageNode.model_validate(data) for lid, data in raw.items()
    }

    # Build the directed graph and human-readable labels
    G = nx.DiGraph()
    labels: dict[str, str] = {}

    for lid, node in nodes.items():
        display_name = node.display_name or node.description or lid[:8]
        parts = [display_name]

        if node.description and node.description != display_name:
            parts.append(node.description)

        parts.append(f"{node.shape[0]} x {node.shape[1]}")

        if node.operations:
            op_names = [op.method.split(".")[-1] for op in node.operations]
            parts.append("ops: " + " -> ".join(op_names))

        G.add_node(lid)
        labels[lid] = "\n".join(parts)

        for parent_lid in node.parents:
            if parent_lid in nodes:
                G.add_edge(parent_lid, lid)

    if len(G.nodes) == 0:
        return

    pos = _topological_layout(G)

    # Scale the figure to the number of nodes
    width = max(10, len(G.nodes) * 2.5)
    height = max(6, len(G.nodes) * 1.2)
    fig, ax = plt.subplots(figsize=(width, height))

    # Invisible nodes — their size controls how far edges shrink from
    # node centres so that arrowheads sit outside the label boxes.
    bbox_size = 3500
    nx.draw_networkx_nodes(
        G, pos, ax=ax, node_size=bbox_size, alpha=0,
    )

    # Directed edges with arrowheads
    nx.draw_networkx_edges(
        G, pos, ax=ax,
        arrows=True,
        arrowstyle="-|>",
        arrowsize=20,
        edge_color="#666666",
        width=2,
        node_size=bbox_size,
        connectionstyle="arc3,rad=0.0",
    )

    # Node labels rendered on top of edges
    nx.draw_networkx_labels(
        G, pos, labels, ax=ax, font_size=8, font_family="monospace",
        bbox=dict(
            boxstyle="round,pad=0.6",
            facecolor="lightblue",
            edgecolor="steelblue",
            alpha=0.9,
        ),
    )

    ax.set_title("Lotus Lineage DAG", fontsize=14, fontweight="bold")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(output_file, dpi=150, bbox_inches="tight")
    plt.close(fig)
