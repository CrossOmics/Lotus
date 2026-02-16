"""
End-to-end test for the lineagetracker package.

Exercises: register, slice, copy, concat, operation recording,
JSON persistence, and PNG visualisation.
"""

import json
import sys
from pathlib import Path

import numpy as np
import anndata

# Ensure the project src is on the path
project_root = Path(__file__).resolve().parent.parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root / "src"))

from lotus.lineagetracker import LineageTracker  # noqa: E402
from lotus.lineagetracker.visualizer import visualize  # noqa: E402


def _make_adata(n_obs: int = 100, n_vars: int = 50) -> anndata.AnnData:
    """Create a small synthetic AnnData for testing."""
    X = np.random.negative_binomial(n=20, p=0.3, size=(n_obs, n_vars)).astype(
        np.float32
    )
    adata = anndata.AnnData(X=X)
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adata.var_names = [f"gene_{i}" for i in range(n_vars)]
    adata.obs["celltype"] = np.random.choice(["T", "B", "NK"], n_obs)
    return adata


def test_lineagetracker():
    # Reset singleton so the test is self-contained
    LineageTracker.reset()
    tracker = LineageTracker.instance()

    # 1. Register root
    adata = _make_adata()
    root_lid = tracker.register(
        adata, parents=[], creation_op="create", description="synthetic 100x50"
    )
    assert tracker.get_lid(adata) == root_lid
    assert root_lid in json.loads(tracker.lineage_file.read_text())
    print(f"[OK] Root registered: {root_lid[:8]}")

    # 2. Record in-place operations
    tracker.record_op(adata, "preprocessing.normalize_total", {"target_sum": 1e4})
    tracker.record_op(adata, "preprocessing.log1p", {})
    node = json.loads(tracker.lineage_file.read_text())[root_lid]
    assert len(node["operations"]) == 2
    assert node["operations"][0]["method"] == "preprocessing.normalize_total"
    print("[OK] Operations recorded on root node")

    # 3. Slice (monkey-patched __getitem__)
    t_cells = adata[adata.obs["celltype"] == "T"]
    t_lid = tracker.get_lid(t_cells)
    assert t_lid is not None
    t_node = json.loads(tracker.lineage_file.read_text())[t_lid]
    assert root_lid in t_node["parents"]
    assert t_node["creation_op"] == "slice"
    print(f"[OK] Slice tracked: {t_lid[:8]}, shape={t_cells.shape}")

    b_cells = adata[adata.obs["celltype"] == "B"]
    b_lid = tracker.get_lid(b_cells)
    assert b_lid is not None
    print(f"[OK] Slice tracked: {b_lid[:8]}, shape={b_cells.shape}")

    # 4. Copy (monkey-patched copy)
    t_copy = t_cells.copy()
    cp_lid = tracker.get_lid(t_copy)
    assert cp_lid is not None
    cp_node = json.loads(tracker.lineage_file.read_text())[cp_lid]
    assert t_lid in cp_node["parents"]
    assert cp_node["creation_op"] == "copy"
    print(f"[OK] Copy tracked: {cp_lid[:8]}")

    # 5. Concat (monkey-patched anndata.concat)
    merged = anndata.concat([t_cells, b_cells])
    m_lid = tracker.get_lid(merged)
    assert m_lid is not None
    m_node = json.loads(tracker.lineage_file.read_text())[m_lid]
    assert set(m_node["parents"]) == {t_lid, b_lid}
    assert m_node["creation_op"] == "concat"
    print(f"[OK] Concat tracked: {m_lid[:8]}, shape={merged.shape}")

    # 6. DAG integrity
    all_nodes = json.loads(tracker.lineage_file.read_text())
    assert len(all_nodes) == 5  # root, t, b, copy, merged
    print(f"[OK] DAG has {len(all_nodes)} nodes")

    # ── 7. Visualisation PNG ─────────────────────────────────────
    visualize(tracker.lineage_file, tracker.graph_file)
    assert tracker.graph_file.exists()
    assert tracker.graph_file.stat().st_size > 0
    print(f"[OK] PNG saved: {tracker.graph_file}")

    # ── Summary ──────────────────────────────────────────────────
    print("\nDAG:")
    for lid, node in all_nodes.items():
        parents = [p[:8] for p in node["parents"]]
        ops = [o["method"].split(".")[-1] for o in node["operations"]]
        ops_str = f"  ops=[{', '.join(ops)}]" if ops else ""
        print(f"  {lid[:8]}  parents={parents}  desc={node['description']}{ops_str}")

    print("\nAll checks passed.")


if __name__ == "__main__":
    test_lineagetracker()
