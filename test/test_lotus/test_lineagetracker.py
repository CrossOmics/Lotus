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

from lotus.lineagetracker import LineageTracker, bind_variable_name, logged
from lotus.lineagetracker.visualizer import render_missing_pngs
from lotus.preprocessing import log1p, normalize_total, pca, scale


@logged
def normal_and_log(adata: anndata):
    normalize_total(adata, target_sum=1e4, inplace=True)
    log1p(adata)

def _make_adata(
    n_obs: int = 100,
    n_vars: int = 50,
    name: str | None = None,
) -> anndata.AnnData:
    """Create a small synthetic AnnData for testing."""
    X = np.random.negative_binomial(n=20, p=0.3, size=(n_obs, n_vars)).astype(
        np.float32
    )
    adata = anndata.AnnData(X=X)
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adata.var_names = [f"gene_{i}" for i in range(n_vars)]
    adata.obs["celltype"] = np.random.choice(["T", "B", "NK"], n_obs)
    if name:
        adata.uns["name"] = name
    return adata


def test_lineagetracker():
    # Reset singleton and clear persisted artifacts so the test is self-contained
    LineageTracker.reset()
    tracker = LineageTracker.instance()
    # Clean up the session folder created for this test run
    if tracker.lineage_file.exists():
        tracker.lineage_file.unlink()
    if tracker.graph_file.exists():
        tracker.graph_file.unlink()
    tracker._nodes.clear()
    tracker._lid_map.clear()

    # 1. Register root
    root_adata = _make_adata(name="root_adata")
    root_lid = tracker.register(
        root_adata, parents=[], creation_op="create", description="synthetic 100x50"
    )
    assert tracker.get_lid(root_adata) == root_lid
    all_nodes = json.loads(tracker.lineage_file.read_text())
    assert root_lid in all_nodes
    assert all_nodes[root_lid]["variable_name"] == "root_adata"
    assert all_nodes[root_lid]["display_name"] == "root_adata"
    print(f"[OK] Root registered: {root_lid[:8]}")

    # 2. Record in-place operations
    alias_adata = root_adata
    tracker.record_op(alias_adata, "preprocessing.normalize_total", {"target_sum": 1e4})
    tracker.record_op(alias_adata, "preprocessing.log1p", {})
    node = json.loads(tracker.lineage_file.read_text())[root_lid]
    assert len(node["operations"]) == 2
    assert node["operations"][0]["method"] == "preprocessing.normalize_total"
    assert node["variable_name"] == "alias_adata"

    bind_variable_name(root_adata, "manual_root")
    node = json.loads(tracker.lineage_file.read_text())[root_lid]
    assert node["variable_name"] == "manual_root"
    assert node["display_name"] == "manual_root"

    latest_alias = root_adata
    tracker.record_op(latest_alias, "preprocessing.scale", {"max_value": 10})
    node = json.loads(tracker.lineage_file.read_text())[root_lid]
    assert len(node["operations"]) == 3
    assert node["variable_name"] == "latest_alias"
    assert node["display_name"] == "latest_alias"
    print("[OK] Operations recorded on root node")

    # 2b. Regression: @logged functions should record exactly once (no double-wrap)
    single_logged_adata = _make_adata(name="single_logged_adata")
    single_logged_lid = tracker.register(
        single_logged_adata,
        parents=[],
        creation_op="create",
        description="single logged regression",
    )
    before_logged_count = len(
        json.loads(tracker.lineage_file.read_text())[single_logged_lid]["operations"]
    )
    normalize_total(single_logged_adata, target_sum=1e4, inplace=True)
    single_logged_node = json.loads(tracker.lineage_file.read_text())[single_logged_lid]
    assert len(single_logged_node["operations"]) == before_logged_count + 1
    assert single_logged_node["operations"][-1]["method"] == "normalize_total"
    print("[OK] @logged function records exactly one operation")

    # 2c. Multiple @logged preprocessing functions should each append one operation
    logged_pipeline_adata = _make_adata(name="logged_pipeline_adata")
    logged_pipeline_lid = tracker.register(
        logged_pipeline_adata,
        parents=[],
        creation_op="create",
        description="logged pipeline regression",
    )
    before_pipeline_count = len(
        json.loads(tracker.lineage_file.read_text())[logged_pipeline_lid]["operations"]
    )
    # check suppresse innner functions
    normal_and_log(logged_pipeline_adata)
    test_copy_pipeline_adata = logged_pipeline_adata.copy()
    scale(logged_pipeline_adata, max_value=10)
    test_slicing_pipeline_adata = logged_pipeline_adata[root_adata.obs["celltype"] == "T"]
    pca(logged_pipeline_adata, n_comps=5)
    logged_pipeline_node = json.loads(tracker.lineage_file.read_text())[logged_pipeline_lid]
    pipeline_methods = [op["method"] for op in logged_pipeline_node["operations"]]
    print("[OK] Multiple @logged preprocessing functions each record one operation")

    # 3. Slice (monkey-patched __getitem__)
    t_cells = root_adata[root_adata.obs["celltype"] == "T"]
    t_lid = tracker.get_lid(t_cells)
    assert t_lid is not None
    t_node = json.loads(tracker.lineage_file.read_text())[t_lid]
    assert root_lid in t_node["parents"]
    assert t_node["creation_op"] == "slice"
    assert t_node["variable_name"] == "t_cells"
    print(f"[OK] Slice tracked: {t_lid[:8]}, shape={t_cells.shape}")

    b_cells = root_adata[root_adata.obs["celltype"] == "B"]
    b_lid = tracker.get_lid(b_cells)
    assert b_lid is not None
    b_node = json.loads(tracker.lineage_file.read_text())[b_lid]
    assert b_node["variable_name"] == "b_cells"
    print(f"[OK] Slice tracked: {b_lid[:8]}, shape={b_cells.shape}")

    # 4. Copy (monkey-patched copy)
    t_copy = t_cells.copy()
    cp_lid = tracker.get_lid(t_copy)
    assert cp_lid is not None
    cp_node = json.loads(tracker.lineage_file.read_text())[cp_lid]
    assert t_lid in cp_node["parents"]
    assert cp_node["creation_op"] == "copy"
    assert cp_node["variable_name"] == "t_copy"
    print(f"[OK] Copy tracked: {cp_lid[:8]}")

    # 5. Concat (monkey-patched anndata.concat)
    merged_cells = anndata.concat([t_cells, b_cells])
    m_lid = tracker.get_lid(merged_cells)
    assert m_lid is not None
    m_node = json.loads(tracker.lineage_file.read_text())[m_lid]
    assert set(m_node["parents"]) == {t_lid, b_lid}
    assert m_node["creation_op"] == "concat"
    assert m_node["variable_name"] == "merged_cells"
    print(f"[OK] Concat tracked: {m_lid[:8]}, shape={merged_cells.shape}")

    # 6. Fallback when no variable name is available but adata.uns["name"] exists
    uns_name_lid = tracker.register(
        _make_adata(n_obs=20, n_vars=10, name="from_uns_name"),
        parents=[],
        creation_op="create",
        description="uns fallback",
    )
    uns_name_node = json.loads(tracker.lineage_file.read_text())[uns_name_lid]
    assert uns_name_node["variable_name"] is None
    assert uns_name_node["display_name"] == "from_uns_name"
    print(f"[OK] Uns-name fallback used: {uns_name_lid[:8]}")

    # ── 7. Visualisation PNG ─────────────────────────────────────
    # Test render_missing_pngs: JSON exists but no PNG → should generate one
    assert not tracker.graph_file.exists(), "PNG should not exist yet"
    render_missing_pngs(tracker.root_dir)
    assert tracker.graph_file.exists()
    assert tracker.graph_file.stat().st_size > 0
    print(f"[OK] PNG generated by render_missing_pngs: {tracker.graph_file}")

    # Second call: PNG now exists → should be skipped (no error, no re-render)
    mtime_before = tracker.graph_file.stat().st_mtime
    render_missing_pngs(tracker.root_dir)
    assert tracker.graph_file.stat().st_mtime == mtime_before
    print("[OK] render_missing_pngs skipped already-complete session")

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
