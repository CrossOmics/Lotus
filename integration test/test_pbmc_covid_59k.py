"""
Integration test: lotus I/O + cplearn core_selection on covid PBMC 59k dataset.

Uses precomputed covid_cplearn.h5ad (load + core_selection already run), so tests
reuse layer and only re-run percentage selection (fast). Generate that file once with:

    python "integration test/make_covid_cplearn.py"

(Requires covid.h5ad in repo root; writes covid_cplearn.h5ad there.)

Run tests:
    cd "integration test" && python test_pbmc_covid_59k.py
    # or: pytest "integration test/test_pbmc_covid_59k.py" -v -s
"""
import sys
import json
import tempfile
import numpy as np
from pathlib import Path

import pytest

project_root = Path(__file__).resolve().parent.parent
src_root = project_root / "src"
for p in (str(project_root), str(src_root)):
    if p not in sys.path:
        sys.path.insert(0, p)

# Precomputed cplearn result (generate with make_covid_cplearn.py)
DATA_PATH = Path("/Users/liuzhenke/Desktop/Lotus/covid_cplearn.h5ad")
PERCENTAGE = 0.30
CPLEARN_COLS = ("cplearn", "cplearn_cluster", "cplearn_is_core", "cplearn_layer", "cplearn_core_selection")


def _ensure_pca(adata):
    """Run a minimal preprocessing pass if X_pca is missing."""
    if "X_pca" in adata.obsm:
        return adata

    print("  [prep] X_pca not found – running lightweight preprocessing ...")
    import lotus.preprocessing as pp

    pp.normalize_total(adata, target_sum=1e4)
    pp.log1p(adata)
    pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
    adata = adata[:, adata.var["highly_variable"]].copy()
    pp.scale(adata, max_value=10)
    pp.pca(adata, n_comps=50)
    print(f"  [prep] X_pca computed: {adata.obsm['X_pca'].shape}")
    return adata


def test_io_and_core_selection():
    """End-to-end: load → (optional prep) → core_selection → assertions (uses covid.h5ad)."""
    import lotus.io as io
    from lotus.core_analysis import core_selection

    if not DATA_PATH.exists():
        pytest.skip(f"Dataset not found: {DATA_PATH} (optional for CI)")

    # ── 1. I/O ──────────────────────────────────────────────────────────────
    print(f"\n[1] Loading {DATA_PATH.name} ...")

    adata = io.standardize_load(DATA_PATH)
    print(f"    Loaded: {adata.n_obs} cells × {adata.n_vars} genes")
    assert adata.n_obs > 0, "AnnData has no cells after loading."
    assert adata.n_vars > 0, "AnnData has no genes after loading."

    # ── 2. Ensure PCA ────────────────────────────────────────────────────────
    adata = _ensure_pca(adata)
    assert "X_pca" in adata.obsm, "X_pca must be present before core_selection."

    # ── 3. core_selection ───────────────────────────────────────────────────
    print(f"\n[3] Running core_selection (percentage={PERCENTAGE}) ...")
    adata, model = core_selection(adata, PERCENTAGE)

    # ── 4. Assertions ────────────────────────────────────────────────────────
    print("\n[4] Running assertions ...")
    n_obs = adata.n_obs

    # 4a. Required columns exist
    for col in CPLEARN_COLS:
        assert col in adata.obs, f"Missing obs column: '{col}'"
    assert "cplearn" in adata.uns, "Missing adata.uns['cplearn']"

    # 4b. core_selection only contains expected values
    vals = set(adata.obs["cplearn_core_selection"].unique())
    assert vals <= {"core", "noncore"}, f"Unexpected categories: {vals}"

    # 4c. Core uses whole layers only: n_core >= target (may overshoot), never split a layer
    n_core = (adata.obs["cplearn_core_selection"] == "core").sum()
    n_noncore = (adata.obs["cplearn_core_selection"] == "noncore").sum()
    target_min = int(np.ceil(PERCENTAGE * n_obs))
    assert n_core >= target_min, (
        f"Core cells {n_core} should be >= target {target_min} (whole-layer policy)."
    )
    assert n_core + n_noncore == n_obs, "core + noncore must equal total cells."

    # 4d. Actual fraction stored in uns; layers_in_core lists whole layers included
    cs_meta = adata.uns["cplearn"]["core_selection"]
    assert cs_meta["n_core"] == int(n_core)
    assert cs_meta["n_noncore"] == int(n_noncore)
    assert abs(cs_meta["actual_fraction_core"] - n_core / n_obs) < 1e-9
    assert "layers_in_core" in cs_meta
    layer_vals = adata.obs["cplearn_layer"].values
    is_core_mask = adata.obs["cplearn_core_selection"].values == "core"
    # Every core cell must be in one of the declared whole layers
    assert np.all(np.isin(layer_vals[is_core_mask], cs_meta["layers_in_core"]))

    # 4f. uns layers_indices_json is valid JSON and non-empty
    layers_json = adata.uns["cplearn"]["layers_indices_json"]
    layers = json.loads(layers_json)
    assert isinstance(layers, list) and len(layers) > 0, (
        "layers_indices_json should be a non-empty list."
    )

    # ── Summary ──────────────────────────────────────────────────────────────
    frac_str = f"{n_core / n_obs:.1%}"
    print(f"\n✅ All checks passed!")
    print(f"   Dataset   : {n_obs:,} cells × {adata.n_vars:,} genes")
    print(f"   Requested : {PERCENTAGE:.0%} core  →  actual {frac_str} ({n_core:,} core / {n_noncore:,} noncore)")
    layer_dist = {
        int(k): int(v)
        for k, v in zip(*np.unique(layer_vals, return_counts=True))
    }
    print(f"   Layer distribution (index → count): {layer_dist}")


def test_io_roundtrip_after_core_selection():
    """Write adata (with cplearn results) to h5ad and read back; assert key columns/uns preserved (uses covid.h5ad)."""
    from anndata import AnnData
    import lotus.io as io
    import lotus.preprocessing as pp
    from lotus.core_analysis import core_selection

    if not DATA_PATH.exists():
        pytest.skip(f"Dataset not found: {DATA_PATH}")

    adata = io.standardize_load(DATA_PATH)
    adata = _ensure_pca(adata)
    adata, _ = core_selection(adata, 0.20)

    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as f:
        path = Path(f.name)
    try:
        io.write(path, adata)
        adata2 = io.read_h5ad(path)
        for col in CPLEARN_COLS:
            assert col in adata2.obs, f"After roundtrip missing: {col}"
        assert "cplearn" in adata2.uns
        assert "core_selection" in adata2.uns["cplearn"]
        assert "layers_indices_json" in adata2.uns["cplearn"]
        n_core = (adata2.obs["cplearn_core_selection"] == "core").sum()
        assert n_core == adata.uns["cplearn"]["core_selection"]["n_core"]
    finally:
        path.unlink(missing_ok=True)


def test_core_percentages_10_30_50():
    """Load covid_cplearn.h5ad, run core_selection at 10%, 30%, 50%; check n_core and layers."""
    import lotus.io as io
    from lotus.core_analysis import core_selection

    if not DATA_PATH.exists():
        pytest.skip(f"Dataset not found: {DATA_PATH}")

    adata = io.standardize_load(DATA_PATH)
    n_obs = adata.n_obs

    for pct in (0.10, 0.30, 0.50):
        adata_run, _ = core_selection(adata, pct)
        n_core = (adata_run.obs["cplearn_core_selection"] == "core").sum()
        target_min = int(np.ceil(pct * n_obs))
        assert n_core >= target_min, f"percentage={pct}: need at least {target_min} core (whole layers), got {n_core}"


if __name__ == "__main__":
    pytest.main([str(Path(__file__).resolve()), "-v", "-s"])

