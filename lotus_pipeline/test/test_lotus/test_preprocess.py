import sys
import unittest
import numpy as np
import pandas as pd
from pathlib import Path
from anndata import AnnData

import lotus_pipeline.lotus.preprocessing as pp

project_root = Path(__file__).resolve().parent.parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))


def test_run_preprocessing():
    """
    Independent unit test function: Test the complete workflow of run_preprocessing.
    """
    print(f"\n[Test] Project root set to: {project_root}")
    print("[Test] Generating synthetic AnnData...")

    # --- A. Prepare Synthetic Data ---
    n_obs, n_vars = 100, 500
    # Generate random count data using negative binomial distribution
    counts = np.random.negative_binomial(n=20, p=0.2, size=(n_obs, n_vars))

    # Construct AnnData object
    adata = AnnData(X=counts.astype(np.float32))
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]

    # Simulate mitochondrial genes (MT-) for QC testing
    mt_genes = [f"MT-gene_{i}" for i in range(10)]
    other_genes = [f"gene_{i}" for i in range(10, n_vars)]
    adata.var_names = mt_genes + other_genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')

    # --- B. Execute Core Function ---
    print("[Test] Running preprocess pipeline...")
    pp.run_preprocessing(
        adata,
        n_pcs=10,
        n_top_genes=100,
        target_sum=1e4,
        n_neighbors=5,
        save_raw=True,
        raw_layer="raw_counts",
        pct_mt_max=20.0,
        hvg_flavor="seurat",
        qc_vars=['mt']
    )

    # --- C. Verify Results (Assertions) ---
    print("[Test] Verifying results...")

    # 1. Verify raw data backup
    assert "raw_counts" in adata.layers, "❌ Failed: Raw counts layer missing."

    # 2. Verify QC metrics generation
    assert "n_genes_by_counts" in adata.obs, "❌ Failed: QC metric n_genes_by_counts missing."
    assert "pct_counts_mt" in adata.obs, "❌ Failed: Mitochondrial percentage missing."

    # 3. Verify normalization status
    # Standard Scanpy workflow: .raw stores normalized + log-transformed data
    assert adata.raw is not None, "❌ Failed: adata.raw should be populated."

    # 4. Verify Highly Variable Genes (HVG)
    n_hvg = sum(adata.var["highly_variable"])
    assert n_hvg == 100, f"❌ Failed: Expected 100 HVGs, got {n_hvg}."

    # 5. Verify Dimensionality Reduction (PCA)
    assert "X_pca" in adata.obsm, "❌ Failed: X_pca missing."
    assert adata.obsm["X_pca"].shape == (adata.n_obs, 10), "❌ Failed: PCA shape mismatch."

    # 6. Verify Neighbor Graph
    assert "neighbors" in adata.uns, "❌ Failed: Neighbors config missing."
    assert "connectivities" in adata.obsp, "❌ Failed: Connectivities matrix missing."

    print("\n✅ All checks passed! run_preprocessing is working correctly.")
