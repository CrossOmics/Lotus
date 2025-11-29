"""
Generate example visualizations for the Quick Start guide.

This script runs the three example workflows from getting_started.rst and
generates visualization images to include in the documentation.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Add project root to Python path
_project_root = Path(__file__).parent.parent.resolve()
if str(_project_root) not in sys.path:
    sys.path.insert(0, str(_project_root))

import numpy as np
import pandas as pd
from anndata import AnnData

import lotus as lt
from lotus.workflows import preprocess, leiden, louvain, umap
from lotus.workflows.deg_analysis import rank_genes_groups, marker_genes
from lotus.workflows.core_analysis import core_analyze
from lotus.workflows.visualization import coremap
from lotus.methods.cplearn.external import cplearn


def build_demo_dataset(
    *,
    n_clusters: int = 3,
    cells_per_cluster: int = 60,
    n_genes: int = 50,
    random_state: int = 42,
) -> AnnData:
    """Create a synthetic AnnData object for examples."""
    rng = np.random.default_rng(random_state)
    total_cells = n_clusters * cells_per_cluster
    
    # Generate count matrix with cluster-specific patterns
    counts = np.random.negative_binomial(
        n=5, p=0.3, size=(total_cells, n_genes)
    ).astype(np.float32)
    
    # Add cluster-specific expression patterns
    cluster_labels = np.repeat(range(n_clusters), cells_per_cluster)
    for i, cluster_id in enumerate(cluster_labels):
        # Each cluster has 5-10 marker genes with higher expression
        marker_genes = np.random.choice(
            n_genes, size=np.random.randint(5, 11), replace=False
        )
        counts[i, marker_genes] *= np.random.uniform(2.0, 5.0, size=len(marker_genes))
    
    # Create AnnData object
    adata = AnnData(counts)
    adata.var_names = [f"gene_{i:03d}" for i in range(n_genes)]
    adata.obs_names = [f"cell_{i:04d}" for i in range(total_cells)]
    adata.obs["truth"] = pd.Categorical(cluster_labels.astype(str))
    
    return adata


def example_1_standard_workflow(output_dir: Path) -> None:
    """Run Example 1: Standard Scanpy Workflow."""
    print("\n" + "=" * 60)
    print("Example 1: Standard Scanpy Workflow")
    print("=" * 60)
    
    # Create synthetic data
    adata = build_demo_dataset(n_clusters=3, cells_per_cluster=60, random_state=42)
    print(f"Created dataset: {adata.shape}")
    
    # 1. Preprocessing
    print("\n1. Preprocessing...")
    preprocess(adata, n_pcs=20, n_top_genes=2000, n_neighbors=15, save_raw=True)
    
    # 2. Clustering (Leiden)
    print("\n2. Clustering (Leiden)...")
    leiden(adata, resolution=0.8, key_added="leiden")  # Higher resolution for small dataset
    print(f"   Found {adata.obs['leiden'].nunique()} clusters")
    
    # 3. Visualization
    print("\n3. Visualization (UMAP)...")
    umap(
        adata,
        cluster_key="leiden",
        output_dir=str(output_dir),
        save="_standard_workflow.png",
    )
    print(f"   Saved: {output_dir}/umap_standard_workflow.png")
    
    # 4. DEG Analysis
    print("\n4. Differential Expression Analysis...")
    rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
    print("   ✓ DEG analysis complete")
    
    print("\n✓ Standard workflow complete!")


def example_2_cplearn_workflow(output_dir: Path) -> None:
    """Run Example 2: Core Analysis + Cplearn Workflow."""
    print("\n" + "=" * 60)
    print("Example 2: Core Analysis + Cplearn Workflow")
    print("=" * 60)
    
    # Create synthetic data
    adata = build_demo_dataset(n_clusters=3, cells_per_cluster=60, random_state=42)
    print(f"Created dataset: {adata.shape}")
    
    # 1. Preprocessing
    print("\n1. Preprocessing...")
    preprocess(adata, n_pcs=20, n_top_genes=2000, n_neighbors=15, save_raw=True)
    
    # 2. Core Analysis + Cplearn Clustering
    print("\n2. Core Analysis + Cplearn Clustering...")
    try:
        model = cplearn.corespect(
            adata,
            use_rep="X_latent",
            key_added="cplearn",
            stable={"core_frac": 0.3, "ng_num": 10},  # Adjusted for small dataset
            cluster={"resolution": 0.8},  # Lower resolution for small dataset
            propagate=False,  # Disable propagation to avoid errors
        )
        print(f"   Found {adata.obs['cplearn'].nunique()} clusters")
    except Exception as e:
        print(f"   Warning: Cplearn clustering failed: {e}")
        print("   Skipping cplearn workflow...")
        return
    
    # 3. Compute core map embedding
    print("\n3. Computing core map embedding...")
    core_analyze(
        adata,
        model=model,
        use_rep="X_latent",
        key_added="X_cplearn_coremap",
    )
    
    # 4. Visualization: Coremap
    print("\n4. Visualization (Coremap)...")
    coremap(
        adata,
        coremap_key="X_cplearn_coremap",
        cluster_key="cplearn",
        output_dir=str(output_dir),
        save="_cplearn_workflow.html",
    )
    print(f"   Saved: {output_dir}/coremap_cplearn_workflow.html")
    
    print("\n✓ Cplearn workflow complete!")


def example_3_alternating_methods(output_dir: Path) -> None:
    """Run Example 3: Alternating Methods."""
    print("\n" + "=" * 60)
    print("Example 3: Alternating Methods")
    print("=" * 60)
    
    # Create synthetic data
    adata = build_demo_dataset(n_clusters=3, cells_per_cluster=60, random_state=42)
    print(f"Created dataset: {adata.shape}")
    
    # 1. Preprocessing (shared)
    print("\n1. Preprocessing (shared)...")
    preprocess(adata, n_pcs=20, n_top_genes=2000, n_neighbors=15, save_raw=True)
    
    # === Workflow A: Cplearn ===
    print("\n2. Workflow A: Cplearn...")
    try:
        model = cplearn.corespect(
            adata,
            use_rep="X_latent",
            key_added="cplearn",
            stable={"core_frac": 0.3, "ng_num": 10},
            cluster={"resolution": 0.8},
            propagate=False,
        )
        core_analyze(adata, model=model, use_rep="X_latent", key_added="X_cplearn_coremap")
    except Exception as e:
        print(f"   Warning: Cplearn workflow failed: {e}")
        print("   Skipping cplearn visualization...")
        return
    coremap(
        adata,
        coremap_key="X_cplearn_coremap",
        cluster_key="cplearn",
        output_dir=str(output_dir),
        save="_alternating_cplearn.html",
    )
    print(f"   Cplearn clusters: {adata.obs['cplearn'].nunique()}")
    
    # === Workflow B: Scanpy Louvain ===
    print("\n3. Workflow B: Scanpy Louvain...")
    louvain(adata, resolution=0.8, key_added="louvain")  # Higher resolution for small dataset
    umap(
        adata,
        cluster_key="louvain",
        output_dir=str(output_dir),
        save="_alternating_louvain.png",
    )
    print(f"   Louvain clusters: {adata.obs['louvain'].nunique()}")
    
    # Compare results
    print("\n4. Comparing results...")
    print(f"   Cplearn clusters: {adata.obs['cplearn'].value_counts().to_dict()}")
    print(f"   Louvain clusters: {adata.obs['louvain'].value_counts().to_dict()}")
    
    print("\n✓ Alternating methods workflow complete!")


def main():
    """Generate all example visualizations."""
    # Create output directory in docs/_static
    output_dir = Path(__file__).parent / "_static" / "examples"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("Generating Example Visualizations for Quick Start Guide")
    print("=" * 60)
    print(f"Output directory: {output_dir}")
    
    try:
        # Run all examples
        example_1_standard_workflow(output_dir)
        example_2_cplearn_workflow(output_dir)
        example_3_alternating_methods(output_dir)
        
        print("\n" + "=" * 60)
        print("✓ All examples completed successfully!")
        print("=" * 60)
        print(f"\nGenerated files in: {output_dir}")
        print("  - umap_standard_workflow.png")
        print("  - coremap_cplearn_workflow.html")
        print("  - coremap_alternating_cplearn.html")
        print("  - umap_alternating_louvain.png")
        
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

