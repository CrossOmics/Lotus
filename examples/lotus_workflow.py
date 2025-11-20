from __future__ import annotations

import argparse
import logging
import sys
from datetime import datetime
from pathlib import Path

# Add project root to Python path automatically
# This allows importing lotus without setting PYTHONPATH manually
_project_root = Path(__file__).parent.parent.resolve()
if str(_project_root) not in sys.path:
    sys.path.insert(0, str(_project_root))

import numpy as np
import pandas as pd
from anndata import AnnData

import lotus as lt
from lotus.workflows import (
    # Preprocess
    preprocess,
    # Visualization
    umap as workflow_umap,
    render_visualizations,
    # Clustering
    clustering,
    # DEG
    marker_genes,
    # CoreSelection
    core_selection,
)


def build_demo_dataset(
    *,
    n_clusters: int = 3,
    cells_per_cluster: int = 60,
    n_genes: int = 50,
    n_latent: int = 12,
    random_state: int = 0,
) -> AnnData:
    """
    Create a synthetic AnnData object mimicking count data prior to preprocessing.
    """
    # Create a controllable Gaussian rate distribution dataset for local workflow reproduction
    rng = np.random.default_rng(random_state)
    total_cells = n_clusters * cells_per_cluster

    latent_repr = np.empty((total_cells, n_latent), dtype=np.float32)
    counts = np.empty((total_cells, n_genes), dtype=np.float32)
    cluster_labels = np.empty(total_cells, dtype=int)

    gene_names = [f"gene_{i:03d}" for i in range(n_genes)]
    obs_names = [f"cell_{i:03d}" for i in range(total_cells)]

    for cluster in range(n_clusters):
        # Generate cluster-specific latent centers and count rates for each cluster to create differential signals
        start = cluster * cells_per_cluster
        end = start + cells_per_cluster

        latent_center = rng.normal(loc=cluster * 2.5, scale=0.4, size=n_latent)
        latent_repr[start:end] = rng.normal(
            loc=latent_center,
            scale=0.35,
            size=(cells_per_cluster, n_latent),
        )

        base_profile = rng.gamma(shape=2.0 + cluster, scale=1.2, size=n_genes)
        noise = rng.normal(loc=0.0, scale=0.6, size=(cells_per_cluster, n_genes))
        rate = np.clip(base_profile + noise, a_min=0.1, a_max=None)
        counts[start:end] = rng.poisson(rate).astype(np.float32)

        cluster_labels[start:end] = cluster

    obs = pd.DataFrame(index=obs_names)
    obs["truth"] = pd.Categorical(cluster_labels)
    var = pd.DataFrame(index=gene_names)

    adata = AnnData(counts, obs=obs, var=var)
    adata.obsm["X_latent_init"] = latent_repr
    return adata


def setup_logging(output_dir: Path) -> logging.Logger:
    """
    Set up logging to both console and file.
    
    Parameters:
        output_dir: Directory to save log file
    
    Returns:
        Logger instance
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create log file with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = output_dir / f"workflow_{timestamp}.log"
    
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode="w", encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    
    logger = logging.getLogger(__name__)
    logger.info(f"Logging initialized. Log file: {log_file}")
    return logger


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Demonstrate lotus workflows on a synthetic dataset.",
    )
    parser.add_argument(
        "--clusters",
        type=int,
        default=3,
        help="Number of ground-truth clusters to simulate (default: 3).",
    )
    parser.add_argument(
        "--cells-per-cluster",
        type=int,
        default=60,
        help="Number of cells per cluster (default: 60).",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Directory to store results. If not specified, creates 'result_YYYYMMDD_HHMMSS' in current directory.",
    )
    args = parser.parse_args()
    
    # Create output directory with timestamp if not specified
    if args.output_dir is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = Path(f"result_{timestamp}")
    else:
        output_dir = Path(args.output_dir)
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Set up logging
    logger = setup_logging(output_dir)

    logger.info("=" * 60)
    logger.info("Lotus Workflow - Starting Analysis")
    logger.info("=" * 60)
    
    adata = build_demo_dataset(
        n_clusters=args.clusters,
        cells_per_cluster=args.cells_per_cluster,
    )
    # Display basic dataset shape for quick validation of simulation scale
    logger.info(f"Generated AnnData with shape {adata.shape} and {adata.n_vars} genes.")

    # ============================================
    # Preprocess: QC → Filtering → Normalization → HVG → Scaling → PCA → Neighbors
    # ============================================
    logger.info("\n" + "=" * 60)
    logger.info("Preprocessing Pipeline")
    logger.info("=" * 60)
    logger.info("Running complete preprocessing pipeline...")
    preprocess(
        adata,
        n_pcs=20,
        target_sum=1e4,
        n_top_genes=min(2000, adata.n_vars),
        n_neighbors=15,
        use_rep="X_pca",
        save_raw=True,
        raw_layer="raw_counts",
    )
    logger.info("✓ Preprocessing complete")
    logger.info(f"  - Data shape: {adata.shape}")
    logger.info(f"  - PCA stored in: `adata.obsm['X_pca']` (shape: {adata.obsm['X_pca'].shape})")
    logger.info(f"  - Latent representation stored in: `adata.obsm['X_latent']` (shape: {adata.obsm['X_latent'].shape})")
    logger.info(f"  - Raw counts saved in: `adata.layers['raw_counts']`")
    logger.info(f"  - Neighbors graph constructed: {adata.obsp.get('distances') is not None}")

    # ============================================
    # Clustering: clustering
    # ============================================
    logger.info("\n" + "=" * 60)
    logger.info("Clustering")
    logger.info("=" * 60)
    logger.info("Performing clustering analysis...")
    model = clustering(
        adata,
        use_rep="X_latent",
        key_added="cplearn_labels",
        stable_core_frac=0.25,
        stable_ng_num=8,
        cluster_resolution=1.2,
        fine_grained=False,
        propagate=True,
        print_summary=True,
    )
    logger.info("✓ Clustering complete")
    logger.info(f"  - Cluster labels stored in: `adata.obs['cplearn_labels']`")
    unique_labels = adata.obs["cplearn_labels"].unique()
    logger.info(f"  - Number of clusters: {len(unique_labels)}")
    logger.info(f"  - Cluster IDs: {sorted(unique_labels.tolist())}")

    # ============================================
    # Visualization: UMAP
    # ============================================
    logger.info("\n" + "=" * 60)
    logger.info("Visualization: UMAP")
    logger.info("=" * 60)
    logger.info("Computing UMAP embedding and generating visualization...")
    workflow_umap(
        adata,
        cluster_key="cplearn_labels",
        truth_key="truth",
        output_dir=output_dir,
        save="_clusters.png",
        show=False,
        compute_umap=True,
    )
    logger.info("✓ UMAP visualization complete")
    logger.info(f"  - UMAP embedding stored in: `adata.obsm['X_umap']` (shape: {adata.obsm['X_umap'].shape})")
    logger.info(f"  - Visualization saved to: {output_dir / 'umap_clusters.png'}")

    # ============================================
    # CoreSelection: Neighbors → CoreSelection
    # ============================================
    logger.info("\n" + "=" * 60)
    logger.info("CoreSelection: Neighbors → CoreSelection")
    logger.info("=" * 60)
    logger.info("Computing core map embedding...")
    core_selection(
        adata,
        model=model,
        use_rep="X_latent",
        key_added="X_cplearn_coremap",
        print_summary=True,
    )
    logger.info("✓ CoreSelection complete")
    embedding = adata.obsm["X_cplearn_coremap"]
    assigned = np.sum(~np.isnan(embedding).any(axis=1))
    logger.info(f"  - Core map embedding stored in: `adata.obsm['X_cplearn_coremap']` (shape: {embedding.shape})")
    logger.info(f"  - Assigned points: {assigned}/{adata.n_obs} ({100*assigned/adata.n_obs:.1f}%)")

    # ============================================
    # DEG: Marker Genes
    # ============================================
    logger.info("\n" + "=" * 60)
    logger.info("DEG: Marker Genes")
    logger.info("=" * 60)
    logger.info("Identifying differentially expressed genes (marker genes)...")
    de_result = marker_genes(
        adata,
        cluster_key="cplearn_labels",
        layer="raw_counts",
        min_detect_pct=0.0,
        min_cells_per_group=5,
        auto_pick_groups=True,
    )
    logger.info("✓ DEG analysis complete")
    logger.info(f"  - Total differentially expressed genes: {len(de_result)}")
    logger.info(f"  - Significant genes (p_adj < 0.05): {np.sum(de_result['p_adj'] < 0.05)}")
    logger.info(f"  - Significant genes (p_adj < 0.01): {np.sum(de_result['p_adj'] < 0.01)}")
    
    logger.info("\nTop 10 differentially expressed genes:")
    cols = [
        "gene",
        "log2fc",
        "z_score",
        "pvalue",
        "p_adj",
        "mean_a",
        "mean_b",
        "pct_expr_a",
        "pct_expr_b",
    ]
    printable = de_result[cols].head(10)
    logger.info("\n" + printable.to_string(index=False))
    
    # Save DEG results to CSV
    deg_output_file = output_dir / "deg_results.csv"
    de_result.to_csv(deg_output_file, index=False)
    logger.info(f"  - DEG results saved to: {deg_output_file}")

    top_markers = (
        de_result["gene"].head(5).astype(str).tolist()
        if not de_result.empty
        else []
    )

    # ============================================
    # Visualization: Render all visualizations
    # ============================================
    logger.info("\n" + "=" * 60)
    logger.info("Visualization: Rendering all visualizations")
    logger.info("=" * 60)
    render_visualizations(adata, top_markers, output_dir)
    
    # Move any images from figures/ to output_dir if they exist
    # (scanpy may save to default figures/ directory despite figdir setting)
    import shutil
    figures_dir = Path("figures")
    if figures_dir.exists():
        for img_file in figures_dir.glob("*.png"):
            if img_file.name.startswith("umap") or img_file.name.startswith("violin"):
                dest_file = output_dir / img_file.name
                if img_file.exists() and not dest_file.exists():
                    shutil.move(str(img_file), str(dest_file))
                    logger.info(f"  - Moved {img_file.name} to output directory")
    
    logger.info("✓ All visualizations complete")
    logger.info(f"  - UMAP clusters plot: {output_dir / 'umap_clusters.png'}")
    if top_markers:
        logger.info(f"  - Marker genes violin plot: {output_dir / 'violin_markers.png'}")
        logger.info(f"  - Top {len(top_markers)} marker genes visualized: {', '.join(top_markers)}")
    
    # Summary
    logger.info("\n" + "=" * 60)
    logger.info("Workflow Summary")
    logger.info("=" * 60)
    logger.info(f"✓ All results saved to: {output_dir.resolve()}")
    logger.info(f"✓ Output directory: {output_dir.name}")
    logger.info(f"  - Log file: {list(output_dir.glob('workflow_*.log'))[0] if list(output_dir.glob('workflow_*.log')) else 'N/A'}")
    logger.info(f"  - DEG results: {output_dir / 'deg_results.csv'}")
    logger.info(f"  - UMAP plot: {output_dir / 'umap_clusters.png'}")
    if top_markers:
        logger.info(f"  - Violin plot: {output_dir / 'violin_markers.png'}")
    logger.info("=" * 60)
    logger.info("Workflow completed successfully!")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
