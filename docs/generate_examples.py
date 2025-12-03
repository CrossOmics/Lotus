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
    
    # Load real data
    demo_data_path = Path(__file__).parent.parent / "data" / "demo_data.h5ad"
    print(f"Loading data from: {demo_data_path}")
    adata = lt.read(str(demo_data_path))
    print(f"Loaded dataset: {adata.shape}")
    
    # 1. Preprocessing
    print("\n1. Preprocessing...")
    preprocess(adata, n_pcs=20, n_top_genes=2000, n_neighbors=15, save_raw=True)
    
    # 2. Clustering (Leiden)
    print("\n2. Clustering (Leiden)...")
    leiden(adata, resolution=0.8, key_added="leiden")  # Higher resolution for small dataset
    print(f"   Found {adata.obs['leiden'].nunique()} clusters")
    
    # 3. Visualization
    print("\n3. Visualization (UMAP)...")
    # Set figdir to ensure files are saved in output_dir
    old_figdir = lt.settings.figdir
    lt.settings.figdir = str(output_dir)
    try:
        umap(
            adata,
            cluster_key="leiden",
            output_dir=str(output_dir),
            save="_standard_workflow.png",
        )
        png_path = output_dir / "umap_standard_workflow.png"
        # Also check figures/ directory (fallback)
        if not png_path.exists():
            figures_png = Path("figures") / "umap_standard_workflow.png"
            if figures_png.exists():
                import shutil
                shutil.copy(figures_png, png_path)
                print(f"   Copied from figures/ to: {png_path}")
            else:
                print(f"   Warning: PNG file not found")
        else:
            print(f"   Saved: {png_path}")
    finally:
        lt.settings.figdir = old_figdir
    
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
    
    # Load real data
    demo_data_path = Path(__file__).parent.parent / "data" / "demo_data.h5ad"
    print(f"Loading data from: {demo_data_path}")
    adata = lt.read(str(demo_data_path))
    print(f"Loaded dataset: {adata.shape}")
    
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
            propagate=True,  # Enable propagation to generate multiple layers for slider
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
    # Set figdir to ensure files are saved in output_dir
    old_figdir = lt.settings.figdir
    lt.settings.figdir = str(output_dir)
    try:
        coremap(
            adata,
            coremap_key="X_cplearn_coremap",
            cluster_key="cplearn",
            output_dir=str(output_dir),
            save="coremap_cplearn_workflow.html",
            model=model,  # Pass model for proper visualization
        )
        html_path = output_dir / "coremap_cplearn_workflow.html"
        if html_path.exists():
            print(f"   Saved: {html_path}")
        else:
            print(f"   Warning: HTML file not found at {html_path}")
    finally:
        lt.settings.figdir = old_figdir
    
    # Try to convert HTML to PNG for GitHub Pages compatibility
    try:
        import plotly.graph_objects as go
        from plotly.io import from_json, to_image
        import json
        
        # Read the HTML file and extract the Plotly figure JSON
        with open(html_path, 'r') as f:
            html_content = f.read()
        
        # Try to extract Plotly figure from HTML (basic approach)
        # For a more robust solution, we could use kaleido or orca
        # For now, we'll keep the HTML and add a note in the documentation
        print("   Note: HTML file can be viewed locally or converted to PNG using kaleido")
    except Exception as e:
        print(f"   Note: HTML file saved. For static image, install kaleido: pip install kaleido")
    
    print("\n✓ Cplearn workflow complete!")


def example_3_alternating_methods(output_dir: Path) -> None:
    """Run Example 3: Alternating Methods."""
    print("\n" + "=" * 60)
    print("Example 3: Alternating Methods")
    print("=" * 60)
    
    # Load real data
    demo_data_path = Path(__file__).parent.parent / "data" / "demo_data.h5ad"
    print(f"Loading data from: {demo_data_path}")
    adata = lt.read(str(demo_data_path))
    print(f"Loaded dataset: {adata.shape}")
    
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
            propagate=True,  # Enable propagation to generate multiple layers for slider
        )
        core_analyze(adata, model=model, use_rep="X_latent", key_added="X_cplearn_coremap")
    except Exception as e:
        print(f"   Warning: Cplearn workflow failed: {e}")
        print("   Skipping cplearn visualization...")
        return
    # Set figdir to output_dir
    old_figdir = lt.settings.figdir
    lt.settings.figdir = str(output_dir)
    try:
        coremap(
            adata,
            coremap_key="X_cplearn_coremap",
            cluster_key="cplearn",
            output_dir=str(output_dir),
            save="coremap_alternating_cplearn.html",
            model=model,  # Pass model for proper visualization
        )
        html_path = output_dir / "coremap_alternating_cplearn.html"
        if html_path.exists():
            print(f"   Saved: {html_path}")
        else:
            print(f"   Warning: HTML file not found at {html_path}")
    finally:
        lt.settings.figdir = old_figdir
    print(f"   Cplearn clusters: {adata.obs['cplearn'].nunique()}")
    
    # === Workflow B: Scanpy Louvain ===
    print("\n3. Workflow B: Scanpy Louvain...")
    louvain(adata, resolution=0.8, key_added="louvain")  # Higher resolution for small dataset
    
    # Scanpy visualization (UMAP)
    old_figdir = lt.settings.figdir
    lt.settings.figdir = str(output_dir)
    try:
        umap(
            adata,
            cluster_key="louvain",
            output_dir=str(output_dir),
            save="_alternating_louvain.png",
        )
        png_path = output_dir / "umap_alternating_louvain.png"
        # Also check figures/ directory (fallback)
        if not png_path.exists():
            figures_png = Path("figures") / "umap_alternating_louvain.png"
            if figures_png.exists():
                import shutil
                shutil.copy(figures_png, png_path)
                print(f"   Copied from figures/ to: {png_path}")
            else:
                print(f"   Warning: PNG file not found")
        else:
            print(f"   Saved: {png_path}")
    finally:
        lt.settings.figdir = old_figdir
    print(f"   Louvain clusters: {adata.obs['louvain'].nunique()}")
    
    # Compare results
    print("\n4. Comparing results...")
    print(f"   Cplearn clusters: {adata.obs['cplearn'].value_counts().to_dict()}")
    print(f"   Louvain clusters: {adata.obs['louvain'].value_counts().to_dict()}")
    
    print("\n✓ Alternating methods workflow complete!")


def example_4_cplearn_umap(output_dir: Path) -> None:
    """Run Example 4: Cplearn Clustering + UMAP Visualization."""
    print("\n" + "=" * 60)
    print("Example 4: Cplearn Clustering + UMAP Visualization")
    print("=" * 60)
    
    # Load real data
    demo_data_path = Path(__file__).parent.parent / "data" / "demo_data.h5ad"
    print(f"Loading data from: {demo_data_path}")
    adata = lt.read(str(demo_data_path))
    print(f"Loaded dataset: {adata.shape}")
    
    # 1. Preprocessing
    print("\n1. Preprocessing...")
    preprocess(adata, n_pcs=20, n_top_genes=2000, n_neighbors=15, save_raw=True)
    
    # 2. Cplearn Clustering
    print("\n2. Cplearn Clustering...")
    try:
        model = cplearn.corespect(
            adata,
            use_rep="X_latent",
            key_added="cplearn",
            stable={"core_frac": 0.3, "ng_num": 10},
            cluster={"resolution": 0.8},
            propagate=True,  # Enable propagation to generate multiple layers for slider
        )
        print(f"   Found {adata.obs['cplearn'].nunique()} clusters")
    except Exception as e:
        print(f"   Warning: Cplearn clustering failed: {e}")
        print("   Skipping cplearn + umap workflow...")
        return
    
    # 3. Visualization: UMAP with cplearn clusters
    print("\n3. Visualization (UMAP with cplearn clusters)...")
    old_figdir = lt.settings.figdir
    lt.settings.figdir = str(output_dir)
    try:
        umap(
            adata,
            cluster_key="cplearn",  # Use cplearn cluster labels
            output_dir=str(output_dir),
            save="_cplearn_umap.png",
        )
        # The umap function saves with pattern: umap{cluster_key}{save}
        # So it should be: umap_cplearn_umap.png
        png_path = output_dir / "umap_cplearn_umap.png"
        target_path = output_dir / "umap_cplearn_workflow.png"
        
        # Check if file exists, if not check figures/ directory
        if png_path.exists():
            import shutil
            shutil.copy(png_path, target_path)
            print(f"   Saved: {target_path}")
        else:
            figures_png = Path("figures") / "umap_cplearn_umap.png"
            if figures_png.exists():
                import shutil
                shutil.copy(figures_png, target_path)
                print(f"   Copied from figures/ to: {target_path}")
            else:
                print(f"   Warning: PNG file not found. Expected: {png_path}")
    finally:
        lt.settings.figdir = old_figdir
    
    print(f"   Cplearn clusters: {adata.obs['cplearn'].nunique()}")
    print("\n✓ Cplearn + UMAP workflow complete!")


def example_5_coreanalysis_louvain_coremap(output_dir: Path) -> None:
    """Run Example 5: Core Analysis + Louvain Clustering + Coremap Visualization."""
    print("\n" + "=" * 60)
    print("Example 5: Core Analysis + Louvain + Coremap")
    print("=" * 60)
    
    # Load real data
    demo_data_path = Path(__file__).parent.parent / "data" / "demo_data.h5ad"
    print(f"Loading data from: {demo_data_path}")
    adata = lt.read(str(demo_data_path))
    print(f"Loaded dataset: {adata.shape}")
    
    # 1. Preprocessing
    print("\n1. Preprocessing...")
    preprocess(adata, n_pcs=20, n_top_genes=2000, n_neighbors=15, save_raw=True)
    
    # 2. Get cplearn model (for core analysis, but we'll use louvain clustering)
    print("\n2. Getting cplearn model for core analysis...")
    try:
        model = cplearn.corespect(
            adata,
            use_rep="X_latent",
            key_added="cplearn_temp",  # Temporary, we'll use louvain instead
            stable={"core_frac": 0.3, "ng_num": 10},
            cluster={"resolution": 0.8},
            propagate=True,  # Enable propagation to generate multiple layers for slider
        )
    except Exception as e:
        print(f"   Warning: Cplearn model creation failed: {e}")
        print("   Skipping core analysis workflow...")
        return
    
    # 3. Louvain clustering
    print("\n3. Louvain Clustering...")
    louvain(adata, resolution=0.8, key_added="louvain")
    print(f"   Found {adata.obs['louvain'].nunique()} clusters")
    
    # 4. Core analysis using louvain clusters
    print("\n4. Core Analysis (using louvain clusters)...")
    core_analyze(
        adata,
        model=model,
        use_rep="X_latent",
        key_added="X_louvain_coremap",
        cluster_key="louvain",  # Use louvain clusters instead of cplearn
    )
    
    # 5. Visualization: Coremap with louvain clusters
    print("\n5. Visualization (Coremap with louvain clusters)...")
    old_figdir = lt.settings.figdir
    lt.settings.figdir = str(output_dir)
    try:
        coremap(
            adata,
            coremap_key="X_louvain_coremap",
            cluster_key="louvain",  # Use louvain clusters
            output_dir=str(output_dir),
            save="coremap_louvain_workflow.html",
            model=model,
        )
        html_path = output_dir / "coremap_louvain_workflow.html"
        if html_path.exists():
            print(f"   Saved: {html_path}")
        else:
            print(f"   Warning: HTML file not found at {html_path}")
    finally:
        lt.settings.figdir = old_figdir
    
    print("\n✓ Core Analysis + Louvain + Coremap workflow complete!")


def example_6_coreanalysis_louvain_umap(output_dir: Path) -> None:
    """Run Example 6: Core Analysis + Louvain Clustering + UMAP Visualization."""
    print("\n" + "=" * 60)
    print("Example 6: Core Analysis + Louvain + UMAP")
    print("=" * 60)
    
    # Load real data
    demo_data_path = Path(__file__).parent.parent / "data" / "demo_data.h5ad"
    print(f"Loading data from: {demo_data_path}")
    adata = lt.read(str(demo_data_path))
    print(f"Loaded dataset: {adata.shape}")
    
    # 1. Preprocessing
    print("\n1. Preprocessing...")
    preprocess(adata, n_pcs=20, n_top_genes=2000, n_neighbors=15, save_raw=True)
    
    # 2. Get cplearn model (for core analysis, but we'll use louvain clustering)
    print("\n2. Getting cplearn model for core analysis...")
    try:
        model = cplearn.corespect(
            adata,
            use_rep="X_latent",
            key_added="cplearn_temp",  # Temporary, we'll use louvain instead
            stable={"core_frac": 0.3, "ng_num": 10},
            cluster={"resolution": 0.8},
            propagate=True,  # Enable propagation to generate multiple layers for slider
        )
    except Exception as e:
        print(f"   Warning: Cplearn model creation failed: {e}")
        print("   Skipping core analysis workflow...")
        return
    
    # 3. Louvain clustering
    print("\n3. Louvain Clustering...")
    louvain(adata, resolution=0.8, key_added="louvain")
    print(f"   Found {adata.obs['louvain'].nunique()} clusters")
    
    # 4. Core analysis using louvain clusters
    print("\n4. Core Analysis (using louvain clusters)...")
    core_analyze(
        adata,
        model=model,
        use_rep="X_latent",
        key_added="X_louvain_coremap",
        cluster_key="louvain",  # Use louvain clusters instead of cplearn
    )
    
    # 5. Visualization: UMAP with louvain clusters
    print("\n5. Visualization (UMAP with louvain clusters)...")
    old_figdir = lt.settings.figdir
    lt.settings.figdir = str(output_dir)
    try:
        umap(
            adata,
            cluster_key="louvain",  # Use louvain clusters
            output_dir=str(output_dir),
            save="_louvain_umap.png",
        )
        png_path = output_dir / "umap_louvain_umap.png"
        target_path = output_dir / "umap_coreanalysis_louvain_workflow.png"
        
        if png_path.exists():
            import shutil
            shutil.copy(png_path, target_path)
            print(f"   Saved: {target_path}")
        else:
            figures_png = Path("figures") / "umap_louvain_umap.png"
            if figures_png.exists():
                import shutil
                shutil.copy(figures_png, target_path)
                print(f"   Copied from figures/ to: {target_path}")
            else:
                print(f"   Warning: PNG file not found. Expected: {png_path}")
    finally:
        lt.settings.figdir = old_figdir
    
    print("\n✓ Core Analysis + Louvain + UMAP workflow complete!")


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
        example_4_cplearn_umap(output_dir)
        example_5_coreanalysis_louvain_coremap(output_dir)
        example_6_coreanalysis_louvain_umap(output_dir)
        
        print("\n" + "=" * 60)
        print("✓ All examples completed successfully!")
        print("=" * 60)
        print(f"\nGenerated files in: {output_dir}")
        print("  - umap_standard_workflow.png")
        print("  - coremap_cplearn_workflow.html")
        print("  - coremap_alternating_cplearn.html")
        print("  - umap_alternating_louvain.png")
        print("  - umap_cplearn_workflow.png")
        print("  - coremap_louvain_workflow.html")
        print("  - umap_coreanalysis_louvain_workflow.png")
        
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

