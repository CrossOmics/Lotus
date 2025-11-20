from __future__ import annotations

from pathlib import Path
from typing import Sequence

from anndata import AnnData

import lotus as lt


def umap(
    adata: AnnData,
    *,
    cluster_key: str = "cplearn_labels",
    truth_key: str | None = "truth",
    output_dir: Path | None = None,
    save: str | bool = "_clusters.png",
    show: bool = False,
    compute_umap: bool = True,
    min_dist: float = 0.5,
    spread: float = 1.0,
    n_components: int = 2,
) -> None:
    """
    UMAP visualization step: Compute UMAP embedding and generate UMAP plot
    
    Parameters:
        adata: AnnData object
        cluster_key: Key name for cluster labels in adata.obs
        truth_key: Key name for truth labels in adata.obs. If None, will not be displayed
        output_dir: Output directory. If None, uses current figdir setting
        save: Save filename or False to not save
        show: Whether to show the plot
        compute_umap: Whether to compute UMAP embedding (if not already computed)
        min_dist: Minimum distance parameter for UMAP
        spread: Spread parameter for UMAP
        n_components: Number of components for UMAP
    """
    # Compute UMAP embedding if needed
    if compute_umap and "X_umap" not in adata.obsm:
        lt.tl.umap(
            adata,
            min_dist=min_dist,
            spread=spread,
            n_components=n_components,
        )
    
    if output_dir is not None:
        output_dir = Path(output_dir)
        # Scanpy requires the directory to exist before saving
        output_dir.mkdir(parents=True, exist_ok=True)
        old_figdir = lt.settings.figdir
        # Use absolute path to ensure scanpy can find the directory
        abs_output_dir = str(output_dir.resolve())
        lt.settings.figdir = abs_output_dir
        try:
            colors = [cluster_key]
            if truth_key and truth_key in adata.obs:
                colors.append(truth_key)
            
            lt.pl.umap(
                adata,
                color=colors,
                wspace=0.35,
                show=show,
                save=save,
            )
        finally:
            lt.settings.figdir = old_figdir
    else:
        colors = [cluster_key]
        if truth_key and truth_key in adata.obs:
            colors.append(truth_key)
        
        lt.pl.umap(
            adata,
            color=colors,
            wspace=0.35,
            show=show,
            save=save,
        )


def render_visualizations(
    adata: AnnData,
    marker_genes: Sequence[str],
    output_dir: Path,
    *,
    cluster_key: str = "cplearn_labels",
    truth_key: str | None = "truth",
) -> None:
    """
    Visualization module: Generate UMAP and marker gene visualizations
    
    Parameters:
        adata: AnnData object
        marker_genes: List of marker genes
        output_dir: Output directory
        cluster_key: Key name for cluster labels in adata.obs
        truth_key: Key name for truth labels in adata.obs. If None, will not be displayed
    """
    output_dir = Path(output_dir)
    # Scanpy requires the directory to exist before saving
    output_dir.mkdir(parents=True, exist_ok=True)
    
    old_figdir = lt.settings.figdir
    # Use absolute path to ensure scanpy can find the directory
    abs_output_dir = str(output_dir.resolve())
    lt.settings.figdir = abs_output_dir
    try:
        # UMAP visualization - pass output_dir=None to use current figdir setting
        umap(adata, cluster_key=cluster_key, truth_key=truth_key, output_dir=None, save="_clusters.png", show=False)
        
        # Violin plot for marker genes
        if marker_genes:
            groupby = truth_key if (truth_key and truth_key in adata.obs) else cluster_key
            # Suppress matplotlib categorical units info messages
            import logging
            import warnings
            # Filter both warnings and logging messages
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message=".*categorical units.*")
                # Temporarily filter logging for matplotlib
                matplotlib_logger = logging.getLogger("matplotlib")
                old_level = matplotlib_logger.level
                matplotlib_logger.setLevel(logging.WARNING)
                try:
                    lt.pl.violin(
                        adata,
                        keys=list(marker_genes),
                        groupby=groupby,
                        layer="raw_counts",
                        multi_panel=True,
                        use_raw=False,
                        show=False,
                        save="_markers.png",
                    )
                finally:
                    matplotlib_logger.setLevel(old_level)
    finally:
        lt.settings.figdir = old_figdir
