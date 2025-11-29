from __future__ import annotations

from pathlib import Path
from typing import Literal, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData
import matplotlib.pyplot as plt

import lotus as lt
from lotus.methods.cplearn.external import cplearn


def umap(
    adata: AnnData,
    *,
    cluster_key: str = "cplearn",
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


def coremap(
    adata: AnnData,
    *,
    coremap_key: str = "X_cplearn_coremap",
    cluster_key: str | None = None,
    core_layer_key: str | None = None,
    truth_key: str | None = "truth",
    output_dir: Path | None = None,
    save: str | bool = "_coremap.html",
    show: bool = False,
    model: cplearn.CorespectModel | None = None,
    use_webgl: bool = False,
    **kwargs,
) -> None:
    """
    CoreMap visualization: Visualize core map embedding from cplearn using cplearn's visualize_coremap
    
    This function visualizes the core map embedding computed by core_analyze(),
    which provides a special dimensionality reduction representation highlighting
    core cell states and cell relationships. It uses cplearn's own visualization
    function which provides interactive Plotly plots with layer sliders.
    
    Parameters:
        adata: AnnData object with coremap embedding in adata.obsm
        coremap_key: Key name for coremap embedding in adata.obsm.
                     Default is "X_cplearn_coremap"
        cluster_key: Key name for cluster labels in adata.obs.
                     If None, auto-detects: "cplearn" > "leiden" > "louvain"
        core_layer_key: Key name for core layer marker in adata.obs.
                        If None, auto-detects: "{coremap_key}_is_core"
        truth_key: Key name for truth labels in adata.obs. If None, will not be displayed
        output_dir: Output directory. If None, uses current figdir setting
        save: Save filename or False to not save. Default is "_coremap.html" (Plotly HTML format)
        show: Whether to show the plot
        model: CorespectModel object from cplearn clustering. If None, will try to reconstruct
               from adata.uns. Required for proper layer visualization.
        use_webgl: Whether to use WebGL for rendering (faster for large datasets)
        **kwargs: Additional arguments (currently unused, kept for compatibility)
    
    Examples:
        >>> # Basic usage with auto-detection
        >>> coremap(adata, model=model, output_dir="./results")
        
        >>> # Explicit keys
        >>> coremap(
        ...     adata,
        ...     coremap_key="X_cplearn_coremap",
        ...     cluster_key="cplearn",
        ...     model=model,
        ...     output_dir="./results"
        ... )
    """
    try:
        from cplearn.coremap.vizualizer import visualize_coremap
    except ImportError:
        import warnings
        warnings.warn(
            "cplearn.coremap.vizualizer.visualize_coremap not found. "
            "Falling back to basic matplotlib visualization.",
            UserWarning,
        )
        # Fallback to basic visualization (could implement matplotlib version here)
        raise ImportError("cplearn.visualize_coremap is required for coremap visualization")
    
    # Check if coremap embedding exists
    if coremap_key not in adata.obsm:
        raise ValueError(
            f"Coremap embedding '{coremap_key}' not found in adata.obsm. "
            f"Available keys: {list(adata.obsm.keys())}. "
            f"Please run core_analyze() first to compute coremap embedding."
        )
    
    embedding = adata.obsm[coremap_key]
    
    # Check if embedding has valid (non-NaN) values
    valid_mask = ~np.isnan(embedding).any(axis=1)
    n_valid = np.sum(valid_mask)
    
    if n_valid == 0:
        import warnings
        warnings.warn(
            f"Coremap embedding '{coremap_key}' contains no valid (non-NaN) values. "
            f"Cannot create visualization.",
            UserWarning,
        )
        return
    
    # Check for history (required for visualize_coremap)
    history_key = f"{coremap_key}_history"
    if history_key not in adata.uns:
        import warnings
        warnings.warn(
            f"Coremap history '{history_key}' not found in adata.uns. "
            f"Cannot use cplearn's visualize_coremap. "
            f"Please ensure core_analyze() was called with a valid model.",
            UserWarning,
        )
        return
    
    # Auto-detect cluster_key if not specified
    if cluster_key is None:
        for key in ["cplearn", "leiden", "louvain"]:
            if key in adata.obs:
                cluster_key = key
                break
    
    if cluster_key is None or cluster_key not in adata.obs:
        import warnings
        warnings.warn(
            f"No valid cluster_key found. Using default labels (all zeros).",
            UserWarning,
        )
        labels = None
    else:
        # Get labels as numpy array
        labels = np.asarray(adata.obs[cluster_key].values, dtype=int)
    
    # Get model and layers_ information
    if model is None:
        # Try to get from adata.uns
        cplearn_info_key = f"{coremap_key}_cplearn"
        if cplearn_info_key in adata.uns:
            # We need the actual model object, not just metadata
            import warnings
            warnings.warn(
                f"Model object not provided. Cannot reconstruct full coremap_obj. "
                f"Please provide the model parameter from cplearn.corespect().",
                UserWarning,
            )
            return
        else:
            raise ValueError(
                f"Model object is required for coremap visualization. "
                f"Please provide the CorespectModel from cplearn.corespect()."
            )
    
    # Get layers_ from model
    if not hasattr(model, "layers_") or model.layers_ is None:
        raise ValueError("Model does not have layers_ attribute. Cannot visualize coremap.")
    
    layers_ = model.layers_
    
    # Reconstruct label_dict from history
    history = adata.uns[history_key]
    label_dict = {}
    for round_id_str, round_data in history.items():
        round_id = int(round_id_str)
        coords = np.asarray(round_data["coordinates"], dtype=float)
        label_dict[round_id] = coords
    
    # Create a simple proxy object for coremap_obj
    class CoremapProxy:
        def __init__(self, label_dict, layers_, X):
            self.label_dict = label_dict
            self.layers_ = layers_
            self.X = X
    
    # Get X from adata
    X = adata.X
    if hasattr(X, "toarray"):  # Sparse matrix
        X = X.toarray()
    X = np.asarray(X, dtype=float)
    
    # Create proxy object
    coremap_obj = CoremapProxy(label_dict, layers_, X)
    
    # Call visualize_coremap
    fig = visualize_coremap(coremap_obj, labels=labels, use_webgl=use_webgl)
    
    # Set up output directory
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save figure
    if save:
        if output_dir is not None:
            save_path = output_dir / save
        else:
            save_path = Path(lt.settings.figdir) / save
        
        # Ensure save path has .html extension for Plotly
        if not str(save_path).endswith('.html'):
            save_path = Path(str(save_path).replace('.png', '.html').replace('.pdf', '.html'))
        
        fig.write_html(str(save_path))
        print(f"Saving coremap visualization to {save_path}")
    
    if show:
        fig.show()


def draw_graph(
    adata: AnnData,
    *,
    cluster_key: str = "cplearn",
    truth_key: str | None = "truth",
    output_dir: Path | None = None,
    save: str | bool = "_draw_graph.png",
    show: bool = False,
    compute_draw_graph: bool = True,
    layout: str = "fa",
    init_pos: str | bool | None = None,
    root: int | None = None,
    random_state: int = 0,
    neighbors_key: str | None = None,
    **kwargs,
) -> None:
    """
    Compute force-directed graph layout and generate plot.
    
    Parameters:
        adata: AnnData object
        cluster_key: Key name for cluster labels in adata.obs. Default: "cplearn"
        truth_key: Key name for truth labels in adata.obs. Default: None
        output_dir: Output directory. Default: None
        save: Save filename or False to not save. Default: "_draw_graph.png"
        show: Whether to show the plot. Default: False
        compute_draw_graph: Compute graph layout if not already computed. Default: True
        layout: Layout algorithm. Default: "fa"
        init_pos: Initial position. Default: None
        root: Root node. Default: None
        random_state: Random seed. Default: 0
        neighbors_key: Key for neighbors in adata.uns. Default: None
        **kwargs: Additional arguments passed to lt.tl.draw_graph
    """
    # Compute graph layout if needed
    basis_key = f"X_draw_graph_{layout}"
    if compute_draw_graph and basis_key not in adata.obsm:
        lt.tl.draw_graph(
            adata,
            layout=layout,
            init_pos=init_pos,
            root=root,
            random_state=random_state,
            neighbors_key=neighbors_key,
            **kwargs,
        )
    
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        old_figdir = lt.settings.figdir
        abs_output_dir = str(output_dir.resolve())
        lt.settings.figdir = abs_output_dir
        try:
            colors = [cluster_key]
            if truth_key and truth_key in adata.obs:
                colors.append(truth_key)
            
            lt.pl.draw_graph(
                adata,
                layout=layout,
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
        
        lt.pl.draw_graph(
            adata,
            layout=layout,
            color=colors,
            wspace=0.35,
            show=show,
            save=save,
        )


def visualization(
    adata: AnnData,
    *,
    method: Literal["umap", "draw_graph", "coremap"] = "umap",
    cluster_key: str = "cplearn",
    truth_key: str | None = "truth",
    output_dir: Path | None = None,
    save: str | bool | None = None,
    show: bool = False,
    **kwargs,
) -> None:
    """
    Unified visualization function supporting multiple embedding methods.
    
    Parameters:
        adata: AnnData object
        method: Visualization method. Default: "umap"
        cluster_key: Key name for cluster labels in adata.obs. Default: "cplearn"
        truth_key: Key name for truth labels in adata.obs. Default: None
        output_dir: Output directory. Default: None
        save: Save filename or False to not save. Default: None (method-specific)
        show: Whether to show the plot. Default: False
        **kwargs: Method-specific parameters:
            - umap: min_dist, spread, n_components, compute_umap
            - draw_graph: layout, init_pos, root, random_state, neighbors_key, compute_draw_graph, **kwargs
            - coremap: coremap_key, model, use_webgl, **kwargs
    """
    # Set default save filename based on method
    if save is None:
        if method == "umap":
            save = "_umap.png"
        elif method == "draw_graph":
            save = "_draw_graph.png"
        elif method == "coremap":
            save = "_coremap.html"
        else:
            save = f"_{method}.png"
    
    # Call method-specific function
    if method == "umap":
        umap(
            adata,
            cluster_key=cluster_key,
            truth_key=truth_key,
            output_dir=output_dir,
            save=save,
            show=show,
            **kwargs,
        )
    elif method == "draw_graph":
        draw_graph(
            adata,
            cluster_key=cluster_key,
            truth_key=truth_key,
            output_dir=output_dir,
            save=save,
            show=show,
            **kwargs,
        )
    elif method == "coremap":
        coremap(
            adata,
            cluster_key=cluster_key,
            truth_key=truth_key,
            output_dir=output_dir,
            save=save,
            show=show,
            **kwargs,
        )
    else:
        raise ValueError(
            f"Unknown visualization method: {method}. "
            f"Supported methods: 'umap', 'draw_graph', 'coremap'"
        )


def render_visualizations(
    adata: AnnData,
    marker_genes: Sequence[str],
    output_dir: Path,
    *,
    cluster_key: str = "cplearn",
    truth_key: str | None = "truth",
    include_coremap: bool = True,
    coremap_key: str = "X_cplearn_coremap",
    model: cplearn.CorespectModel | None = None,
) -> None:
    """
    Visualization module: Generate UMAP, CoreMap, and marker gene visualizations
    
    Parameters:
        adata: AnnData object
        marker_genes: List of marker genes
        output_dir: Output directory
        cluster_key: Key name for cluster labels in adata.obs
        truth_key: Key name for truth labels in adata.obs. If None, will not be displayed
        include_coremap: Whether to include coremap visualization if available
        coremap_key: Key name for coremap embedding in adata.obsm
        model: CorespectModel object from cplearn clustering. Required for coremap visualization.
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
        
        # CoreMap visualization (if available)
        if include_coremap and coremap_key in adata.obsm:
            try:
                coremap(
                    adata,
                    coremap_key=coremap_key,
                    cluster_key=cluster_key,
                    truth_key=truth_key,
                    output_dir=None,
                    save="_coremap.html",
                    show=False,
                    model=model,
                    use_webgl=False,
                )
            except Exception as e:
                import warnings
                warnings.warn(
                    f"Failed to generate coremap visualization: {type(e).__name__}: {e}",
                    UserWarning,
                )
        
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

__all__ = [
    "umap",
    "draw_graph",
    "coremap",
    "visualization",
    "render_visualizations",
]