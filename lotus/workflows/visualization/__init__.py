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
    
    Returns:
        None. Updates `adata.obsm` with UMAP embedding in `adata.obsm['X_umap']` if `compute_umap=True` and embedding not already present.
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
    fast_view: bool = False,
    **kwargs,
) -> None:
    """
    CoreMap visualization: Visualize core map embedding from cplearn using cplearn's visualize_coremap
    
    This function visualizes the core map embedding computed by core_analyze(),
    which provides a special dimensionality reduction representation highlighting
    core cell states and cell relationships. It uses cplearn's own visualization
    function which provides interactive Plotly plots with layer sliders (sidebar).
    
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
               from adata.uns. Required for proper layer visualization and sidebar (layer slider).
        use_webgl: Whether to use WebGL for rendering (faster for large datasets)
        fast_view: Whether to use a faster, less interactive view for large datasets. Default: False.
        **kwargs: Additional arguments (currently unused, kept for compatibility)
    
    Returns:
        None. Generates interactive Plotly visualization (saved to file if `save` is specified).
    
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
    
    # Step 1: Generate UMAP skeleton if not exists (same as API version)
    if 'X_umap' not in adata.obsm or adata.obsm['X_umap'].shape[1] != 2:
        print("[COREMAP] Computing UMAP skeleton...")
        # Use the representation for UMAP (get from kwargs if provided)
        use_rep = kwargs.get('use_rep', None)
        if use_rep and use_rep in adata.obsm:
            X_for_umap = adata.obsm[use_rep]
        elif 'X_pca' in adata.obsm:
            X_for_umap = adata.obsm['X_pca']
        elif 'X_latent' in adata.obsm:
            X_for_umap = adata.obsm['X_latent']
        else:
            X_for_umap = adata.X
            if hasattr(X_for_umap, 'toarray'):
                X_for_umap = X_for_umap.toarray()
        
        try:
            import umap
            reducer = umap.UMAP(n_components=2)
            X_umap = reducer.fit_transform(X_for_umap)
            adata.obsm['X_umap_skeleton'] = X_umap
            print(f"[COREMAP] UMAP skeleton computed: shape {X_umap.shape}")
        except ImportError:
            raise ImportError("UMAP is required for coremap visualization with sidebar. Please install umap-learn.")
    else:
        X_umap = adata.obsm['X_umap']
        print("[COREMAP] Using existing UMAP embedding")
    
    print("[COREMAP] Initializing Coremap...")
    try:
        from cplearn.coremap import Coremap
        cmap = Coremap(model, global_umap=X_umap, fast_view=fast_view)
    except Exception as e:
        raise RuntimeError(
            f"Failed to initialize Coremap object: {str(e)}. "
            f"This is required for sidebar (layer slider) functionality."
        ) from e
    
    # Step 3: Prepare labels for visualization (same as API version)
    labels_to_use = None
    if cluster_key and cluster_key in adata.obs:
        labels_to_use = np.asarray(adata.obs[cluster_key].values, dtype=int)
        print(f"[COREMAP] Using cluster labels from '{cluster_key}'")
    elif hasattr(model, 'labels_') and model is not None:
        labels_to_use = np.asarray(model.labels_, dtype=int)
        print("[COREMAP] Using labels from model")
    else:
        print("[COREMAP] Warning: No labels found, visualization may not show clusters")
    
    # Step 4: Generate visualization (same as API version)
    print("[COREMAP] Generating layer-wise visualization...")
    fig = visualize_coremap(cmap, labels_to_use, use_webgl=use_webgl)
    
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
        
        # Use same method as API version
        fig.write_html(str(save_path))
        print(f"Saving coremap visualization to {save_path}")
    
    if show:
        fig.show()


def tsne(
    adata: AnnData,
    *,
    cluster_key: str = "cplearn",
    truth_key: str | None = "truth",
    output_dir: Path | None = None,
    save: str | bool = "_tsne.png",
    show: bool = False,
    compute_tsne: bool = True,
    n_pcs: int | None = None,
    use_rep: str | None = None,
    perplexity: float = 30,
    random_state: int = 0,
    early_exaggeration: float = 12,
    learning_rate: float = 1000,
) -> None:
    """
    Compute t-SNE embedding and generate t-SNE plot.
    
    Parameters:
        adata: AnnData object
        cluster_key: Key name for cluster labels in adata.obs. Default: "cplearn"
        truth_key: Key name for truth labels in adata.obs. Default: None
        output_dir: Output directory. Default: None
        save: Save filename or False to not save. Default: "_tsne.png"
        show: Whether to show the plot. Default: False
        compute_tsne: Compute t-SNE embedding if not already computed. Default: True
        n_pcs: Number of PCs to use. Default: None
        use_rep: Representation to use. Default: None
        perplexity: Perplexity parameter. Default: 30
        random_state: Random seed. Default: 0
        early_exaggeration: Early exaggeration parameter. Default: 12
        learning_rate: Learning rate. Default: 1000
    
    Returns:
        None. Updates `adata.obsm` with t-SNE embedding in `adata.obsm['X_tsne']` if `compute_tsne=True` and embedding not already present.
    """
    # Compute t-SNE embedding if needed
    if compute_tsne and "X_tsne" not in adata.obsm:
        lt.tl.tsne(
            adata,
            n_pcs=n_pcs,
            use_rep=use_rep,
            perplexity=perplexity,
            random_state=random_state,
            early_exaggeration=early_exaggeration,
            learning_rate=learning_rate,
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
            
            lt.pl.tsne(
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
        
        lt.pl.tsne(
            adata,
            color=colors,
            wspace=0.35,
            show=show,
            save=save,
        )


def diffmap(
    adata: AnnData,
    *,
    cluster_key: str = "cplearn",
    truth_key: str | None = "truth",
    output_dir: Path | None = None,
    save: str | bool = "_diffmap.png",
    show: bool = False,
    compute_diffmap: bool = True,
    n_comps: int = 15,
    neighbors_key: str | None = None,
    random_state: int = 0,
) -> None:
    """
    Compute diffusion map embedding and generate diffusion map plot.
    
    Parameters:
        adata: AnnData object
        cluster_key: Key name for cluster labels in adata.obs. Default: "cplearn"
        truth_key: Key name for truth labels in adata.obs. Default: None
        output_dir: Output directory. Default: None
        save: Save filename or False to not save. Default: "_diffmap.png"
        show: Whether to show the plot. Default: False
        compute_diffmap: Compute diffusion map embedding if not already computed. Default: True
        n_comps: Number of components. Default: 15
        neighbors_key: Key for neighbors in adata.uns. Default: None
        random_state: Random seed. Default: 0
    
    Returns:
        None. Updates `adata.obsm` with diffusion map embedding in `adata.obsm['X_diffmap']` if `compute_diffmap=True` and embedding not already present.
    """
    # Compute diffusion map embedding if needed
    if compute_diffmap and "X_diffmap" not in adata.obsm:
        lt.tl.diffmap(
            adata,
            n_comps=n_comps,
            neighbors_key=neighbors_key,
            random_state=random_state,
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
            
            lt.pl.diffmap(
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
        
        lt.pl.diffmap(
            adata,
            color=colors,
            wspace=0.35,
            show=show,
            save=save,
        )


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
    
    Returns:
        None. Updates `adata.obsm` with graph layout in `adata.obsm[f'X_draw_graph_{layout}']` if `compute_draw_graph=True` and layout not already present.
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
    method: Literal["umap", "tsne", "diffmap", "draw_graph", "coremap"] = "umap",
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
            - tsne: n_pcs, use_rep, perplexity, random_state, early_exaggeration, learning_rate, compute_tsne
            - diffmap: n_comps, neighbors_key, random_state, compute_diffmap
            - draw_graph: layout, init_pos, root, random_state, neighbors_key, compute_draw_graph, **kwargs
            - coremap: coremap_key, model, use_webgl, **kwargs
    
    Returns:
        None. Updates `adata.obsm` with embedding (method-specific key) if computed, and generates visualization plots.
    """
    # Set default save filename based on method
    if save is None:
        if method == "umap":
            save = "_umap.png"
        elif method == "tsne":
            save = "_tsne.png"
        elif method == "diffmap":
            save = "_diffmap.png"
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
    elif method == "tsne":
        tsne(
            adata,
            cluster_key=cluster_key,
            truth_key=truth_key,
            output_dir=output_dir,
            save=save,
            show=show,
            **kwargs,
        )
    elif method == "diffmap":
        diffmap(
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
            f"Supported methods: 'umap', 'tsne', 'diffmap', 'draw_graph', 'coremap'"
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
    
    Returns:
        None. Generates and saves visualization plots (UMAP, CoreMap, marker gene violin plots) to `output_dir`.
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
    "tsne",
    "diffmap",
    "draw_graph",
    "coremap",
    "visualization",
    "render_visualizations",
]