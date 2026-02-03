from __future__ import annotations

from typing import Any
import numpy as np
from anndata import AnnData
from cplearn.coremap import Coremap
from cplearn.coremap.vizualizer import visualize_coremap


def coremap(
        model: Any,
        adata: AnnData | None = None,
        *,
        # --- Coremap Configuration ---
        fast_view: bool = True,
        anchor_finding_mode: str = 'default',
        anchor_reach: int | None = None,
        # --- Visualization Parameters ---
        labels: np.ndarray | None = None,
        use_webgl: bool = True,
        # --- Standard Plotting Parameters ---
        show: bool = True,
        return_fig: bool = False,
) -> Any | None:
    """
    Visualize the CoreSpect analysis results using Coremap.

    This function wraps `cplearn.coremap.Coremap` and `visualize_coremap` to produce
    an interactive, hierarchical visualization of the core structure.

    Parameters:
        model: The trained CorespectModel object (returned by `lotus.tl.corespect` with copy=True
               or accessible via other means).
        adata: Optional AnnData object. If provided, the global UMAP embedding will be
               extracted from `adata.obsm['X_umap']` to ensure visual consistency.
        fast_view: If True, uses the existing global UMAP coordinates for layout.
                   If False, computes an anchored map layout (slower but better preserves hierarchy).
        anchor_finding_mode: Mode for finding anchors in fine-grained visualization ('default', etc.).
        anchor_reach: Reach parameter for anchor finding (radius).
        labels: Custom cluster labels. If None, uses `model.labels_` from the model.
        use_webgl: Whether to use WebGL for rendering. Recommended for large datasets (>10k cells).
        show: Whether to display the figure immediately.
        return_fig: Whether to return the Plotly figure object.

    Returns:
        If return_fig=True, returns the Plotly figure object.
        Otherwise, returns None.
    """

    # 1. Prepare Global Embedding (UMAP)
    # If adata is provided, we try to reuse its calculated UMAP to ensure
    # the Coremap visualization matches other standard Scanpy plots.
    global_umap = None
    if adata is not None:
        if "X_umap" in adata.obsm:
            global_umap = adata.obsm["X_umap"]
        else:
            print("Warning: 'X_umap' not found in adata.obsm. Coremap will compute its own UMAP layout.")

    # 2. Initialize Coremap Object
    # Bridges the raw CorespectModel data with the visualization layout
    cmap = Coremap(
        corespect_obj=model,
        global_umap=global_umap,
        fast_view=fast_view,
        anchor_finding_mode=anchor_finding_mode,
        anchor_reach=anchor_reach
    )

    # 3. Determine Labels
    # Default to the labels found by the model, but allow user override
    if labels is None:
        labels = model.labels_

    # 4. Generate Visualization
    # Creates the interactive Plotly figure
    fig = visualize_coremap(
        cmap,
        labels=labels,
        use_webgl=use_webgl
    )

    # 5. Output Handling
    if show:
        fig.show()

    if return_fig:
        return fig

    return None
