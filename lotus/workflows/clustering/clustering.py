from __future__ import annotations

from collections import Counter
from typing import Iterable, Literal

import numpy as np
import pandas as pd
from anndata import AnnData

from lotus.methods.cplearn import cplearn
from lotus.methods.scanpy import tl as sc_tl


def summarize_clusters(labels: Iterable[int]) -> str:
    """Summarize cluster label statistics"""
    histogram = Counter(labels)
    pieces = [f"{label}: {count}" for label, count in sorted(histogram.items())]
    return ", ".join(pieces)


def cluster(
    adata: AnnData,
    *,
    method: Literal["cplearn", "leiden", "louvain"] = "leiden",
    use_rep: str | None = None,
    key_added: str | None = None,
    cluster_resolution: float = 1.2,
    # cplearn-specific parameters (deprecated, use cplearn API directly)
    stable_core_frac: float = 0.25,
    stable_ng_num: int = 8,
    fine_grained: bool = False,
    propagate: bool = True,
    # scanpy-specific parameters
    neighbors_key: str | None = None,
    obsp: str | None = None,
    # common parameters
    print_summary: bool = True,
) -> cplearn.CorespectModel | None:
    """
    Clustering step: Perform clustering analysis using scanpy methods (leiden, louvain)
    
    This function primarily implements scanpy clustering methods.
    For cplearn clustering, please use the cplearn API directly:
    
    .. code-block:: python
    
        from lotus.methods.cplearn import cplearn
        model = cplearn.corespect(adata, use_rep="X_pca", key_added="cplearn")
    
    Supports clustering methods:
    - "leiden" (default): Scanpy Leiden algorithm
    - "louvain": Scanpy Louvain algorithm
    - "cplearn": Deprecated - use cplearn API directly (see above)
    
    Compatible with scanpy workflow:
    - Accepts scanpy standard representations (X_pca, X_umap, etc.)
    - Works with scanpy neighbors graph (stored in adata.obsp)
    - Outputs cluster labels compatible with scanpy format
    
    Parameters:
        adata: AnnData object (compatible with scanpy AnnData)
        method: Clustering method to use. Options: "leiden" (default), "louvain", "cplearn" (deprecated)
        use_rep: Representation to use for clustering (deprecated for cplearn, use cplearn API directly)
        key_added: Key name for cluster labels in adata.obs.
                   If None, uses method-specific default: "leiden", "louvain", or "cplearn"
        cluster_resolution: Clustering resolution (applies to all methods)
        stable_core_frac: Stable core fraction (deprecated, use cplearn API directly)
        stable_ng_num: Number of neighbors for stable core (deprecated, use cplearn API directly)
        fine_grained: Whether to use fine-grained clustering (deprecated, use cplearn API directly)
        propagate: Whether to propagate labels (deprecated, use cplearn API directly)
        neighbors_key: Key in adata.uns where neighbors are stored (scanpy methods only)
        obsp: Key in adata.obsp to use as adjacency matrix (scanpy methods only)
        print_summary: Whether to print cluster summary
    
    Returns:
        None (for scanpy methods) or CorespectModel (for cplearn, deprecated)
    """
    # Set default key_added based on method
    if key_added is None:
        if method == "cplearn":
            key_added = "cplearn"
        elif method == "leiden":
            key_added = "leiden"
        elif method == "louvain":
            key_added = "louvain"
    
    # Handle cplearn method - redirect to cplearn API
    if method == "cplearn":
        import warnings
        warnings.warn(
            "Using method='cplearn' in cluster() is deprecated. "
            "Please use cplearn API directly:\n"
            "  from lotus.methods.cplearn import cplearn\n"
            "  model = cplearn.corespect(adata, use_rep='X_pca', key_added='cplearn')\n"
            "For core layer analysis, use:\n"
            "  from lotus.workflows import core_analyze\n"
            "  core_analyze(adata, model=model)",
            DeprecationWarning,
            stacklevel=2,
        )
        # Still support it for backward compatibility, but warn
        # Auto-detect representation if not specified
        if use_rep is None:
            if "X_latent" in adata.obsm:
                use_rep = "X_latent"
            elif "X_pca" in adata.obsm:
                use_rep = "X_pca"
            else:
                use_rep = "X"
        
        # Ensure neighbors graph exists (required for cplearn)
        if "neighbors" not in adata.uns:
            warnings.warn(
                "Neighbors graph not found. Please run preprocessing with neighbors() first, "
                "or use scanpy's sc.pp.neighbors() to compute neighbors graph.",
                UserWarning,
            )
        
        model = cplearn.corespect(
            adata,
            use_rep=use_rep,
            stable={
                "auto_select_core_frac": False,
                "core_frac": stable_core_frac,
                "ng_num": stable_ng_num,
                "densification": "k-nn",
            },
            cluster={
                "resolution": cluster_resolution,
                "auto_select_resolution": False,
                "densification": False,
            },
            fine_grained=fine_grained,
            propagate=propagate,
            key_added=key_added,
        )
        
        # Ensure output is compatible with scanpy format (categorical)
        cplearn_key = "cplearn"
        if cplearn_key in adata.obs:
            if not isinstance(adata.obs[cplearn_key].dtype, pd.CategoricalDtype):
                adata.obs[cplearn_key] = adata.obs[cplearn_key].astype("category")
        
        if print_summary:
            if cplearn_key in adata.obs:
                print(
                    "Cluster summary:",
                    summarize_clusters(adata.obs[cplearn_key]),
                )
        
        return model
    
    # Handle scanpy methods (leiden, louvain) - primary implementation
    if method in ("leiden", "louvain"):
        # Ensure neighbors graph exists (required for scanpy clustering)
        if "neighbors" not in adata.uns and neighbors_key is None:
            import warnings
            warnings.warn(
                "Neighbors graph not found. Please run preprocessing with neighbors() first, "
                "or use scanpy's sc.pp.neighbors() to compute neighbors graph.",
                UserWarning,
            )
        
        # Call scanpy clustering function
        if method == "leiden":
            sc_tl.leiden(
                adata,
                resolution=cluster_resolution,
                key_added=key_added,
                neighbors_key=neighbors_key,
                obsp=obsp,
            )
        elif method == "louvain":
            sc_tl.louvain(
                adata,
                resolution=cluster_resolution,
                key_added=key_added,
                neighbors_key=neighbors_key,
                obsp=obsp,
            )
        
        # Ensure output is categorical (scanpy already does this, but ensure it)
        if key_added in adata.obs:
            if not isinstance(adata.obs[key_added].dtype, pd.CategoricalDtype):
                adata.obs[key_added] = adata.obs[key_added].astype("category")
        
        # No need to copy to clustering_labels - use the key_added directly
        
        if print_summary:
            print(
                "Cluster summary:",
                summarize_clusters(adata.obs[key_added]),
            )
            if key_added in adata.obs:
                print(f"Stored clustering labels: `adata.obs['{key_added}']`")
        
        return None
    
    # Should not reach here, but handle unknown method
    raise ValueError(f"Unknown clustering method: {method}. Supported methods: 'leiden', 'louvain', 'cplearn'")


# Alias for backward compatibility
run_clustering = cluster

__all__ = [
    "cluster",
    "run_clustering",
]