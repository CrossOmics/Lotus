from __future__ import annotations

from collections import Counter
from typing import Iterable, Literal

import pandas as pd
from anndata import AnnData

from lotus.methods.cplearn.external import cplearn
from lotus.methods.scanpy import tl as sc_tl


def summarize_clusters(labels: Iterable[int]) -> str:
    """Summarize cluster label statistics"""
    histogram = Counter(labels)
    pieces = [f"{label}: {count}" for label, count in sorted(histogram.items())]
    return ", ".join(pieces)


def clustering(
    adata: AnnData,
    *,
    method: Literal["cplearn", "leiden", "louvain"] = "cplearn",
    use_rep: str | None = None,
    key_added: str | None = None,
    cluster_resolution: float = 1.2,
    # cplearn-specific parameters
    stable_core_frac: float = 0.25,
    stable_ng_num: int = 8,
    fine_grained: bool = False,
    propagate: bool = True,
    # scanpy-specific parameters
    random_state: int = 0,
    neighbors_key: str | None = None,
    obsp: str | None = None,
    # common parameters
    print_summary: bool = True,
) -> cplearn.CorespectModel | None:
    """
    Clustering step: Perform clustering analysis using cplearn, leiden, or louvain
    
    Supports multiple clustering methods:
    - "cplearn" (default): Lotus cplearn clustering algorithm
    - "leiden": Scanpy Leiden algorithm
    - "louvain": Scanpy Louvain algorithm
    
    Compatible with scanpy workflow:
    - Accepts scanpy standard representations (X_pca, X_umap, etc.)
    - Works with scanpy neighbors graph (stored in adata.obsp)
    - Outputs cluster labels compatible with scanpy format
    
    Parameters:
        adata: AnnData object (compatible with scanpy AnnData)
        method: Clustering method to use. Options: "cplearn" (default), "leiden", "louvain"
        use_rep: Representation to use for clustering (cplearn only).
                 If None, auto-detects: "X_latent" > "X_pca" > "X"
                 Default is None (auto-detect)
        key_added: Key name for cluster labels in adata.obs.
                   If None, uses method-specific default: "cplearn_labels", "leiden", or "louvain"
        cluster_resolution: Clustering resolution (applies to all methods)
        stable_core_frac: Stable core fraction (cplearn only)
        stable_ng_num: Number of neighbors for stable core (cplearn only)
        fine_grained: Whether to use fine-grained clustering (cplearn only)
        propagate: Whether to propagate labels (cplearn only)
        random_state: Random seed for reproducibility (scanpy methods only)
        neighbors_key: Key in adata.uns where neighbors are stored (scanpy methods only)
        obsp: Key in adata.obsp to use as adjacency matrix (scanpy methods only)
        print_summary: Whether to print cluster summary
    
    Returns:
        CorespectModel object (for cplearn) or None (for scanpy methods)
    """
    # Set default key_added based on method
    if key_added is None:
        if method == "cplearn":
            key_added = "cplearn_labels"
        elif method == "leiden":
            key_added = "leiden"
        elif method == "louvain":
            key_added = "louvain"
    
    # Handle scanpy methods (leiden, louvain)
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
                random_state=random_state,
                neighbors_key=neighbors_key,
                obsp=obsp,
            )
        elif method == "louvain":
            sc_tl.louvain(
                adata,
                resolution=cluster_resolution,
                key_added=key_added,
                random_state=random_state,
                neighbors_key=neighbors_key,
                obsp=obsp,
            )
        
        # Ensure output is categorical (scanpy already does this, but ensure it)
        if key_added in adata.obs:
            if not isinstance(adata.obs[key_added].dtype, pd.CategoricalDtype):
                adata.obs[key_added] = adata.obs[key_added].astype("category")
        
        if print_summary:
            print(
                "Cluster summary:",
                summarize_clusters(adata.obs[key_added]),
            )
        
        return None
    
    # Handle cplearn method (default)
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
        import warnings
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
    if key_added in adata.obs:
        if not isinstance(adata.obs[key_added].dtype, pd.CategoricalDtype):
            adata.obs[key_added] = adata.obs[key_added].astype("category")
    
    if print_summary:
        print(
            "Cluster summary:",
            summarize_clusters(adata.obs[key_added]),
        )
    
    return model


# Alias for backward compatibility
run_clustering = clustering
