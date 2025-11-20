from __future__ import annotations

from collections import Counter
from typing import Iterable

import pandas as pd
from anndata import AnnData

from lotus.methods.cplearn.external import cplearn


def summarize_clusters(labels: Iterable[int]) -> str:
    """Summarize cluster label statistics"""
    histogram = Counter(labels)
    pieces = [f"{label}: {count}" for label, count in sorted(histogram.items())]
    return ", ".join(pieces)


def clustering(
    adata: AnnData,
    *,
    use_rep: str | None = None,
    key_added: str = "cplearn_labels",
    stable_core_frac: float = 0.25,
    stable_ng_num: int = 8,
    cluster_resolution: float = 1.2,
    fine_grained: bool = False,
    propagate: bool = True,
    print_summary: bool = True,
) -> cplearn.CorespectModel:
    """
    Clustering step: Perform clustering analysis using cplearn
    
    Compatible with scanpy workflow:
    - Accepts scanpy standard representations (X_pca, X_umap, etc.)
    - Works with scanpy neighbors graph (stored in adata.obsp)
    - Outputs cluster labels compatible with scanpy format
    
    Parameters:
        adata: AnnData object (compatible with scanpy AnnData)
        use_rep: Representation to use for clustering. 
                 If None, auto-detects: "X_latent" > "X_pca" > "X"
                 Default is None (auto-detect)
        key_added: Key name for cluster labels in adata.obs
                   Compatible with scanpy keys like "leiden", "louvain"
        stable_core_frac: Stable core fraction
        stable_ng_num: Number of neighbors for stable core
        cluster_resolution: Clustering resolution
        fine_grained: Whether to use fine-grained clustering
        propagate: Whether to propagate labels
        print_summary: Whether to print cluster summary
    
    Returns:
        CorespectModel object
    """
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
