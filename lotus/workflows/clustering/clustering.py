from __future__ import annotations

from collections import Counter
from typing import Iterable, Literal

import numpy as np
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
    random_state: int = 0,
    neighbors_key: str | None = None,
    obsp: str | None = None,
    # common parameters
    print_summary: bool = True,
) -> cplearn.CorespectModel | None:
    """
    Perform clustering analysis.
    
    Parameters:
        adata (AnnData): AnnData object
        method (Literal["cplearn", "leiden", "louvain"]): Clustering method. Default: "leiden"
        use_rep (str | None): Representation to use. Default: None (auto-detect)
        key_added (str | None): Key name for cluster labels in adata.obs. Default: None (method-specific)
        cluster_resolution (float): Clustering resolution. Default: 1.2
        stable_core_frac (float): Stable core fraction (cplearn only). Default: 0.25
        stable_ng_num (int): Number of neighbors for stable core (cplearn only). Default: 8
        fine_grained (bool): Use fine-grained clustering (cplearn only). Default: False
        propagate (bool): Propagate labels (cplearn only). Default: True
        random_state (int): Random seed (scanpy methods only). Default: 0
        neighbors_key (str | None): Key in adata.uns for neighbors (scanpy methods only). Default: None
        obsp (str | None): Key in adata.obsp for adjacency matrix (scanpy methods only). Default: None
        print_summary (bool): Print cluster summary. Default: True
    
    Returns:
        None | cplearn.CorespectModel: None for scanpy methods, CorespectModel for cplearn
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
            "Using method='cplearn' in clustering() is deprecated. "
            "Please use cplearn API directly:\n"
            "  from lotus.methods.cplearn.external import cplearn\n"
            "  model = cplearn.corespect(adata, use_rep='X_pca', key_added='cplearn')\n"
            "For core layer analysis, use:\n"
            "  from lotus.workflows import core_analysis\n"
            "  core_analysis(adata, model=model)",
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
        
        # No need to copy to clustering_labels - use the key_added directly
        
        if print_summary:
            print(
                "Cluster summary:",
                summarize_clusters(adata.obs[key_added]),
            )
            if clustering_key in adata.obs:
                print(f"Stored clustering labels: `adata.obs['{clustering_key}']`")
        
        return None
    
    # Should not reach here, but handle unknown method
    raise ValueError(f"Unknown clustering method: {method}. Supported methods: 'leiden', 'louvain', 'cplearn'")


# Alias for backward compatibility
run_clustering = clustering
