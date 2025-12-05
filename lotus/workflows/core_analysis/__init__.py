from __future__ import annotations

from typing import Any, Mapping

import numpy as np
import pandas as pd
from anndata import AnnData

from lotus.methods.cplearn import cplearn


def core_analyze(
    adata: AnnData,
    model: cplearn.CorespectModel | None = None,
    *,
    use_rep: str | None = None,
    key_added: str = "X_cplearn_coremap",
    cluster_key: str | None = None,
    print_summary: bool = True,
    # corespect parameters (only used when model=None)
    flowrank: Mapping[str, Any] | None = None,
    stable: Mapping[str, Any] | None = None,
    fine: Mapping[str, Any] | None = None,
    cluster: Mapping[str, Any] | None = None,
    propagate_cfg: Mapping[str, Any] | None = None,
    fine_grained: bool = False,
    propagate: bool = True,
    cplearn_key: str = "cplearn",
    **kwargs,
) -> cplearn.CorespectModel | None:
    """
    CoreAnalysis step: Compute core map embedding and extract core layer information
    
    This function computes core map embedding and stores it in adata.obsm and extracts and saves core layer information (core_layers_) for visualization
    In addition, it is compatible with scanpy workflow because it works with scanpy neighbors graph (stored in adata.obsp),
    accepts scanpy standard representations (X_pca, X_umap, etc.), outputs embedding in adata.obsm compatible with scanpy format and saves core layer information in adata.obs for scanpy plotting
    
    If `model` is None, automatically calls `cplearn.corespect()` internally to generate the model.
    This allows for a simplified workflow: core_analyze -> clustering.
    
    Parameters:
        adata: AnnData object (compatible with scanpy AnnData)
        model: CorespectModel object from clustering step. If None, will automatically
               call `cplearn.corespect()` internally.
        use_rep: Representation to use for computation.
                 If None, auto-detects: "X_latent" > "X_pca" > "X"
                 Default is None (auto-detect)
        key_added: Key name for embedding results in adata.obsm
                   Compatible with scanpy embedding keys
        cluster_key: Key name for cluster labels in adata.obs.
                     If None, auto-detects: "cplearn" > "leiden" > "louvain"
                     Default is None (auto-detect)
        print_summary: Whether to print assignment summary
        
        # corespect parameters (only used when model=None)
        flowrank: FlowRank configuration dict
        stable: Stable core configuration dict
        fine: Fine-grained configuration dict
        cluster: Clustering configuration dict
        propagate_cfg: Propagation configuration dict
        fine_grained: Whether to run fine-grained refinement
        propagate: Whether to run propagation (default: True for slider support)
        cplearn_key: Key name for cplearn labels in adata.obs (default: "cplearn")
        **kwargs: Additional arguments passed to cplearn.corespect()
    
    Returns:
        The CorespectModel object (useful for visualization with coremap())
    """
    # Auto-detect representation if not specified
    if use_rep is None:
        if "X_latent" in adata.obsm:
            use_rep = "X_latent"
        elif "X_pca" in adata.obsm:
            use_rep = "X_pca"
        else:
            use_rep = "X"
    
    # If model is None, automatically call corespect to generate model
    if model is None:
        if print_summary:
            print("Running cplearn.corespect() to generate clustering model...")
        
        model = cplearn.corespect(
            adata,
            use_rep=use_rep,
            flowrank=flowrank,
            stable=stable,
            fine=fine,
            cluster=cluster,
            propagate_cfg=propagate_cfg,
            fine_grained=fine_grained,
            propagate=propagate,
            key_added=cplearn_key,
            **kwargs,
        )
        
        # Save marker: core_analyze has been done
        if 'core_analyze' not in adata.uns:
            adata.uns['core_analyze'] = {}
        
        adata.uns['core_analyze']['done'] = True
        adata.uns['core_analyze']['cplearn_clustering_done'] = True
        adata.uns['core_analyze']['cplearn_key'] = cplearn_key
        adata.uns['core_analyze']['use_rep'] = use_rep
        adata.uns['core_analyze']['stable'] = stable or {}
        adata.uns['core_analyze']['cluster'] = cluster or {}
        
        if print_summary:
            print(f"Saved core_analyze marker in adata.uns['core_analyze']")
    
    # Auto-detect cluster key if not specified
    if cluster_key is None:
        for key in ["cplearn", "leiden", "louvain"]:
            if key in adata.obs:
                cluster_key = key
                break
        # If not found, try to get from model
        if cluster_key is None and hasattr(model, "labels_"):
            cluster_key = "cplearn"
            # Save labels from model if not already saved
            labels = np.asarray(model.labels_, dtype=int)
            adata.obs[cluster_key] = pd.Categorical(labels)
    
    # Check for neighbors graph (required)
    if "neighbors" not in adata.uns:
        import warnings
        warnings.warn(
            "Neighbors graph not found. Please run preprocessing with neighbors() first, "
            "or use scanpy's sc.pp.neighbors() to compute neighbors graph.",
            UserWarning,
        )
    
    # Check if model has valid clustering results (at least 2 clusters)
    has_valid_clusters = False
    if hasattr(model, "labels_"):
        unique_labels = np.unique(model.labels_)
        # Filter out -1 (unassigned) labels
        valid_labels = unique_labels[unique_labels >= 0]
        has_valid_clusters = len(valid_labels) >= 2
    
    # Compute Anchored Map embedding (only if we have valid clusters)
    if has_valid_clusters:
        try:
            cplearn.coremap_embedding(
                adata,
                model=model,
                use_rep=use_rep,
                key_added=key_added,
            )
        except Exception as e:
            import warnings
            warnings.warn(
                f"Failed to compute coremap embedding: {type(e).__name__}: {e}. "
                f"Skipping embedding computation but will still save core layer information.",
                UserWarning,
            )
            # Create empty embedding as placeholder
            adata.obsm[key_added] = np.full((adata.n_obs, 2), np.nan, dtype=np.float32)
    else:
        import warnings
        warnings.warn(
            f"Model has insufficient clusters ({len(np.unique(model.labels_)) if hasattr(model, 'labels_') else 0}). "
            f"Skipping coremap embedding computation but will still save core layer information.",
            UserWarning,
        )
        # Create empty embedding as placeholder
        adata.obsm[key_added] = np.full((adata.n_obs, 2), np.nan, dtype=np.float32)
    
    # Extract and save core layer information
    if hasattr(model, "core_layers_") and model.core_layers_ is not None:
        # Get core layer indices (Layer 0 = core layer)
        core_layers = model.core_layers_
        if len(core_layers) > 0:
            # Flatten all core layer indices into a single set
            core_indices = set()
            for core_layer in core_layers:
                core_layer_arr = np.array(core_layer).flatten()
                core_indices.update(core_layer_arr.astype(int))
            
            # Create is_core marker in obs for easy visualization
            is_core = np.zeros(adata.n_obs, dtype=bool)
            is_core[list(core_indices)] = True
            adata.obs[f"{key_added}_is_core"] = pd.Categorical(is_core)
            
            # Also save core layer information in uns for detailed access
            if f"{key_added}_cplearn" not in adata.uns:
                adata.uns[f"{key_added}_cplearn"] = {}
            
            adata.uns[f"{key_added}_cplearn"]["core_layers"] = [
                np.array(layer).astype(int).tolist() 
                for layer in core_layers
            ]
            adata.uns[f"{key_added}_cplearn"]["core_indices"] = sorted(list(core_indices))
            adata.uns[f"{key_added}_cplearn"]["n_core_cells"] = len(core_indices)
    
    # Ensure cplearn labels are saved (if from model)
    cplearn_key = "cplearn"
    if cplearn_key not in adata.obs and hasattr(model, "labels_"):
        labels = np.asarray(model.labels_, dtype=int)
        adata.obs[cplearn_key] = pd.Categorical(labels)
    
    # Save cluster key info for reference
    if cluster_key:
        if f"{key_added}_cplearn" not in adata.uns:
            adata.uns[f"{key_added}_cplearn"] = {}
        adata.uns[f"{key_added}_cplearn"]["cluster_key"] = cluster_key
        adata.uns[f"{key_added}_cplearn"]["cplearn_key"] = cplearn_key
    
    if print_summary:
        embedding = adata.obsm[key_added]
        assigned = np.sum(~np.isnan(embedding).any(axis=1))
        print(
            f"Stored anchored map embedding in `adata.obsm['{key_added}']` "
            f"({assigned}/{adata.n_obs} points assigned)."
        )
        
        # Print core layer summary
        if f"{key_added}_is_core" in adata.obs:
            n_core = np.sum(adata.obs[f"{key_added}_is_core"] == True)
            print(
                f"Stored core layer information: "
                f"`adata.obs['{key_added}_is_core']` ({n_core}/{adata.n_obs} core cells)."
            )
        
        if cluster_key:
            print(f"Using cluster labels from `adata.obs['{cluster_key}']` for visualization.")
        
        # Print cplearn labels info if available
        if cplearn_key in adata.obs:
            n_clusters = adata.obs[cplearn_key].nunique()
            print(f"Stored cplearn labels: `adata.obs['{cplearn_key}']` ({n_clusters} clusters).")
    
    return model


__all__ = [
    "core_analyze",
]