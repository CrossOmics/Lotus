from __future__ import annotations

import numpy as np
from anndata import AnnData

from lotus.methods.cplearn.external import cplearn


def core_selection(
    adata: AnnData,
    model: cplearn.CorespectModel,
    *,
    use_rep: str | None = None,
    key_added: str = "X_cplearn_coremap",
    print_summary: bool = True,
) -> None:
    """
    CoreSelection step: Compute core map embedding (Neighbors → CoreSelection)
    
    Compatible with scanpy workflow:
    - Works with scanpy neighbors graph (stored in adata.obsp)
    - Accepts scanpy standard representations (X_pca, X_umap, etc.)
    - Outputs embedding in adata.obsm compatible with scanpy format
    
    Note: Neighbors should be computed in preprocessing step before calling this function.
    Can use either lotus preprocessing or scanpy's sc.pp.neighbors().
    
    Parameters:
        adata: AnnData object (compatible with scanpy AnnData)
        model: CorespectModel object from clustering step
        use_rep: Representation to use for computation.
                 If None, auto-detects: "X_latent" > "X_pca" > "X"
                 Default is None (auto-detect)
        key_added: Key name for embedding results in adata.obsm
                   Compatible with scanpy embedding keys
        print_summary: Whether to print assignment summary
    """
    # Auto-detect representation if not specified
    if use_rep is None:
        if "X_latent" in adata.obsm:
            use_rep = "X_latent"
        elif "X_pca" in adata.obsm:
            use_rep = "X_pca"
        else:
            use_rep = "X"
    
    # Check for neighbors graph (required)
    if "neighbors" not in adata.uns:
        import warnings
        warnings.warn(
            "Neighbors graph not found. Please run preprocessing with neighbors() first, "
            "or use scanpy's sc.pp.neighbors() to compute neighbors graph.",
            UserWarning,
        )
    
    # Compute Anchored Map embedding
    cplearn.coremap_embedding(
        adata,
        model=model,
        use_rep=use_rep,
        key_added=key_added,
    )
    
    if print_summary:
        embedding = adata.obsm[key_added]
        assigned = np.sum(~np.isnan(embedding).any(axis=1))
        print(
            f"Stored anchored map embedding in `adata.obsm['{key_added}']` "
            f"({assigned}/{adata.n_obs} points assigned)."
        )


# Alias for backward compatibility
compute_coremap_embedding = core_selection
