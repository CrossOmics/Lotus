from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd
from anndata import AnnData

from lotus.methods.cplearn.external import cplearn


def pick_groups(labels: Sequence[int]) -> tuple[set[int], set[int]]:
    """
    Select two groups for differential expression analysis
    
    Parameters:
        labels: Sequence of cluster labels
    
    Returns:
        Two sets of group labels
    """
    unique = [label for label in sorted(set(labels)) if label != -1]
    if len(unique) < 2:
        raise ValueError("At least two non-negative cluster labels are required for differential expression analysis")
    return {unique[0]}, {unique[1]}


def marker_genes(
    adata: AnnData,
    *,
    cluster_key: str | None = None,
    groups_a: set[int] | None = None,
    groups_b: set[int] | None = None,
    layer: str | None = None,
    min_detect_pct: float = 0.0,
    min_cells_per_group: int = 5,
    auto_pick_groups: bool = True,
) -> pd.DataFrame:
    """
    Identify differentially expressed genes (marker genes) between groups.
    
    Parameters:
        adata (AnnData): AnnData object
        cluster_key (str | None): Key name for cluster labels in adata.obs. Default: None (auto-detect)
        groups_a (set[int] | None): Set of cluster labels for first group. Default: None
        groups_b (set[int] | None): Set of cluster labels for second group. Default: None
        layer (str | None): Layer to use for analysis. Default: None (auto-detect)
        min_detect_pct (float): Minimum detection percentage. Default: 0.0
        min_cells_per_group (int): Minimum number of cells per group. Default: 5
        auto_pick_groups (bool): Auto-select first two non-negative clusters. Default: True
    
    Returns:
        pd.DataFrame: Differential expression analysis results
    """
    # Auto-detect cluster key if not specified
    if cluster_key is None:
        for key in ["cplearn", "leiden", "louvain"]:
            if key in adata.obs:
                cluster_key = key
                break
        if cluster_key is None:
            raise ValueError(
                "No cluster key found. Please specify cluster_key or run clustering first. "
                "Compatible keys: 'cplearn', 'leiden', 'louvain'"
            )
    
    # Auto-detect layer if not specified
    if layer is None:
        for layer_name in ["raw_counts", "raw", "counts"]:
            if layer_name in adata.layers:
                layer = layer_name
                break
        # If no layer found, use None (will use adata.X)
        if layer is None:
            layer = None
    
    # Handle categorical cluster labels (scanpy format)
    labels_raw = adata.obs[cluster_key]
    if isinstance(labels_raw.dtype, pd.CategoricalDtype):
        # Convert categorical to numeric for comparison
        labels = np.asarray([int(x) if str(x).isdigit() else -1 for x in labels_raw])
    else:
        labels = np.asarray(labels_raw.astype(int))
    
    if auto_pick_groups and (groups_a is None or groups_b is None):
        groups_a, groups_b = pick_groups(labels)
    
    if groups_a is None or groups_b is None:
        raise ValueError("Must specify groups_a and groups_b, or set auto_pick_groups=True")
    
    groups_a = {int(x) for x in groups_a}
    groups_b = {int(x) for x in groups_b}
    print(f"Comparing groups {groups_a} vs {groups_b}")
    
    # Check if cplearn layers info exists, if not create default
    layers_key = f"{cluster_key}_cplearn"
    if layers_key not in adata.uns:
        # Create default layers structure (all cells in one layer)
        # This allows DEG analysis to work even without cplearn clustering
        print(f"[DEG] Creating default layers structure for {layers_key}")
        adata.uns[layers_key] = {
            'layers': [[i for i in range(adata.n_obs)]]
        }
    
    # For scanpy methods (leiden, louvain), labels are strings ('0', '1', etc.)
    # but de_from_adata expects integer labels. Create a temporary integer label column
    temp_key = f"{cluster_key}_int"
    needs_update = False
    
    if temp_key not in adata.obs:
        needs_update = True
    elif adata.obs[temp_key].dtype != 'int64':
        needs_update = True
    else:
        # Check if the integer labels match the current string labels
        # Convert current labels to integers and compare
        current_labels_int = np.asarray([int(x) if str(x).isdigit() else -1 for x in adata.obs[cluster_key]])
        stored_labels_int = adata.obs[temp_key].values
        
        # Check if they match (same length and same values)
        if len(current_labels_int) != len(stored_labels_int) or not np.array_equal(current_labels_int, stored_labels_int):
            needs_update = True
    
    if needs_update:
        # Convert labels to integers and store temporarily
        adata.obs[temp_key] = labels
        use_key = temp_key
    else:
        use_key = temp_key  # Use existing integer labels
    
    options = cplearn.DEOptions(
        min_detect_pct=min_detect_pct,
        min_cells_per_group=min_cells_per_group,
    )
    
    # Use the integer label key for de_from_adata
    result = cplearn.de_from_adata(
        adata,
        key=use_key,  # Use integer label key
        layers_key=layers_key,  # Explicitly pass layers_key to use original cluster_key
        groups_a=groups_a,
        groups_b=groups_b,
        layer=layer,
        options=options,
    )
    
    # Clean up temporary integer label column if we created it
    if use_key == temp_key and temp_key in adata.obs:
        # Keep it for now in case it's needed, but could remove it
        pass
    
    # Ensure result has the expected column names and order
    # cplearn.de_from_adata returns columns: gene, statistic, pvalue, log2fc, effect_r, 
    # mean_a, mean_b, pct_expr_a, pct_expr_b, z_score, p_adj, ...
    # Reorder to match expected format: gene, log2fc, z_score, pvalue, p_adj, mean_a, mean_b, pct_expr_a, pct_expr_b
    expected_cols = ['gene', 'log2fc', 'z_score', 'pvalue', 'p_adj', 'mean_a', 'mean_b', 'pct_expr_a', 'pct_expr_b']
    available_cols = [col for col in expected_cols if col in result.columns]
    other_cols = [col for col in result.columns if col not in expected_cols]
    
    # Reorder columns: expected columns first, then others
    result = result[available_cols + other_cols]
    
    return result


# Alias for backward compatibility
run_differential_expression = marker_genes

__all__ = [
    "marker_genes",
    "run_differential_expression",
]
