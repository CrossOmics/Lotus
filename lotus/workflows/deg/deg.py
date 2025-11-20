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
    Marker Genes step: Identify differentially expressed genes (marker genes) between groups
    
    Compatible with scanpy workflow:
    - Accepts scanpy cluster keys (e.g., "leiden", "louvain")
    - Works with scanpy layers (e.g., "raw", "counts", or custom layers)
    - Outputs DataFrame compatible with scanpy's rank_genes_groups format
    
    Parameters:
        adata: AnnData object (compatible with scanpy AnnData)
        cluster_key: Key name for cluster labels in adata.obs.
                    If None, auto-detects: "cplearn_labels" > "leiden" > "louvain"
                    Default is None (auto-detect)
        groups_a: Set of cluster labels for first group. 
                  If None and auto_pick_groups=True, will be auto-selected
        groups_b: Set of cluster labels for second group.
                  If None and auto_pick_groups=True, will be auto-selected
        layer: Layer to use for analysis.
               If None, auto-detects: "raw_counts" > "raw" > "counts" > None (use X)
               Default is None (auto-detect)
        min_detect_pct: Minimum detection percentage
        min_cells_per_group: Minimum number of cells per group
        auto_pick_groups: Whether to automatically select the first two non-negative clusters as comparison groups
    
    Returns:
        DataFrame with differential expression analysis results
        Compatible with scanpy's rank_genes_groups format
    """
    # Auto-detect cluster key if not specified
    if cluster_key is None:
        for key in ["cplearn_labels", "leiden", "louvain"]:
            if key in adata.obs:
                cluster_key = key
                break
        if cluster_key is None:
            raise ValueError(
                "No cluster key found. Please specify cluster_key or run clustering first. "
                "Compatible keys: 'cplearn_labels', 'leiden', 'louvain'"
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
    
    options = cplearn.DEOptions(
        min_detect_pct=min_detect_pct,
        min_cells_per_group=min_cells_per_group,
    )
    result = cplearn.de_from_adata(
        adata,
        groups_a=groups_a,
        groups_b=groups_b,
        layer=layer,
        options=options,
    )
    return result


# Alias for backward compatibility
run_differential_expression = marker_genes
