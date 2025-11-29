from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd
from anndata import AnnData

from lotus.methods.cplearn.external import cplearn
from lotus.methods.scanpy import tl as sc_tl


def pick_groups(labels: Sequence[int]) -> tuple[set[int], set[int]]:
    """
    Select two groups for differential expression analysis
    
    Parameters:
        labels: Sequence of cluster labels
    
    Returns:
        Two sets of group labels (groups_a, groups_b)
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


def rank_genes_groups(
    adata: AnnData,
    groupby: str,
    *,
    groups: str | Sequence[str] = "all",
    reference: str = "rest",
    n_genes: int | None = None,
    method: str | None = None,
    use_raw: bool | None = None,
    layer: str | None = None,
    key_added: str | None = None,
    **kwargs,
) -> None:
    """
    Rank genes for characterizing groups (differential expression analysis).
    
    This function performs differential expression analysis using scanpy's rank_genes_groups,
    which supports multiple statistical methods (wilcoxon, t-test, logreg, etc.) and can
    compare multiple groups simultaneously.
    
    Parameters:
        adata: AnnData object
        groupby: Key in adata.obs that contains group labels
        groups: Subset of groups to compare. Default: "all"
        reference: Reference group for comparison. Default: "rest" (compare each group to all others)
        n_genes: Number of genes to rank. Default: None (all genes)
        method: Statistical method to use. Options: "wilcoxon" (default), "t-test", "logreg", etc.
        use_raw: Whether to use raw data. Default: None (auto-detect)
        layer: Layer to use for analysis. Default: None (use X)
        key_added: Key name for results in adata.uns. Default: None ("rank_genes_groups")
        **kwargs: Additional arguments passed to scanpy's rank_genes_groups
    
    Returns:
        None (results stored in adata.uns)
    """
    sc_tl.rank_genes_groups(
        adata,
        groupby=groupby,
        groups=groups,
        reference=reference,
        n_genes=n_genes,
        method=method,
        use_raw=use_raw,
        layer=layer,
        key_added=key_added,
        **kwargs,
    )


def filter_rank_genes_groups(
    adata: AnnData,
    *,
    key: str | None = None,
    groupby: str | None = None,
    use_raw: bool | None = None,
    key_added: str = "rank_genes_groups_filtered",
    min_in_group_fraction: float = 0.25,
    min_fold_change: float = 1,
    max_out_group_fraction: float = 0.5,
    compare_abs: bool = False,
) -> None:
    """
    Filter ranked genes groups based on expression criteria.
    
    This function filters the results from rank_genes_groups based on expression
    fraction, fold change, and other criteria to identify high-quality marker genes.
    
    Parameters:
        adata: AnnData object
        key: Key in adata.uns containing rank_genes_groups results. Default: None ("rank_genes_groups")
        groupby: Key in adata.obs that contains group labels. Default: None (auto-detect from key)
        use_raw: Whether to use raw data. Default: None (auto-detect)
        key_added: Key name for filtered results in adata.uns. Default: "rank_genes_groups_filtered"
        min_in_group_fraction: Minimum fraction of cells expressing gene in group. Default: 0.25
        min_fold_change: Minimum fold change. Default: 1
        max_out_group_fraction: Maximum fraction of cells expressing gene outside group. Default: 0.5
        compare_abs: Whether to compare absolute values. Default: False
    
    Returns:
        None (filtered results stored in adata.uns)
    """
    sc_tl.filter_rank_genes_groups(
        adata,
        key=key,
        groupby=groupby,
        use_raw=use_raw,
        key_added=key_added,
        min_in_group_fraction=min_in_group_fraction,
        min_fold_change=min_fold_change,
        max_out_group_fraction=max_out_group_fraction,
        compare_abs=compare_abs,
    )


def marker_gene_overlap(
    adata: AnnData,
    reference_markers: dict[str, set[int] | list[int]],
    *,
    key: str = "rank_genes_groups",
    method: str = "overlap_count",
    normalize: str | None = None,
    top_n_markers: int | None = None,
    adj_pval_threshold: float | None = None,
    key_added: str = "marker_gene_overlap",
    inplace: bool = False,
) -> pd.DataFrame | None:
    """
    Analyze marker gene overlap between data and reference markers.
    
    This function compares marker genes identified in the data with reference marker
    genes (e.g., from literature or other datasets) to assess similarity and validate
    cell type annotations.
    
    Parameters:
        adata: AnnData object
        reference_markers: Dictionary mapping group names to sets/lists of marker gene names
        key: Key in adata.uns containing rank_genes_groups results. Default: "rank_genes_groups"
        method: Method for computing overlap. Options: "overlap_count", "jaccard", etc. Default: "overlap_count"
        normalize: Normalization method. Options: "reference", "data", or None. Default: None
        top_n_markers: Number of top markers to consider. Default: None (all)
        adj_pval_threshold: Adjusted p-value threshold for filtering. Default: None
        key_added: Key name for results in adata.uns. Default: "marker_gene_overlap"
        inplace: Whether to modify adata in place. Default: False
    
    Returns:
        pd.DataFrame | None: Overlap analysis results (if inplace=False)
    """
    return sc_tl.marker_gene_overlap(
        adata,
        reference_markers,
        key=key,
        method=method,
        normalize=normalize,
        top_n_markers=top_n_markers,
        adj_pval_threshold=adj_pval_threshold,
        key_added=key_added,
        inplace=inplace,
    )


__all__ = [
    "marker_genes",
    "run_differential_expression",
    "pick_groups",
    "rank_genes_groups",
    "filter_rank_genes_groups",
    "marker_gene_overlap",
]