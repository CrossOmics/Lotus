from typing import Literal, Sequence, Union, Iterable, Any, Optional
import pandas as pd
import scanpy as sc
from anndata import AnnData


def rank_genes_groups(
        adata: AnnData,
        groupby: str,
        *,
        mask_var: Any | None = None,
        use_raw: bool | None = None,
        groups: Literal['all'] | Iterable[str] = 'all',
        reference: str = 'rest',
        n_genes: int | None = None,
        rankby_abs: bool = False,
        pts: bool = False,
        key_added: str | None = None,
        copy: bool = False,
        method: Literal['logreg', 't-test', 'wilcoxon', 't-test_overestim_var'] | None = None,
        corr_method: Literal['benjamini-hochberg', 'bonferroni'] = 'benjamini-hochberg',
        tie_correct: bool = False,
        layer: str | None = None,
        **kwds,
) -> AnnData | None:
    """
    Rank genes for characterizing groups.

    Wraps scanpy.tl.rank_genes_groups.

    Parameters:
        adata: Annotated data matrix.
        groupby: The key of the observations grouping to consider.
        mask_var: Select subset of genes to use in statistical tests.
        use_raw: Use raw attribute of adata if present.
        groups: Subset of groups to which comparison shall be restricted.
        reference: If 'rest', compare each group to the union of the rest.
        n_genes: The number of genes that appear in the returned tables.
        method: Statistical method ('logreg', 't-test', 'wilcoxon', 't-test_overestim_var').
        corr_method: p-value correction method.
        pts: Compute the fraction of cells expressing the genes.
        key_added: The key in adata.uns information is saved to.
        copy: Whether to copy adata or modify it inplace.

    Returns:
        Returns None if copy=False, else returns an AnnData object.
    """
    return sc.tl.rank_genes_groups(
        adata,
        groupby,
        mask_var=mask_var,
        use_raw=use_raw,
        groups=groups,
        reference=reference,
        n_genes=n_genes,
        rankby_abs=rankby_abs,
        pts=pts,
        key_added=key_added,
        copy=copy,
        method=method,
        corr_method=corr_method,
        tie_correct=tie_correct,
        layer=layer,
        **kwds
    )


def filter_rank_genes_groups(
        adata: AnnData,
        *,
        key: str | None = None,
        groupby: str | None = None,
        use_raw: bool | None = None,
        key_added: str = 'rank_genes_groups_filtered',
        min_in_group_fraction: float = 0.25,
        min_fold_change: float = 1.0,
        max_out_group_fraction: float = 0.5,
        compare_abs: bool = False,
) -> None:
    """
    Filter out genes based on fold change and fraction of cells expressing.

    Wraps scanpy.tl.filter_rank_genes_groups.

    Parameters:
        adata: Annotated data matrix.
        key: Key in adata.uns where rank_genes_groups result is stored.
        groupby: The key of the observations grouping.
        key_added: Key to store filtered results.
        min_in_group_fraction: Min fraction of cells expressing gene in group.
        min_fold_change: Min fold change.
        max_out_group_fraction: Max fraction of cells expressing gene out of group.
        compare_abs: Compare absolute values of log fold change.

    Returns:
        None.
    """
    sc.tl.filter_rank_genes_groups(
        adata,
        key=key,
        groupby=groupby,
        use_raw=use_raw,
        key_added=key_added,
        min_in_group_fraction=min_in_group_fraction,
        min_fold_change=min_fold_change,
        max_out_group_fraction=max_out_group_fraction,
        compare_abs=compare_abs
    )


def marker_gene_overlap(
        adata: AnnData,
        reference_markers: dict[str, set] | dict[str, list],
        *,
        key: str = 'rank_genes_groups',
        method: Literal['overlap_count', 'overlap_coef', 'jaccard'] = 'overlap_count',
        normalize: Literal['reference', 'data'] | None = None,
        top_n_markers: int | None = None,
        adj_pval_threshold: float | None = None,
        key_added: str = 'marker_gene_overlap',
        inplace: bool = False,
) -> pd.DataFrame | AnnData | None:
    """
    Calculate an overlap score between data-derived marker genes and provided markers.

    Wraps scanpy.tl.marker_gene_overlap.

    Parameters:
        adata: Annotated data matrix.
        reference_markers: Dictionary of marker genes (e.g. {'CellType': {'GeneA', 'GeneB'}}).
        key: Key in adata.uns where rank_genes_groups output is stored.
        method: Method to calculate overlap ('overlap_count', 'overlap_coef', 'jaccard').
        normalize: Normalization option for overlap_count.
        top_n_markers: Number of top markers to use.
        adj_pval_threshold: Significance threshold on adjusted p-values.
        key_added: Name of .uns field to store scores.
        inplace: Return a dataframe or store it inplace.

    Returns:
        pandas.DataFrame if inplace=False, else returns AnnData (or None if modifying input inplace).
    """
    return sc.tl.marker_gene_overlap(
        adata,
        reference_markers,
        key=key,
        method=method,
        normalize=normalize,
        top_n_markers=top_n_markers,
        adj_pval_threshold=adj_pval_threshold,
        key_added=key_added,
        inplace=inplace
    )


def score_genes(
        adata: AnnData,
        gene_list: Sequence[str],
        *,
        ctrl_as_ref: bool = True,
        ctrl_size: int = 50,
        gene_pool: Sequence[str] | None = None,
        n_bins: int = 25,
        score_name: str = 'score',
        random_state: int | Any = 0,
        copy: bool = False,
        use_raw: bool | None = None,
        layer: str | None = None,
) -> AnnData | None:
    """
    Score a set of genes.

    Wraps scanpy.tl.score_genes.

    Parameters:
        adata: Annotated data matrix.
        gene_list: The list of gene names used for score calculation.
        ctrl_as_ref: Allow the algorithm to use the control genes as reference.
        ctrl_size: Number of reference genes to be sampled from each bin.
        gene_pool: Genes for sampling the reference set.
        n_bins: Number of expression level bins.
        score_name: Name of the field to be added in .obs.
        random_state: Random seed.
        copy: Copy adata or modify it inplace.
        use_raw: Whether to use raw attribute.
        layer: Key from adata.layers.

    Returns:
        Returns None if copy=False, else returns an AnnData object.
    """
    return sc.tl.score_genes(
        adata,
        gene_list,
        ctrl_as_ref=ctrl_as_ref,
        ctrl_size=ctrl_size,
        gene_pool=gene_pool,
        n_bins=n_bins,
        score_name=score_name,
        random_state=random_state,
        copy=copy,
        use_raw=use_raw,
        layer=layer
    )


def score_genes_cell_cycle(
        adata: AnnData,
        *,
        s_genes: Sequence[str],
        g2m_genes: Sequence[str],
        copy: bool = False,
        **kwargs,
) -> AnnData | None:
    """
    Score cell cycle genes.

    Wraps scanpy.tl.score_genes_cell_cycle.

    Parameters:
        adata: Annotated data matrix.
        s_genes: List of genes associated with S phase.
        g2m_genes: List of genes associated with G2M phase.
        copy: Copy adata or modify it inplace.
        **kwargs: Passed to score_genes().

    Returns:
        Returns None if copy=False, else returns an AnnData object.
        Sets 'S_score', 'G2M_score', and 'phase' in adata.obs.
    """
    return sc.tl.score_genes_cell_cycle(
        adata,
        s_genes=s_genes,
        g2m_genes=g2m_genes,
        copy=copy,
        **kwargs
    )


def rank_genes_groups_df(
        adata: AnnData,
        group: Union[str, Iterable[str], None] = None,
        *,
        key: str = "rank_genes_groups",
        pval_cutoff: Optional[float] = None,
        log2fc_min: Optional[float] = None,
        log2fc_max: Optional[float] = None,
        gene_symbols: Optional[str] = None,
) -> pd.DataFrame:
    """
    Extract differential expression results into a Pandas DataFrame.

    This function retrieves the results stored in `adata.uns[key]` (usually generated by
    `sc.tl.rank_genes_groups`) and formats them into a clean, easy-to-export table.

    Parameters:
        adata: The annotated data matrix.
        group: Which group(s) to extract results for. If None, returns all groups.
        key: The key in `adata.uns` where DEG results are stored.
        pval_cutoff: Filter results to include only genes with adjusted p-value below this threshold.
        log2fc_min: Filter results to include only genes with log2 fold change above this threshold.
        log2fc_max: Filter results to include only genes with log2 fold change below this threshold.
        gene_symbols: Key in `adata.var` to use for gene symbols. If None, uses index.

    Returns:
        A pandas DataFrame containing the differential expression statistics
        (names, scores, logfoldchanges, pvals, pvals_adj).
    """
    return sc.get.rank_genes_groups_df(
        adata,
        group=group,
        key=key,
        pval_cutoff=pval_cutoff,
        log2fc_min=log2fc_min,
        log2fc_max=log2fc_max,
        gene_symbols=gene_symbols
    )
