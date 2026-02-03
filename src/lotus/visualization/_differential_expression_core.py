from typing import Literal, Sequence, Mapping, Union, Iterable, Any
import scanpy as sc
from anndata import AnnData
from matplotlib.axes import Axes


def rank_genes_groups(
        adata: AnnData,
        groups: str | Sequence[str] | None = None,
        *,
        n_genes: int = 20,
        gene_symbols: str | None = None,
        key: str = 'rank_genes_groups',
        fontsize: int = 8,
        ncols: int = 4,
        sharey: bool = True,
        show: bool | None = None,
        save: bool | None = None,
        ax: Axes | None = None,
        **kwds
) -> list[Axes] | None:
    """
    Plot ranking of genes.

    Wraps scanpy.pl.rank_genes_groups.

    Parameters:
        adata: Annotated data matrix.
        groups: The groups for which to show the gene ranking.
        n_genes: Number of genes to show.
        gene_symbols: Key for field in .var that stores gene symbols.
        key: Key used to store the ranking results in adata.uns.
        fontsize: Fontsize for gene names.
        ncols: Number of panels shown per row.
        sharey: Controls if the y-axis of each panels should be shared.
        show: Show the plot.
        save: Save the figure.
        ax: A matplotlib axes object.

    Returns:
        List of each group’s matplotlib axis or None if show=True.
    """
    return sc.pl.rank_genes_groups(
        adata, groups=groups, n_genes=n_genes, gene_symbols=gene_symbols,
        key=key, fontsize=fontsize, ncols=ncols, sharey=sharey, show=show,
        save=save, ax=ax, **kwds
    )


def rank_genes_groups_violin(
        adata: AnnData,
        groups: Sequence[str] | None = None,
        *,
        n_genes: int = 20,
        gene_names: Iterable[str] | None = None,
        gene_symbols: str | None = None,
        use_raw: bool | None = None,
        key: str | None = None,
        split: bool = True,
        density_norm: Literal['area', 'count', 'width'] = 'width',
        strip: bool = True,
        jitter: float | bool = True,
        size: int = 1,
        ax: Axes | None = None,
        show: bool | None = None,
        save: bool | None = None,
        scale: str = 'width'
) -> None:
    """
    Plot ranking of genes for all tested comparisons using violin plots.

    Wraps scanpy.pl.rank_genes_groups_violin.

    Parameters:
        adata: Annotated data matrix.
        groups: List of group names.
        n_genes: Number of genes to show.
        gene_names: List of genes to plot.
        use_raw: Use raw attribute of adata.
        split: Whether to split the violins or not.
        strip: Show a strip plot on top of the violin plot.
        jitter: Add jitter to the stripplot.
        show: Show the plot.
        save: Save the figure.

    Returns:
        None.
    """
    sc.pl.rank_genes_groups_violin(
        adata, groups=groups, n_genes=n_genes, gene_names=gene_names,
        gene_symbols=gene_symbols, use_raw=use_raw, key=key, split=split,
        density_norm=density_norm, strip=strip, jitter=jitter, size=size,
        ax=ax, show=show, save=save, scale=scale
    )


def rank_genes_groups_stacked_violin(
        adata: AnnData,
        groups: str | Sequence[str] | None = None,
        *,
        n_genes: int | None = None,
        groupby: str | None = None,
        gene_symbols: str | None = None,
        var_names: Sequence[str] | None = None,
        min_logfoldchange: float | None = None,
        key: str | None = None,
        show: bool | None = None,
        save: bool | None = None,
        return_fig: bool = False,
        ax: Axes | None = None,
        **kwds
) -> Any | dict | None:
    """
    Plot ranking of genes using stacked_violin plot.

    Wraps scanpy.pl.rank_genes_groups_stacked_violin.

    Parameters:
        adata: Annotated data matrix.
        groups: The groups for which to show the gene ranking.
        n_genes: Number of genes to show.
        groupby: The key of the observation grouping to consider.
        min_logfoldchange: Value to filter genes.
        return_fig: Returns StackedViolin object.
        show: Show the plot.
        save: Save the figure.

    Returns:
        If return_fig is True, returns a StackedViolin object, else if show is false, return axes dict.
    """
    return sc.pl.rank_genes_groups_stacked_violin(
        adata, groups=groups, n_genes=n_genes, groupby=groupby,
        gene_symbols=gene_symbols, var_names=var_names,
        min_logfoldchange=min_logfoldchange, key=key, show=show, save=save,
        return_fig=return_fig, ax=ax, **kwds
    )


def rank_genes_groups_heatmap(
        adata: AnnData,
        groups: str | Sequence[str] | None = None,
        *,
        n_genes: int | None = None,
        groupby: str | None = None,
        gene_symbols: str | None = None,
        var_names: Sequence[str] | None = None,
        min_logfoldchange: float | None = None,
        key: str | None = None,
        show: bool | None = None,
        save: bool | None = None,
        ax: Axes | None = None,
        **kwds
) -> None:
    """
    Plot ranking of genes using heatmap plot.

    Wraps scanpy.pl.rank_genes_groups_heatmap.

    Parameters:
        adata: Annotated data matrix.
        groups: The groups for which to show the gene ranking.
        n_genes: Number of genes to show.
        groupby: The key of the observation grouping to consider.
        min_logfoldchange: Value to filter genes.
        show: Show the plot.
        save: Save the figure.

    Returns:
        None.
    """
    sc.pl.rank_genes_groups_heatmap(
        adata, groups=groups, n_genes=n_genes, groupby=groupby,
        gene_symbols=gene_symbols, var_names=var_names,
        min_logfoldchange=min_logfoldchange, key=key, show=show, save=save,
        ax=ax, **kwds
    )


def rank_genes_groups_dotplot(
        adata: AnnData,
        groups: str | Sequence[str] | None = None,
        *,
        n_genes: int | None = None,
        groupby: str | None = None,
        values_to_plot: Literal[
                            'scores', 'logfoldchanges', 'pvals', 'pvals_adj', 'log10_pvals', 'log10_pvals_adj'] | None = None,
        var_names: Sequence[str] | Mapping[str, Sequence[str]] | None = None,
        gene_symbols: str | None = None,
        min_logfoldchange: float | None = None,
        key: str | None = None,
        show: bool | None = None,
        save: bool | None = None,
        return_fig: bool = False,
        ax: Axes | None = None,
        **kwds
) -> Any | dict | None:
    """
    Plot ranking of genes using dotplot.


    Wraps scanpy.pl.rank_genes_groups_dotplot.

    Parameters:
        adata: Annotated data matrix.
        groups: The groups for which to show the gene ranking.
        n_genes: Number of genes to show.
        values_to_plot: Instead of mean gene value, plot computed statistics.
        return_fig: Returns DotPlot object.
        show: Show the plot.
        save: Save the figure.

    Returns:
        If return_fig is True, returns a DotPlot object, else if show is false, return axes dict.
    """
    return sc.pl.rank_genes_groups_dotplot(
        adata, groups=groups, n_genes=n_genes, groupby=groupby,
        values_to_plot=values_to_plot, var_names=var_names,
        gene_symbols=gene_symbols, min_logfoldchange=min_logfoldchange,
        key=key, show=show, save=save, return_fig=return_fig, ax=ax, **kwds
    )


def rank_genes_groups_matrixplot(
        adata: AnnData,
        groups: str | Sequence[str] | None = None,
        *,
        n_genes: int | None = None,
        groupby: str | None = None,
        values_to_plot: Literal[
                            'scores', 'logfoldchanges', 'pvals', 'pvals_adj', 'log10_pvals', 'log10_pvals_adj'] | None = None,
        var_names: Sequence[str] | Mapping[str, Sequence[str]] | None = None,
        gene_symbols: str | None = None,
        min_logfoldchange: float | None = None,
        key: str | None = None,
        show: bool | None = None,
        save: bool | None = None,
        return_fig: bool = False,
        ax: Axes | None = None,
        **kwds
) -> Any | dict | None:
    """
    Plot ranking of genes using matrixplot.

    Wraps scanpy.pl.rank_genes_groups_matrixplot.

    Parameters:
        adata: Annotated data matrix.
        groups: The groups for which to show the gene ranking.
        n_genes: Number of genes to show.
        values_to_plot: Instead of mean gene value, plot computed statistics.
        return_fig: Returns MatrixPlot object.
        show: Show the plot.
        save: Save the figure.

    Returns:
        If return_fig is True, returns a MatrixPlot object, else if show is false, return axes dict.
    """
    return sc.pl.rank_genes_groups_matrixplot(
        adata, groups=groups, n_genes=n_genes, groupby=groupby,
        values_to_plot=values_to_plot, var_names=var_names,
        gene_symbols=gene_symbols, min_logfoldchange=min_logfoldchange,
        key=key, show=show, save=save, return_fig=return_fig, ax=ax, **kwds
    )


def rank_genes_groups_tracksplot(
        adata: AnnData,
        groups: str | Sequence[str] | None = None,
        *,
        n_genes: int | None = None,
        groupby: str | None = None,
        var_names: Sequence[str] | None = None,
        gene_symbols: str | None = None,
        min_logfoldchange: float | None = None,
        key: str | None = None,
        show: bool | None = None,
        save: bool | None = None,
        ax: Axes | None = None,
        **kwds
) -> None:
    """
    Plot ranking of genes using tracksplot.

    Wraps scanpy.pl.rank_genes_groups_tracksplot.

    Parameters:
        adata: Annotated data matrix.
        groups: The groups for which to show the gene ranking.
        n_genes: Number of genes to show.
        groupby: The key of the observation grouping to consider.
        show: Show the plot.
        save: Save the figure.

    Returns:
        None.
    """
    sc.pl.rank_genes_groups_tracksplot(
        adata, groups=groups, n_genes=n_genes, groupby=groupby,
        var_names=var_names, gene_symbols=gene_symbols,
        min_logfoldchange=min_logfoldchange, key=key, show=show, save=save,
        ax=ax, **kwds
    )
