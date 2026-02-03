from typing import Literal, Sequence, Mapping, Union, Iterable, Any
import scanpy as sc
from anndata import AnnData
from matplotlib.axes import Axes
from matplotlib.figure import Figure


def scrublet_score_distribution(
        adata: AnnData,
        *,
        scale_hist_obs: Literal['linear', 'log', 'symlog', 'logit'] | str = 'log',
        scale_hist_sim: Literal['linear', 'log', 'symlog', 'logit'] | str = 'linear',
        figsize: tuple[float | int, float | int] = (8, 3),
        return_fig: bool = False,
        show: bool | None = None,
        save: str | bool | None = None
) -> Figure | Sequence[tuple[Axes, Axes]] | tuple[Axes, Axes] | None:
    """
    Plot histogram of doublet scores for observed transcriptomes and simulated doublets.

    Wraps scanpy.pl.scrublet_score_distribution.

    Parameters:
        adata: AnnData object resulting from scrublet().
        scale_hist_obs: Set y axis scale transformation for observed transcriptomes.
        scale_hist_sim: Set y axis scale transformation for simulated doublets.
        figsize: Figure size (width, height).
        return_fig: Return the figure object.
        show: Show the plot, do not return axis.
        save: If True or a str, save the figure.

    Returns:
        If return_fig is True, a Figure. If show==False a list of Axes.
    """
    return sc.pl.scrublet_score_distribution(
        adata, scale_hist_obs=scale_hist_obs, scale_hist_sim=scale_hist_sim,
        figsize=figsize, return_fig=return_fig, show=show, save=save
    )


def highly_variable_genes(
        adata_or_result: AnnData,
        *,
        log: bool = False,
        show: bool | None = None,
        save: bool | str | None = None,
        highly_variable_genes: bool = True
) -> None:
    """
    Plot dispersions or normalized variance versus means for genes.

    Wraps scanpy.pl.highly_variable_genes.

    Parameters:
        adata_or_result: Result of highly_variable_genes().
        log: Plot on logarithmic axes.
        show: Show the plot, do not return axis.
        save: If True or a str, save the figure.

    Returns:
        None.
    """
    sc.pl.highly_variable_genes(
        adata_or_result, log=log, show=show, save=save,
        highly_variable_genes=highly_variable_genes
    )


def highest_expr_genes(
        adata: AnnData,
        n_top: int = 30,
        *,
        layer: str | None = None,
        gene_symbols: str | None = None,
        log: bool = False,
        show: bool | None = None,
        save: str | bool | None = None,
        ax: Axes | None = None,
        **kwds
) -> Axes | None:
    """
    Fraction of counts assigned to each gene over all cells.

    Wraps scanpy.pl.highest_expr_genes.

    Parameters:
        adata: Annotated data matrix.
        n_top: Number of top genes to show.
        layer: Layer from which to pull data.
        gene_symbols: Key for field in .var that stores gene symbols.
        log: Plot x-axis in log scale.
        show: Show the plot, do not return axis.
        save: If True or a str, save the figure.
        ax: A matplotlib axes object.
        **kwds: Passed to boxplot().

    Returns:
        If show==False a Axes object.
    """
    return sc.pl.highest_expr_genes(
        adata, n_top=n_top, layer=layer, gene_symbols=gene_symbols,
        log=log, show=show, save=save, ax=ax, **kwds
    )
