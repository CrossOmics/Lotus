from typing import Literal, Sequence, Mapping, Union, Iterable, Any
import scanpy as sc
from anndata import AnnData
from matplotlib.axes import Axes
from matplotlib.colors import Normalize
from pandas import DataFrame


def DotPlot(
        adata: AnnData,
        var_names: str | Sequence[str] | Mapping[str, str | Sequence[str]],
        groupby: str | Sequence[str],
        *,
        use_raw: bool | None = None,
        log: bool = False,
        num_categories: int = 7,
        categories_order: Sequence[str] | None = None,
        title: str | None = None,
        figsize: tuple[float, float] | None = None,
        gene_symbols: str | None = None,
        var_group_positions: Sequence[tuple[int, int]] | None = None,
        var_group_labels: Sequence[str] | None = None,
        var_group_rotation: float | None = None,
        layer: str | None = None,
        expression_cutoff: float = 0.0,
        mean_only_expressed: bool = False,
        standard_scale: Literal['var', 'group'] | None = None,
        dot_color_df: DataFrame | None = None,
        dot_size_df: DataFrame | None = None,
        ax: Axes | None = None,
        vmin: float | None = None,
        vmax: float | None = None,
        vcenter: float | None = None,
        norm: Normalize | None = None,
        **kwds
) -> Any:
    """
    Initialize a DotPlot object for fine-tuning.

    Wraps scanpy.pl.DotPlot class constructor.

    Parameters:
        adata: Annotated data matrix.
        var_names: Genes to plot.
        groupby: Key of the observation grouping.
        use_raw: Use raw attribute.
        expression_cutoff: Threshold for binarizing expression.
        standard_scale: Standardize dimension between 0 and 1.
        ax: A matplotlib axes object.

    Returns:
        scanpy.plotting.DotPlot object. Call .show() or .make_figure() on it.
    """
    return sc.pl.DotPlot(
        adata, var_names, groupby, use_raw=use_raw, log=log,
        num_categories=num_categories, categories_order=categories_order,
        title=title, figsize=figsize, gene_symbols=gene_symbols,
        var_group_positions=var_group_positions, var_group_labels=var_group_labels,
        var_group_rotation=var_group_rotation, layer=layer,
        expression_cutoff=expression_cutoff, mean_only_expressed=mean_only_expressed,
        standard_scale=standard_scale, dot_color_df=dot_color_df,
        dot_size_df=dot_size_df, ax=ax, vmin=vmin, vmax=vmax,
        vcenter=vcenter, norm=norm, **kwds
    )


def MatrixPlot(
        adata: AnnData,
        var_names: str | Sequence[str] | Mapping[str, str | Sequence[str]],
        groupby: str | Sequence[str],
        *,
        use_raw: bool | None = None,
        log: bool = False,
        num_categories: int = 7,
        categories_order: Sequence[str] | None = None,
        title: str | None = None,
        figsize: tuple[float, float] | None = None,
        gene_symbols: str | None = None,
        var_group_positions: Sequence[tuple[int, int]] | None = None,
        var_group_labels: Sequence[str] | None = None,
        var_group_rotation: float | None = None,
        layer: str | None = None,
        standard_scale: Literal['var', 'group'] | None = None,
        ax: Axes | None = None,
        values_df: DataFrame | None = None,
        vmin: float | None = None,
        vmax: float | None = None,
        vcenter: float | None = None,
        norm: Normalize | None = None,
        **kwds
) -> Any:
    """
    Initialize a MatrixPlot object for fine-tuning.

    Wraps scanpy.pl.MatrixPlot class constructor.

    Parameters:
        adata: Annotated data matrix.
        var_names: Genes to plot.
        groupby: Key of the observation grouping.
        use_raw: Use raw attribute.
        standard_scale: Standardize dimension between 0 and 1.
        values_df: Optional dataframe with values to plot.
        ax: A matplotlib axes object.

    Returns:
        scanpy.plotting.MatrixPlot object. Call .show() or .make_figure() on it.
    """
    return sc.pl.MatrixPlot(
        adata, var_names, groupby, use_raw=use_raw, log=log,
        num_categories=num_categories, categories_order=categories_order,
        title=title, figsize=figsize, gene_symbols=gene_symbols,
        var_group_positions=var_group_positions, var_group_labels=var_group_labels,
        var_group_rotation=var_group_rotation, layer=layer,
        standard_scale=standard_scale, ax=ax, values_df=values_df,
        vmin=vmin, vmax=vmax, vcenter=vcenter, norm=norm, **kwds
    )


def StackedViolin(
        adata: AnnData,
        var_names: str | Sequence[str] | Mapping[str, str | Sequence[str]],
        groupby: str | Sequence[str],
        *,
        use_raw: bool | None = None,
        log: bool = False,
        num_categories: int = 7,
        categories_order: Sequence[str] | None = None,
        title: str | None = None,
        figsize: tuple[float, float] | None = None,
        gene_symbols: str | None = None,
        var_group_positions: Sequence[tuple[int, int]] | None = None,
        var_group_labels: Sequence[str] | None = None,
        var_group_rotation: float | None = None,
        layer: str | None = None,
        standard_scale: Literal['var', 'group'] | None = None,
        ax: Axes | None = None,
        vmin: float | None = None,
        vmax: float | None = None,
        vcenter: float | None = None,
        norm: Normalize | None = None,
        stripplot: bool = False,
        jitter: float | bool = False,
        size: float = 1,
        order: Sequence[str] | None = None,
        density_norm: Literal['area', 'count', 'width'] = 'width',
        row_palette: str | None = None,
        swap_axes: bool = False,
        **kwds
) -> Any:
    """
    Initialize a StackedViolin object for fine-tuning.

    Wraps scanpy.pl.StackedViolin class constructor.

    Parameters:
        adata: Annotated data matrix.
        var_names: Genes to plot.
        groupby: Key of the observation grouping.
        stripplot: Add a stripplot on top of the violin plot.
        jitter: Add jitter to the stripplot.
        density_norm: Method to scale the width of each violin.
        row_palette: Palette for stacked violins.
        swap_axes: Swap x and y axes.

    Returns:
        scanpy.plotting.StackedViolin object. Call .show() or .make_figure() on it.
    """
    return sc.pl.StackedViolin(
        adata, var_names, groupby, use_raw=use_raw, log=log,
        num_categories=num_categories, categories_order=categories_order,
        title=title, figsize=figsize, gene_symbols=gene_symbols,
        var_group_positions=var_group_positions, var_group_labels=var_group_labels,
        var_group_rotation=var_group_rotation, layer=layer,
        standard_scale=standard_scale, ax=ax, vmin=vmin, vmax=vmax,
        vcenter=vcenter, norm=norm, stripplot=stripplot, jitter=jitter,
        size=size, order=order, density_norm=density_norm,
        row_palette=row_palette, swap_axes=swap_axes, **kwds
    )
