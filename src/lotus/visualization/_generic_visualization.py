from typing import Literal, Sequence, Mapping, Union, Iterable, Any
import scanpy as sc
from anndata import AnnData
from matplotlib.axes import Axes
from matplotlib.colors import Colormap, Normalize


def scatter(
        adata: AnnData,
        x: str | None = None,
        y: str | None = None,
        *,
        color: str | Sequence[str] | None = None,
        use_raw: bool | None = None,
        layers: str | Sequence[str] | None = None,
        sort_order: bool = True,
        alpha: float | None = None,
        basis: Literal['pca', 'tsne', 'umap', 'diffmap', 'draw_graph_fr'] | None = None,
        groups: str | Iterable[str] | None = None,
        components: str | Sequence[str] | None = None,
        projection: Literal['2d', '3d'] = '2d',
        legend_loc: str = 'right margin',
        legend_fontsize: float | str | None = None,
        legend_fontweight: int | str | None = None,
        legend_fontoutline: float | None = None,
        color_map: str | Colormap | None = None,
        palette: Any | None = None,
        frameon: bool | None = None,
        size: float | Sequence[float] | None = None,
        title: str | Sequence[str] | None = None,
        show: bool | None = None,
        save: str | bool | None = None,
        ax: Axes | None = None,
        **kwargs
) -> Axes | list[Axes] | None:
    """
    Scatter plot along observations or variables axes.

    Wraps scanpy.pl.scatter.

    Parameters:
        adata: Annotated data matrix.
        x: x coordinate.
        y: y coordinate.
        color: Keys for annotations of observations/cells or variables/genes.
        use_raw: Whether to use raw attribute of adata.
        layers: Specify the layer for x, y and color.
        basis: String that denotes a plotting tool that computed coordinates.
        sort_order: Plot data points with higher values on top of others.
        groups: Restrict to a few categories in categorical observation annotation.
        components: For instance, ['1,2', '2,3'].
        projection: Projection of plot ('2d' or '3d').
        legend_loc: Location of legend.
        size: Point size.
        title: Provide title for panels.
        show: Show the plot, do not return axis.
        save: If True or a str, save the figure.
        ax: A matplotlib axes object.

    Returns:
        If show==False a Axes or a list of it.
    """
    return sc.pl.scatter(
        adata, x=x, y=y, color=color, use_raw=use_raw, layers=layers,
        sort_order=sort_order, alpha=alpha, basis=basis, groups=groups,
        components=components, projection=projection, legend_loc=legend_loc,
        legend_fontsize=legend_fontsize, legend_fontweight=legend_fontweight,
        legend_fontoutline=legend_fontoutline, color_map=color_map, palette=palette,
        frameon=frameon, size=size, title=title, show=show, save=save, ax=ax, **kwargs
    )


def heatmap(
        adata: AnnData,
        var_names: str | Sequence[str] | Mapping[str, str | Sequence[str]],
        groupby: str | Sequence[str],
        *,
        use_raw: bool | None = None,
        log: bool = False,
        num_categories: int = 7,
        dendrogram: bool | str = False,
        gene_symbols: str | None = None,
        var_group_positions: Sequence[tuple[int, int]] | None = None,
        var_group_labels: Sequence[str] | None = None,
        var_group_rotation: float | None = None,
        layer: str | None = None,
        standard_scale: Literal['var', 'obs'] | None = None,
        swap_axes: bool = False,
        show_gene_labels: bool | None = None,
        show: bool | None = None,
        save: str | bool | None = None,
        figsize: tuple[float, float] | None = None,
        vmin: float | None = None,
        vmax: float | None = None,
        vcenter: float | None = None,
        norm: Normalize | None = None,
        **kwds
) -> dict[str, Axes] | None:
    """
    Heatmap of the expression values of genes.

    Wraps scanpy.pl.heatmap.

    Parameters:
        adata: Annotated data matrix.
        var_names: Subset of adata.var_names.
        groupby: The key of the observation grouping to consider.
        use_raw: Use raw attribute of adata if present.
        log: Plot on logarithmic axis.
        dendrogram: If True, adds a dendrogram based on hierarchical clustering.
        standard_scale: Whether to standardize that dimension between 0 and 1.
        swap_axes: If True, x are groupby categories and y are var_names.
        show: Show the plot.
        save: Save the figure.

    Returns:
        Dict of Axes if show is False.
    """
    return sc.pl.heatmap(
        adata, var_names, groupby, use_raw=use_raw, log=log,
        num_categories=num_categories, dendrogram=dendrogram,
        gene_symbols=gene_symbols, var_group_positions=var_group_positions,
        var_group_labels=var_group_labels, var_group_rotation=var_group_rotation,
        layer=layer, standard_scale=standard_scale, swap_axes=swap_axes,
        show_gene_labels=show_gene_labels, show=show, save=save, figsize=figsize,
        vmin=vmin, vmax=vmax, vcenter=vcenter, norm=norm, **kwds
    )


def dotplot(
        adata: AnnData,
        var_names: str | Sequence[str] | Mapping[str, str | Sequence[str]],
        groupby: str | Sequence[str],
        *,
        use_raw: bool | None = None,
        log: bool = False,
        num_categories: int = 7,
        categories_order: Sequence[str] | None = None,
        expression_cutoff: float = 0.0,
        mean_only_expressed: bool = False,
        standard_scale: Literal['var', 'group'] | None = None,
        title: str | None = None,
        colorbar_title: str = 'Mean expression\nin group',
        size_title: str = 'Fraction of cells\nin group (%)',
        figsize: tuple[float, float] | None = None,
        dendrogram: bool | str = False,
        gene_symbols: str | None = None,
        var_group_positions: Sequence[tuple[int, int]] | None = None,
        var_group_labels: Sequence[str] | None = None,
        var_group_rotation: float | None = None,
        layer: str | None = None,
        swap_axes: bool = False,
        show: bool | None = None,
        save: str | bool | None = None,
        ax: Axes | None = None,
        return_fig: bool = False,
        vmin: float | None = None,
        vmax: float | None = None,
        vcenter: float | None = None,
        norm: Normalize | None = None,
        cmap: str | Colormap = 'Reds',
        dot_max: float | None = None,
        dot_min: float | None = None,
        smallest_dot: float = 0.0,
        **kwds
) -> Any | dict | None:
    """
    Make a dot plot of the expression values of var_names.

    Wraps scanpy.pl.dotplot.

    Parameters:
        adata: Annotated data matrix.
        var_names: Subset of adata.var_names.
        groupby: Key of the observation grouping.
        use_raw: Use raw attribute.
        expression_cutoff: Cutoff for binarizing gene expression.
        mean_only_expressed: If True, average only over expressing cells.
        standard_scale: Standardize dimension between 0 and 1.
        return_fig: Returns DotPlot object.
        show: Show the plot.
        save: Save the figure.

    Returns:
        If return_fig is True, returns a DotPlot object, else if show is false, return axes dict.
    """
    return sc.pl.dotplot(
        adata, var_names, groupby, use_raw=use_raw, log=log,
        num_categories=num_categories, categories_order=categories_order,
        expression_cutoff=expression_cutoff, mean_only_expressed=mean_only_expressed,
        standard_scale=standard_scale, title=title, colorbar_title=colorbar_title,
        size_title=size_title, figsize=figsize, dendrogram=dendrogram,
        gene_symbols=gene_symbols, var_group_positions=var_group_positions,
        var_group_labels=var_group_labels, var_group_rotation=var_group_rotation,
        layer=layer, swap_axes=swap_axes, show=show, save=save, ax=ax,
        return_fig=return_fig, vmin=vmin, vmax=vmax, vcenter=vcenter, norm=norm,
        cmap=cmap, dot_max=dot_max, dot_min=dot_min, smallest_dot=smallest_dot, **kwds
    )


def tracksplot(
        adata: AnnData,
        var_names: str | Sequence[str] | Mapping[str, str | Sequence[str]],
        groupby: str,
        *,
        use_raw: bool | None = None,
        log: bool = False,
        dendrogram: bool | str = False,
        gene_symbols: str | None = None,
        var_group_positions: Sequence[tuple[int, int]] | None = None,
        var_group_labels: Sequence[str] | None = None,
        layer: str | None = None,
        show: bool | None = None,
        save: str | bool | None = None,
        figsize: tuple[float, float] | None = None,
        **kwds
) -> dict[str, Axes] | None:
    """
    Compact plot of expression of a list of genes.

    Wraps scanpy.pl.tracksplot.

    Parameters:
        adata: Annotated data matrix.
        var_names: Subset of adata.var_names.
        groupby: Key of the observation grouping.
        use_raw: Use raw attribute.
        log: Plot on logarithmic axis.
        dendrogram: Add dendrogram based on hierarchical clustering.
        show: Show the plot.
        save: Save the figure.

    Returns:
        A list of Axes.
    """
    return sc.pl.tracksplot(
        adata, var_names, groupby, use_raw=use_raw, log=log,
        dendrogram=dendrogram, gene_symbols=gene_symbols,
        var_group_positions=var_group_positions, var_group_labels=var_group_labels,
        layer=layer, show=show, save=save, figsize=figsize, **kwds
    )


def violin(
        adata: AnnData,
        keys: str | Sequence[str],
        groupby: str | None = None,
        *,
        log: bool = False,
        use_raw: bool | None = None,
        stripplot: bool = True,
        jitter: float | bool = True,
        size: int = 1,
        layer: str | None = None,
        density_norm: Literal['area', 'count', 'width'] = 'width',
        order: Sequence[str] | None = None,
        multi_panel: bool | None = None,
        xlabel: str = '',
        ylabel: str | Sequence[str] | None = None,
        rotation: float | None = None,
        show: bool | None = None,
        save: bool | str | None = None,
        ax: Axes | None = None,
        **kwds
) -> Axes | Any | None:
    """
    Violin plot. Wraps seaborn.violinplot() for AnnData.

    Wraps scanpy.pl.violin.

    Parameters:
        adata: Annotated data matrix.
        keys: Keys for accessing variables or fields of .obs.
        groupby: Key of the observation grouping.
        stripplot: Add a stripplot on top of the violin plot.
        jitter: Add jitter to the stripplot.
        layer: Layer to plot.
        show: Show the plot.
        save: Save the figure.

    Returns:
        A Axes object if ax is None else None.
    """
    return sc.pl.violin(
        adata, keys, groupby=groupby, log=log, use_raw=use_raw,
        stripplot=stripplot, jitter=jitter, size=size, layer=layer,
        density_norm=density_norm, order=order, multi_panel=multi_panel,
        xlabel=xlabel, ylabel=ylabel, rotation=rotation, show=show,
        save=save, ax=ax, **kwds
    )


def stacked_violin(
        adata: AnnData,
        var_names: str | Sequence[str] | Mapping[str, str | Sequence[str]],
        groupby: str | Sequence[str],
        *,
        log: bool = False,
        use_raw: bool | None = None,
        num_categories: int = 7,
        title: str | None = None,
        colorbar_title: str = 'Median expression\nin group',
        figsize: tuple[float, float] | None = None,
        dendrogram: bool | str = False,
        gene_symbols: str | None = None,
        var_group_positions: Sequence[tuple[int, int]] | None = None,
        var_group_labels: Sequence[str] | None = None,
        standard_scale: Literal['var', 'group'] | None = None,
        var_group_rotation: float | None = None,
        layer: str | None = None,
        categories_order: Sequence[str] | None = None,
        swap_axes: bool = False,
        show: bool | None = None,
        save: bool | str | None = None,
        return_fig: bool = False,
        ax: Axes | None = None,
        vmin: float | None = None,
        vmax: float | None = None,
        vcenter: float | None = None,
        norm: Normalize | None = None,
        cmap: str | Colormap = 'Blues',
        stripplot: bool = False,
        jitter: float | bool = False,
        size: float = 1,
        row_palette: str | None = None,
        yticklabels: bool = False,
        **kwds
) -> Any | dict | None:
    """
    Stacked violin plots.

    Wraps scanpy.pl.stacked_violin.

    Parameters:
        adata: Annotated data matrix.
        var_names: Subset of adata.var_names.
        groupby: Key of the observation grouping.
        stacked_violin: Makes a compact image composed of individual violin plots.
        return_fig: Returns StackedViolin object.
        show: Show the plot.
        save: Save the figure.

    Returns:
        If return_fig is True, returns a StackedViolin object, else if show is false, return axes dict.
    """
    return sc.pl.stacked_violin(
        adata, var_names, groupby, log=log, use_raw=use_raw,
        num_categories=num_categories, title=title, colorbar_title=colorbar_title,
        figsize=figsize, dendrogram=dendrogram, gene_symbols=gene_symbols,
        var_group_positions=var_group_positions, var_group_labels=var_group_labels,
        standard_scale=standard_scale, var_group_rotation=var_group_rotation,
        layer=layer, categories_order=categories_order, swap_axes=swap_axes,
        show=show, save=save, return_fig=return_fig, ax=ax, vmin=vmin, vmax=vmax,
        vcenter=vcenter, norm=norm, cmap=cmap, stripplot=stripplot, jitter=jitter,
        size=size, row_palette=row_palette, yticklabels=yticklabels, **kwds
    )


def matrixplot(
        adata: AnnData,
        var_names: str | Sequence[str] | Mapping[str, str | Sequence[str]],
        groupby: str | Sequence[str],
        *,
        use_raw: bool | None = None,
        log: bool = False,
        num_categories: int = 7,
        categories_order: Sequence[str] | None = None,
        figsize: tuple[float, float] | None = None,
        dendrogram: bool | str = False,
        title: str | None = None,
        cmap: str | Colormap = 'viridis',
        colorbar_title: str = 'Mean expression\nin group',
        gene_symbols: str | None = None,
        var_group_positions: Sequence[tuple[int, int]] | None = None,
        var_group_labels: Sequence[str] | None = None,
        var_group_rotation: float | None = None,
        layer: str | None = None,
        standard_scale: Literal['var', 'group'] | None = None,
        values_df: Any | None = None,
        swap_axes: bool = False,
        show: bool | None = None,
        save: str | bool | None = None,
        ax: Axes | None = None,
        return_fig: bool = False,
        vmin: float | None = None,
        vmax: float | None = None,
        vcenter: float | None = None,
        norm: Normalize | None = None,
        **kwds
) -> Any | dict | None:
    """
    Create a heatmap of the mean expression values per group.

    Wraps scanpy.pl.matrixplot.

    Parameters:
        adata: Annotated data matrix.
        var_names: Subset of adata.var_names.
        groupby: Key of the observation grouping.
        return_fig: Returns MatrixPlot object.
        show: Show the plot.
        save: Save the figure.

    Returns:
        If return_fig is True, returns a MatrixPlot object, else if show is false, return axes dict.
    """
    return sc.pl.matrixplot(
        adata, var_names, groupby, use_raw=use_raw, log=log,
        num_categories=num_categories, categories_order=categories_order,
        figsize=figsize, dendrogram=dendrogram, title=title, cmap=cmap,
        colorbar_title=colorbar_title, gene_symbols=gene_symbols,
        var_group_positions=var_group_positions, var_group_labels=var_group_labels,
        var_group_rotation=var_group_rotation, layer=layer,
        standard_scale=standard_scale, values_df=values_df, swap_axes=swap_axes,
        show=show, save=save, ax=ax, return_fig=return_fig, vmin=vmin, vmax=vmax,
        vcenter=vcenter, norm=norm, **kwds
    )


def clustermap(
        adata: AnnData,
        obs_keys: str | None = None,
        *,
        use_raw: bool | None = None,
        show: bool | None = None,
        save: bool | str | None = None,
        **kwds
) -> Any | None:
    """
    Hierarchically-clustered heatmap. Wraps seaborn.clustermap() for AnnData.

    Wraps scanpy.pl.clustermap.

    Parameters:
        adata: Annotated data matrix.
        obs_keys: Categorical annotation to plot with a different color map.
        use_raw: Whether to use raw attribute of adata.
        show: Show the plot.
        save: Save the figure.

    Returns:
        If show is False, a ClusterGrid object.
    """
    return sc.pl.clustermap(
        adata, obs_keys=obs_keys, use_raw=use_raw, show=show, save=save, **kwds
    )


def ranking(
        adata: AnnData,
        attr: Literal['var', 'obs', 'uns', 'varm', 'obsm'],
        keys: str | Sequence[str],
        *,
        dictionary: dict | None = None,
        indices: Sequence[int] | None = None,
        labels: Sequence[str] | None = None,
        color: str = 'black',
        n_points: int = 30,
        log: bool = False,
        include_lowest: bool = False,
        show: bool | None = None,
) -> Any | None:
    """
    Plot rankings.

    Wraps scanpy.pl.ranking.

    Parameters:
        adata: The data.
        attr: The attribute of AnnData that contains the score.
        keys: The scores to look up an array from the attribute of adata.
        n_points: Number of points to plot.
        show: Show the plot.

    Returns:
        Returns matplotlib gridspec with access to the axes.
    """
    return sc.pl.ranking(
        adata, attr, keys, dictionary=dictionary, indices=indices, labels=labels,
        color=color, n_points=n_points, log=log, include_lowest=include_lowest,
        show=show
    )


def dendrogram(
        adata: AnnData,
        groupby: str,
        *,
        dendrogram_key: str | None = None,
        orientation: Literal['top', 'bottom', 'left', 'right'] = 'top',
        remove_labels: bool = False,
        show: bool | None = None,
        save: str | bool | None = None,
        ax: Axes | None = None
) -> Axes | None:
    """
    Plot a dendrogram of the categories defined in groupby.

    Wraps scanpy.pl.dendrogram.

    Parameters:
        adata: Annotated data matrix.
        groupby: Categorical data column used to create the dendrogram.
        dendrogram_key: Key under which the dendrogram information was stored.
        orientation: Origin of the tree.
        remove_labels: Don't draw labels.
        show: Show the plot.
        save: Save the figure.

    Returns:
        matplotlib.axes.Axes
    """
    return sc.pl.dendrogram(
        adata, groupby, dendrogram_key=dendrogram_key, orientation=orientation,
        remove_labels=remove_labels, show=show, save=save, ax=ax
    )
