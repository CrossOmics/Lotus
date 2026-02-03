from typing import Literal, Sequence, Mapping, Union, Iterable, Any, Callable
import scanpy as sc
from anndata import AnnData
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.colors import Colormap, Normalize


def pca(
        adata: AnnData,
        *,
        color: str | Sequence[str] | None = None,
        mask_obs: Any | None = None,
        gene_symbols: str | None = None,
        use_raw: bool | None = None,
        sort_order: bool = True,
        edges: bool = False,
        edges_width: float = 0.1,
        edges_color: str | Sequence[float] | Sequence[str] = 'grey',
        neighbors_key: str | None = None,
        arrows: bool = False,
        arrows_kwds: Mapping[str, Any] | None = None,
        groups: str | Sequence[str] | None = None,
        components: str | Sequence[str] | None = None,
        dimensions: tuple[int, int] | Sequence[tuple[int, int]] | None = None,
        layer: str | None = None,
        projection: Literal['2d', '3d'] = '2d',
        scale_factor: float | None = None,
        color_map: Colormap | str | None = None,
        cmap: Colormap | str | None = None,
        palette: str | Sequence[str] | Any | None = None,
        na_color: str | tuple[float, ...] = 'lightgray',
        na_in_legend: bool = True,
        size: float | Sequence[float] | None = None,
        frameon: bool | None = None,
        legend_fontsize: float | str | None = None,
        legend_fontweight: int | str = 'bold',
        legend_loc: str = 'right margin',
        legend_fontoutline: int | None = None,
        colorbar_loc: str | None = 'right',
        vmax: str | float | Callable | Sequence | None = None,
        vmin: str | float | Callable | Sequence | None = None,
        vcenter: str | float | Callable | Sequence | None = None,
        norm: Normalize | None = None,
        add_outline: bool | None = False,
        outline_width: tuple[float, float] = (0.3, 0.05),
        outline_color: tuple[str, str] = ('black', 'white'),
        ncols: int = 4,
        hspace: float = 0.25,
        wspace: float | None = None,
        title: str | Sequence[str] | None = None,
        show: bool | None = None,
        save: bool | str | None = None,
        ax: Axes | None = None,
        return_fig: bool | None = None,
        marker: str = '.',
        annotate_var_explained: bool = False,
        **kwargs
) -> Figure | Axes | list[Axes] | None:
    """
    Scatter plot in PCA coordinates.

    Wraps scanpy.pl.pca.

    Parameters:
        adata: Annotated data matrix.
        color: Keys for annotations of observations/cells or variables/genes.
        edges: Show edges.
        groups: Restrict to a few categories in categorical observation annotation.
        components: For instance, ['1,2', '2,3'].
        projection: Projection of plot ('2d' or '3d').
        legend_loc: Location of legend.
        palette: Colors to use for plotting categorical annotation groups.
        frameon: Draw a frame around the scatter plot.
        title: Provide title for panels.
        annotate_var_explained: Annotate the explained variance.
        show: Show the plot, do not return axis.
        save: If True or a str, save the figure.
        ax: A matplotlib axes object.
        return_fig: Return the matplotlib figure.

    Returns:
        If show==False a Axes or a list of it.
    """
    return sc.pl.pca(
        adata, color=color, mask_obs=mask_obs, gene_symbols=gene_symbols,
        use_raw=use_raw, sort_order=sort_order, edges=edges,
        edges_width=edges_width, edges_color=edges_color,
        neighbors_key=neighbors_key, arrows=arrows, arrows_kwds=arrows_kwds,
        groups=groups, components=components, dimensions=dimensions,
        layer=layer, projection=projection, scale_factor=scale_factor,
        color_map=color_map, cmap=cmap, palette=palette, na_color=na_color,
        na_in_legend=na_in_legend, size=size, frameon=frameon,
        legend_fontsize=legend_fontsize, legend_fontweight=legend_fontweight,
        legend_loc=legend_loc, legend_fontoutline=legend_fontoutline,
        colorbar_loc=colorbar_loc, vmax=vmax, vmin=vmin, vcenter=vcenter,
        norm=norm, add_outline=add_outline, outline_width=outline_width,
        outline_color=outline_color, ncols=ncols, hspace=hspace,
        wspace=wspace, title=title, show=show, save=save, ax=ax,
        return_fig=return_fig, marker=marker,
        annotate_var_explained=annotate_var_explained, **kwargs
    )


def pca_loadings(
        adata: AnnData,
        components: str | Sequence[int] | None = None,
        *,
        include_lowest: bool = True,
        n_points: int | None = None,
        show: bool | None = None,
        save: str | bool | None = None
) -> None:
    """
    Rank genes according to contributions to PCs.

    Wraps scanpy.pl.pca_loadings.

    Parameters:
        adata: Annotated data matrix.
        components: For example, '1,2,3' means [1, 2, 3].
        include_lowest: Whether to show the variables with both highest and lowest loadings.
        n_points: Number of variables to plot for each component.
        show: Show the plot, do not return axis.
        save: If True or a str, save the figure.

    Returns:
        None.
    """
    sc.pl.pca_loadings(
        adata, components=components, include_lowest=include_lowest,
        n_points=n_points, show=show, save=save
    )


def pca_variance_ratio(
        adata: AnnData,
        n_pcs: int = 30,
        *,
        log: bool = False,
        show: bool | None = None,
        save: bool | str | None = None
) -> None:
    """
    Plot the variance ratio.

    Wraps scanpy.pl.pca_variance_ratio.

    Parameters:
        adata: Annotated data matrix.
        n_pcs: Number of PCs to show.
        log: Plot on logarithmic scale.
        show: Show the plot, do not return axis.
        save: If True or a str, save the figure.

    Returns:
        None.
    """
    sc.pl.pca_variance_ratio(
        adata, n_pcs=n_pcs, log=log, show=show, save=save
    )


def pca_overview(
        adata: AnnData,
        *,
        color: str | Sequence[str] | None = None,
        use_raw: bool | None = None,
        sort_order: bool = True,
        groups: str | Sequence[str] | None = None,
        dimensions: tuple[int, int] | Sequence[tuple[int, int]] | None = None,
        components: str | Sequence[str] | None = None,
        projection: Literal['2d', '3d'] = '2d',
        legend_loc: str = 'right margin',
        legend_fontsize: float | str | None = None,
        legend_fontweight: int | str = 'bold',
        legend_fontoutline: int | None = None,
        colorbar_loc: str | None = 'right',
        size: float | Sequence[float] | None = None,
        color_map: Colormap | str | None = None,
        palette: str | Sequence[str] | Any | None = None,
        na_color: str | tuple[float, ...] = 'lightgray',
        na_in_legend: bool = True,
        frameon: bool | None = None,
        title: str | Sequence[str] | None = None,
        vmin: str | float | Callable | Sequence | None = None,
        vmax: str | float | Callable | Sequence | None = None,
        vcenter: str | float | Callable | Sequence | None = None,
        add_outline: bool | None = False,
        outline_color: tuple[str, str] = ('black', 'white'),
        outline_width: tuple[float, float] = (0.3, 0.05),
        ncols: int | None = None,
        wspace: float | None = None,
        hspace: float | None = None,
        return_fig: bool | None = None,
        show: bool | None = None,
        save: bool | str | None = None,
        **kwargs
) -> Figure | Axes | list[Axes] | None:
    """
    Plot PCA results (scatter plot overview).

    Wraps scanpy.pl.pca_overview.

    Parameters:
        adata: Annotated data matrix.
        color: Keys for observation/cell annotation.
        use_raw: Use raw attribute of adata.
        sort_order: Plot data points with higher values on top.
        groups: Restrict to a few categories.
        components: For instance, ['1,2', '2,3'].
        projection: Projection of plot.
        legend_loc: Location of legend.
        palette: Colors to use for plotting categorical annotation groups.
        frameon: Draw a frame around the scatter plot.
        title: Provide title for panels.
        show: Show the plot, do not return axis.
        save: If True or a str, save the figure.
        return_fig: Return the matplotlib figure.

    Returns:
        If show==False a Axes or a list of it.
    """
    return sc.pl.pca_overview(
        adata, color=color, use_raw=use_raw, sort_order=sort_order,
        groups=groups, dimensions=dimensions, components=components,
        projection=projection, legend_loc=legend_loc,
        legend_fontsize=legend_fontsize, legend_fontweight=legend_fontweight,
        legend_fontoutline=legend_fontoutline, colorbar_loc=colorbar_loc,
        size=size, color_map=color_map, palette=palette, na_color=na_color,
        na_in_legend=na_in_legend, frameon=frameon, title=title,
        vmin=vmin, vmax=vmax, vcenter=vcenter, add_outline=add_outline,
        outline_color=outline_color, outline_width=outline_width,
        ncols=ncols, wspace=wspace, hspace=hspace, return_fig=return_fig,
        show=show, save=save, **kwargs
    )
