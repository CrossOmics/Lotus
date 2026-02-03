from typing import Literal, Sequence, Mapping, Union, Iterable, Any, Callable
import scanpy as sc
from anndata import AnnData
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.colors import Colormap, Normalize


def tsne(
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
        **kwargs
) -> Figure | Axes | list[Axes] | None:
    """
    Scatter plot in tSNE basis.

    Wraps scanpy.pl.tsne.

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
        show: Show the plot, do not return axis.
        save: If True or a str, save the figure.
        ax: A matplotlib axes object.
        return_fig: Return the matplotlib figure.

    Returns:
        If show==False a Axes or a list of it.
    """
    return sc.pl.tsne(
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
        return_fig=return_fig, marker=marker, **kwargs
    )


def umap(
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
        **kwargs
) -> Figure | Axes | list[Axes] | None:
    """
    Scatter plot in UMAP basis.

    Wraps scanpy.pl.umap.

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
        show: Show the plot, do not return axis.
        save: If True or a str, save the figure.
        ax: A matplotlib axes object.
        return_fig: Return the matplotlib figure.

    Returns:
        If show==False a Axes or a list of it.
    """
    return sc.pl.umap(
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
        return_fig=return_fig, marker=marker, **kwargs
    )


def diffmap(
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
        **kwargs
) -> Figure | Axes | list[Axes] | None:
    """
    Scatter plot in Diffusion Map basis.

    Wraps scanpy.pl.diffmap.

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
        show: Show the plot, do not return axis.
        save: If True or a str, save the figure.
        ax: A matplotlib axes object.
        return_fig: Return the matplotlib figure.

    Returns:
        If show==False a Axes or a list of it.
    """
    return sc.pl.diffmap(
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
        return_fig=return_fig, marker=marker, **kwargs
    )


def draw_graph(
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
        layout: Literal['fr', 'drl', 'kk', 'grid_fr', 'lgl', 'rt', 'rt_circular', 'fa'] | None = None,
        **kwargs
) -> Figure | Axes | list[Axes] | None:
    """
    Scatter plot in graph-drawing basis.

    Wraps scanpy.pl.draw_graph.

    Parameters:
        adata: Annotated data matrix.
        color: Keys for annotations of observations/cells or variables/genes.
        layout: One of the draw_graph() layouts.
        edges: Show edges.
        groups: Restrict to a few categories in categorical observation annotation.
        components: For instance, ['1,2', '2,3'].
        projection: Projection of plot ('2d' or '3d').
        legend_loc: Location of legend.
        palette: Colors to use for plotting categorical annotation groups.
        frameon: Draw a frame around the scatter plot.
        title: Provide title for panels.
        show: Show the plot, do not return axis.
        save: If True or a str, save the figure.
        ax: A matplotlib axes object.
        return_fig: Return the matplotlib figure.

    Returns:
        If show==False a Axes or a list of it.
    """
    return sc.pl.draw_graph(
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
        return_fig=return_fig, marker=marker, layout=layout, **kwargs
    )


def spatial(
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
        na_color: str | tuple[float, ...] | None = None,
        na_in_legend: bool = True,
        size: float = 1.0,
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
        basis: str = 'spatial',
        img: Any | None = None,
        img_key: str | None = None,
        library_id: str | None = None,
        crop_coord: tuple[int, int, int, int] | None = None,
        alpha_img: float = 1.0,
        bw: bool | None = False,
        spot_size: float | None = None,
        **kwargs
) -> Figure | Axes | list[Axes] | None:
    """
    Scatter plot in spatial coordinates.

    Wraps scanpy.pl.spatial. Note: Deprecated since Scanpy 1.11.0, use squidpy.pl.spatial_scatter instead.

    Parameters:
        adata: Annotated data matrix.
        color: Keys for annotations of observations/cells or variables/genes.
        library_id: library_id for Visium data.
        img_key: Key for image data.
        img: image data to plot.
        scale_factor: Scaling factor used to map from coordinate space to pixel space.
        spot_size: Diameter of spot.
        crop_coord: Coordinates to use for cropping the image.
        alpha_img: Alpha value for image.
        bw: Plot image data in gray scale.
        basis: Name of the basis.
        show: Show the plot, do not return axis.
        save: If True or a str, save the figure.
        ax: A matplotlib axes object.
        return_fig: Return the matplotlib figure.

    Returns:
        If show==False a Axes or a list of it.
    """
    return sc.pl.spatial(
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
        return_fig=return_fig, marker=marker, basis=basis, img=img,
        img_key=img_key, library_id=library_id, crop_coord=crop_coord,
        alpha_img=alpha_img, bw=bw, spot_size=spot_size, **kwargs
    )


def embedding(
        adata: AnnData,
        basis: str,
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
        **kwargs
) -> Figure | Axes | list[Axes] | None:
    """
    Scatter plot for user specified embedding basis.

    Wraps scanpy.pl.embedding.

    Parameters:
        adata: Annotated data matrix.
        basis: Name of the obsm basis to use.
        color: Keys for annotations of observations/cells or variables/genes.
        edges: Show edges.
        groups: Restrict to a few categories in categorical observation annotation.
        components: For instance, ['1,2', '2,3'].
        projection: Projection of plot ('2d' or '3d').
        legend_loc: Location of legend.
        palette: Colors to use for plotting categorical annotation groups.
        frameon: Draw a frame around the scatter plot.
        title: Provide title for panels.
        show: Show the plot, do not return axis.
        save: If True or a str, save the figure.
        ax: A matplotlib axes object.
        return_fig: Return the matplotlib figure.

    Returns:
        If show==False a Axes or a list of it.
    """
    return sc.pl.embedding(
        adata, basis, color=color, mask_obs=mask_obs, gene_symbols=gene_symbols,
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
        return_fig=return_fig, marker=marker, **kwargs
    )


def embedding_density(
        adata: AnnData,
        basis: str = 'umap',
        *,
        key: str | None = None,
        groupby: str | None = None,
        group: str | Sequence[str] | None = 'all',
        color_map: Colormap | str = 'YlOrRd',
        bg_dotsize: int | None = 80,
        fg_dotsize: int | None = 180,
        vmax: int | None = 1,
        vmin: int | None = 0,
        vcenter: int | None = None,
        norm: Normalize | None = None,
        ncols: int = 4,
        hspace: float = 0.25,
        wspace: float | None = None,
        title: str | Sequence[str] | None = None,
        show: bool | None = None,
        save: bool | str | None = None,
        ax: Axes | None = None,
        return_fig: bool | None = None,
        **kwargs
) -> Figure | Axes | None:
    """
    Plot the density of cells in an embedding.

    Wraps scanpy.pl.embedding_density.

    Parameters:
        adata: Annotated data matrix.
        basis: The embedding over which the density was calculated.
        key: Name of the .obs covariate that contains the density estimates.
        groupby: Name of the condition used in tl.embedding_density.
        group: The category in the categorical observation annotation to be plotted.
        color_map: Matplolib color map to use for density plotting.
        bg_dotsize: Dot size for background data points.
        fg_dotsize: Dot size for foreground data points.
        show: Show the plot, do not return axis.
        save: If True or a str, save the figure.
        ax: A matplotlib axes object.
        return_fig: Return the matplotlib figure.

    Returns:
        Figure or Axes or None.
    """
    return sc.pl.embedding_density(
        adata, basis=basis, key=key, groupby=groupby, group=group,
        color_map=color_map, bg_dotsize=bg_dotsize, fg_dotsize=fg_dotsize,
        vmax=vmax, vmin=vmin, vcenter=vcenter, norm=norm, ncols=ncols,
        hspace=hspace, wspace=wspace, title=title, show=show, save=save,
        ax=ax, return_fig=return_fig, **kwargs
    )
