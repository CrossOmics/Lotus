from typing import Literal, Sequence, Mapping, Union, Iterable, Any
import scanpy as sc
from anndata import AnnData
from matplotlib.axes import Axes
from matplotlib.colors import Colormap, Normalize
from pandas import DataFrame

def dpt_groups_pseudotime(
        adata: AnnData,
        *,
        color_map: str | Colormap | None = None,
        palette: Sequence[str] | Any | None = None,
        show: bool | None = None,
        save: bool | str | None = None,
        marker: str | Sequence[str] = '.'
) -> None:
    """
    Plot groups and pseudotime.

    Wraps scanpy.pl.dpt_groups_pseudotime.

    Parameters:
        adata: Annotated data matrix.
        color_map: Color map to use for continous variables.
        palette: Colors to use for plotting categorical annotation groups.
        show: Show the plot, do not return axis.
        save: If True or a str, save the figure.
        marker: Marker style.

    Returns:
        None.
    """
    sc.pl.dpt_groups_pseudotime(
        adata, color_map=color_map, palette=palette, show=show,
        save=save, marker=marker
    )


def dpt_timeseries(
        adata: AnnData,
        *,
        color_map: str | Colormap | None = None,
        show: bool | None = None,
        save: bool | str | None = None,
        as_heatmap: bool = True,
        marker: str = '.'
) -> None:
    """
    Heatmap of pseudotime series.

    Wraps scanpy.pl.dpt_timeseries.

    Parameters:
        adata: Annotated data matrix.
        color_map: Color map.
        as_heatmap: Plot the timeseries as heatmap.
        show: Show the plot.
        save: Save the figure.
        marker: Marker style.

    Returns:
        None.
    """
    sc.pl.dpt_timeseries(
        adata, color_map=color_map, show=show, save=save,
        as_heatmap=as_heatmap, marker=marker
    )


def paga(
        adata: AnnData,
        *,
        threshold: float | None = None,
        color: str | Mapping[str | int, Mapping[Any, float]] | None = None,
        layout: Literal['fr', 'drl', 'kk', 'grid_fr', 'lgl', 'rt', 'rt_circular', 'fa', 'eq_tree'] | None = None,
        layout_kwds: Mapping[str, Any] = dict(),
        init_pos: Any | None = None,
        root: int | str | Sequence[int] | None = 0,
        labels: str | Sequence[str] | Mapping[str, str] | None = None,
        single_component: bool = False,
        solid_edges: str = 'connectivities',
        dashed_edges: str | None = None,
        transitions: str | None = None,
        fontsize: int | None = None,
        fontweight: str = 'bold',
        fontoutline: int | None = None,
        text_kwds: Mapping[str, Any] = dict(),
        node_size_scale: float = 1.0,
        node_size_power: float = 0.5,
        edge_width_scale: float = 1.0,
        min_edge_width: float | None = None,
        max_edge_width: float | None = None,
        arrowsize: int = 30,
        title: str | None = None,
        left_margin: float = 0.01,
        random_state: int | None = 0,
        pos: Any | None = None,
        normalize_to_color: bool = False,
        cmap: str | Colormap | None = None,
        cax: Axes | None = None,
        colorbar: bool | None = None,
        cb_kwds: Mapping[str, Any] = dict(),
        frameon: bool | None = None,
        add_pos: bool = True,
        export_to_gexf: bool = False,
        use_raw: bool = True,
        colors: Any | None = None,
        groups: Any | None = None,
        plot: bool = True,
        show: bool | None = None,
        save: bool | str | None = None,
        ax: Axes | None = None
) -> Axes | list[Axes] | None:
    """
    Plot the PAGA graph through thresholding low-connectivity edges.


    Wraps scanpy.pl.paga.

    Parameters:
        adata: Annotated data matrix.
        threshold: Do not draw edges for weights below this threshold.
        color: Gene name or obs annotation defining the node colors.
        layout: Plotting layout that computes positions.
        root: Index of the root node for tree layouts.
        node_size_scale: Increase or decrease the size of the nodes.
        edge_width_scale: Edge width scale.
        show: Show the plot.
        save: Save the figure.

    Returns:
        Axes or list of Axes if show is False.
    """
    return sc.pl.paga(
        adata, threshold=threshold, color=color, layout=layout,
        layout_kwds=layout_kwds, init_pos=init_pos, root=root, labels=labels,
        single_component=single_component, solid_edges=solid_edges,
        dashed_edges=dashed_edges, transitions=transitions, fontsize=fontsize,
        fontweight=fontweight, fontoutline=fontoutline, text_kwds=text_kwds,
        node_size_scale=node_size_scale, node_size_power=node_size_power,
        edge_width_scale=edge_width_scale, min_edge_width=min_edge_width,
        max_edge_width=max_edge_width, arrowsize=arrowsize, title=title,
        left_margin=left_margin, random_state=random_state, pos=pos,
        normalize_to_color=normalize_to_color, cmap=cmap, cax=cax,
        colorbar=colorbar, cb_kwds=cb_kwds, frameon=frameon, add_pos=add_pos,
        export_to_gexf=export_to_gexf, use_raw=use_raw, colors=colors,
        groups=groups, plot=plot, show=show, save=save, ax=ax
    )


def paga_path(
        adata: AnnData,
        nodes: Sequence[str | int],
        keys: Sequence[str],
        *,
        use_raw: bool = True,
        annotations: Sequence[str] = ('dpt_pseudotime',),
        color_map: str | Colormap | None = None,
        color_maps_annotations: Mapping[str, str | Colormap] = dict(dpt_pseudotime='Greys'),
        palette_groups: Sequence[str] | None = None,
        n_avg: int = 1,
        groups_key: str | None = None,
        xlim: tuple[float | None, float | None] = (None, None),
        title: str | None = None,
        left_margin: float | None = None,
        ytick_fontsize: int | None = None,
        title_fontsize: int | None = None,
        show_node_names: bool = True,
        show_yticks: bool = True,
        show_colorbar: bool = True,
        legend_fontsize: int | None = None,
        legend_fontweight: str | None = None,
        normalize_to_zero_one: bool = False,
        as_heatmap: bool = True,
        return_data: bool = False,
        show: bool | None = None,
        save: bool | str | None = None,
        ax: Axes | None = None
) -> tuple[Axes, DataFrame] | Axes | DataFrame | None:
    """
    Gene expression and annotation changes along paths in the abstracted graph.

    Wraps scanpy.pl.paga_path.

    Parameters:
        adata: Annotated data matrix.
        nodes: A path through nodes of the abstracted graph.
        keys: Genes or annotations to plot.
        as_heatmap: Plot the timeseries as heatmap.
        return_data: Return the timeseries data in addition to the axes.
        show: Show the plot.
        save: Save the figure.

    Returns:
        Axes object, or tuple of (Axes, DataFrame) if return_data is True.
    """
    return sc.pl.paga_path(
        adata, nodes, keys, use_raw=use_raw, annotations=annotations,
        color_map=color_map, color_maps_annotations=color_maps_annotations,
        palette_groups=palette_groups, n_avg=n_avg, groups_key=groups_key,
        xlim=xlim, title=title, left_margin=left_margin,
        ytick_fontsize=ytick_fontsize, title_fontsize=title_fontsize,
        show_node_names=show_node_names, show_yticks=show_yticks,
        show_colorbar=show_colorbar, legend_fontsize=legend_fontsize,
        legend_fontweight=legend_fontweight,
        normalize_to_zero_one=normalize_to_zero_one, as_heatmap=as_heatmap,
        return_data=return_data, show=show, save=save, ax=ax
    )


def paga_compare(
        adata: AnnData,
        basis: str | None = None,
        *,
        edges: bool = False,
        color: str | None = None,
        alpha: float | None = None,
        groups: str | None = None,
        components: str | None = None,
        projection: str = '2d',
        legend_loc: str = 'on data',
        legend_fontsize: int | None = None,
        legend_fontweight: str = 'bold',
        legend_fontoutline: int | None = None,
        color_map: str | Colormap | None = None,
        palette: str | None = None,
        frameon: bool = False,
        size: float | None = None,
        title: str | None = None,
        right_margin: float | None = None,
        left_margin: float = 0.05,
        show: bool | None = None,
        save: bool | str | None = None,
        title_graph: str | None = None,
        groups_graph: str | None = None,
        pos: Any | None = None,
        **paga_graph_params
) -> list[Axes] | None:
    """
    Scatter and PAGA graph side-by-side.

    Wraps scanpy.pl.paga_compare.

    Parameters:
        adata: Annotated data matrix.
        basis: Basis to use for scatter plot.
        edges: Show edges in scatter plot.
        color: Coloring of nodes/cells.
        show: Show the plot.
        save: Save the figure.

    Returns:
        A list of Axes if show is False.
    """
    return sc.pl.paga_compare(
        adata, basis=basis, edges=edges, color=color, alpha=alpha,
        groups=groups, components=components, projection=projection,
        legend_loc=legend_loc, legend_fontsize=legend_fontsize,
        legend_fontweight=legend_fontweight,
        legend_fontoutline=legend_fontoutline, color_map=color_map,
        palette=palette, frameon=frameon, size=size, title=title,
        right_margin=right_margin, left_margin=left_margin, show=show,
        save=save, title_graph=title_graph, groups_graph=groups_graph,
        pos=pos, **paga_graph_params
    )


def correlation_matrix(
        adata: AnnData,
        groupby: str,
        *,
        show_correlation_numbers: bool = False,
        dendrogram: bool | str | None = None,
        figsize: tuple[float, float] | None = None,
        show: bool | None = None,
        save: bool | str | None = None,
        ax: Axes | None = None,
        vmin: float | None = None,
        vmax: float | None = None,
        vcenter: float | None = None,
        norm: Normalize | None = None,
        **kwds
) -> list[Axes] | None:
    """
    Plot the correlation matrix computed as part of dendrogram.


    Wraps scanpy.pl.correlation_matrix.

    Parameters:
        adata: Annotated data matrix.
        groupby: Categorical data column used to create the dendrogram.
        show_correlation_numbers: Plot the correlation on top of each cell.
        dendrogram: Add dendrogram based on hierarchical clustering.
        show: Show the plot.
        save: Save the figure.

    Returns:
        If show=False, returns a list of Axes objects.
    """
    return sc.pl.correlation_matrix(
        adata, groupby, show_correlation_numbers=show_correlation_numbers,
        dendrogram=dendrogram, figsize=figsize, show=show, save=save,
        ax=ax, vmin=vmin, vmax=vmax, vcenter=vcenter, norm=norm, **kwds
    )
