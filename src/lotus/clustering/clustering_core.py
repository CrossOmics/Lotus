from typing import Literal, Sequence, Mapping, Any, Tuple, Dict, Optional, Union, Iterable
import numpy as np
import scanpy as sc
from anndata import AnnData
from scipy.sparse import csr_matrix, csc_matrix
from lotus.lineagetracker import logged


@logged
def tsne(
        adata: AnnData,
        n_pcs: int | None = None,
        *,
        use_rep: str | None = None,
        perplexity: float = 30,
        metric: str = 'euclidean',
        early_exaggeration: float = 12,
        learning_rate: float = 1000,
        random_state: int | np.random.RandomState | None = 0,
        use_fast_tsne: bool = False,
        n_jobs: int | None = None,
        key_added: str | None = None,
        copy: bool = False,
) -> AnnData | None:
    """
    Compute a t-SNE embedding.

    Wraps scanpy.tl.tsne.
    """
    return sc.tl.tsne(
        adata,
        n_pcs=n_pcs,
        use_rep=use_rep,
        perplexity=perplexity,
        metric=metric,
        early_exaggeration=early_exaggeration,
        learning_rate=learning_rate,
        random_state=random_state,
        use_fast_tsne=use_fast_tsne,
        n_jobs=n_jobs,
        key_added=key_added,
        copy=copy,
    )
@logged
def leiden(
        adata: AnnData,
        resolution: float = 1,
        *,
        restrict_to: Tuple[str, Sequence[str]] | None = None,
        random_state: int | np.random.RandomState | None = 0,
        key_added: str = 'leiden',
        adjacency: csr_matrix | csc_matrix | None = None,
        directed: bool | None = None,
        use_weights: bool = True,
        n_iterations: int = -1,
        partition_type: Any | None = None,
        neighbors_key: str | None = None,
        obsp: str | None = None,
        copy: bool = False,
        flavor: Literal['leidenalg', 'igraph'] = 'leidenalg',
        **clustering_args,
) -> AnnData | None:
    """
    Cluster cells into subgroups using the Leiden algorithm.

    This requires having run `neighbors()` or `bbknn()` first.

    Parameters:
        adata: The annotated data matrix.
        resolution: A parameter value controlling the coarseness of the clustering.
                    Higher values lead to more clusters.
        restrict_to: Restrict the clustering to the categories within the key for sample annotation.
        random_state: Change the initialization of the optimization.
        key_added: adata.obs key under which to add the cluster labels.
        adjacency: Sparse adjacency matrix of the graph, defaults to neighbors connectivities.
        directed: Whether to treat the graph as directed or undirected.
        use_weights: If True, edge weights from the graph are used in the computation.
        n_iterations: How many iterations of the Leiden clustering algorithm to perform.
        partition_type: Type of partition to use.
        neighbors_key: Use neighbors connectivities as adjacency.
        obsp: Use .obsp[obsp] as adjacency.
        copy: Whether to copy adata or modify it inplace.
        flavor: Which package’s implementation to use ('leidenalg', 'igraph').
        **clustering_args: Any further arguments to pass to find_partition().

    Returns:
        Returns None if copy=False, else returns an AnnData object.
        Sets `adata.obs[key_added]` (cluster labels).
    """
    return sc.tl.leiden(
        adata,
        resolution=resolution,
        restrict_to=restrict_to,
        random_state=random_state,
        key_added=key_added,
        adjacency=adjacency,
        directed=directed,
        use_weights=use_weights,
        n_iterations=n_iterations,
        partition_type=partition_type,
        neighbors_key=neighbors_key,
        obsp=obsp,
        copy=copy,
        flavor=flavor,
        **clustering_args
    )

@logged
def louvain(
        adata: AnnData,
        resolution: float | None = None,
        *,
        random_state: int | np.random.RandomState | None = 0,
        restrict_to: Tuple[str, Sequence[str]] | None = None,
        key_added: str = 'louvain',
        adjacency: csr_matrix | csc_matrix | None = None,
        flavor: Literal['vtraag', 'igraph', 'rapids'] = 'vtraag',
        directed: bool = True,
        use_weights: bool = False,
        partition_type: Any | None = None,
        partition_kwargs: Mapping[str, Any] = {},
        neighbors_key: str | None = None,
        obsp: str | None = None,
        copy: bool = False,
) -> AnnData | None:
    """
    Cluster cells into subgroups using the Louvain algorithm.

    This requires having run `neighbors()` or `bbknn()` first.

    Parameters:
        adata: The annotated data matrix.
        resolution: Resolution parameter (higher means more clusters). Defaults to 1.0.
        random_state: Change the initialization of the optimization.
        restrict_to: Restrict the clustering to the categories within the key for sample annotation.
        key_added: Key under which to add the cluster labels.
        adjacency: Sparse adjacency matrix of the graph.
        flavor: Choose between packages for computing the clustering ('vtraag', 'igraph', 'rapids').
        directed: Interpret the adjacency matrix as directed graph?
        use_weights: Use weights from knn graph.
        partition_type: Type of partition to use.
        partition_kwargs: Key word arguments to pass to partitioning.
        neighbors_key: Use neighbors connectivities as adjacency.
        obsp: Use .obsp[obsp] as adjacency.
        copy: Copy adata or modify it inplace.

    Returns:
        Returns None if copy=False, else returns an AnnData object.
        Sets `adata.obs[key_added]` (cluster labels).
    """
    return sc.tl.louvain(
        adata,
        resolution=resolution,
        random_state=random_state,
        restrict_to=restrict_to,
        key_added=key_added,
        adjacency=adjacency,
        flavor=flavor,
        directed=directed,
        use_weights=use_weights,
        partition_type=partition_type,
        partition_kwargs=partition_kwargs,
        neighbors_key=neighbors_key,
        obsp=obsp,
        copy=copy
    )

@logged
def dendrogram(
        adata: AnnData,
        groupby: str,
        *,
        n_pcs: int | None = None,
        use_rep: str | None = None,
        var_names: Sequence[str] | None = None,
        use_raw: bool | None = None,
        cor_method: Literal['pearson', 'kendall', 'spearman'] = 'pearson',
        linkage_method: str = 'complete',
        optimal_ordering: bool = False,
        key_added: str | None = None,
        inplace: bool = True,
) -> Dict[str, Any] | None:
    """
    Compute a hierarchical clustering for the given `groupby` categories.

    By default, the PCA representation is used unless `.X` has less than 50 variables.
    Alternatively, a list of `var_names` (e.g. genes) can be given.
    Average values of either `var_names` or components are used to compute a correlation matrix.

    Parameters:
        adata: The annotated data matrix.
        groupby: The key of the observation grouping to consider.
        n_pcs: Use this many PCs. If n_pcs==0 use .X if use_rep is None.
        use_rep: Use the indicated representation. 'X' or any key for .obsm is valid.
                 If None, the representation is chosen automatically.
        var_names: List of var_names to use for computing the hierarchical clustering.
                   If var_names is given, then use_rep and n_pcs are ignored.
        use_raw: Only when var_names is not None. Use raw attribute of adata if present.
        cor_method: Correlation method to use. Options are 'pearson', 'kendall', and 'spearman'.
        linkage_method: Linkage method to use. See scipy.cluster.hierarchy.linkage() for more information.
        optimal_ordering: Reorders the linkage matrix so that the distance between successive leaves is minimal.
        key_added: By default, the dendrogram information is added to .uns[f'dendrogram_{groupby}'].
        inplace: If True, adds dendrogram information to adata.uns[key_added], else returns the information.

    Returns:
        Returns None if inplace=True, else returns a dict with dendrogram information.
        Sets `adata.uns[f'dendrogram_{groupby}' | key_added]` if inplace=True.
    """
    return sc.tl.dendrogram(
        adata,
        groupby,
        n_pcs=n_pcs,
        use_rep=use_rep,
        var_names=var_names,
        use_raw=use_raw,
        cor_method=cor_method,
        linkage_method=linkage_method,
        optimal_ordering=optimal_ordering,
        key_added=key_added,
        inplace=inplace
    )

@logged
def umap(
        adata: AnnData,
        *,
        min_dist: float = 0.5,
        spread: float = 1.0,
        n_components: int = 2,
        maxiter: Optional[int] = None,
        alpha: float = 1.0,
        gamma: float = 1.0,
        negative_sample_rate: int = 5,
        init_pos: Union[Literal['paga', 'spectral', 'random'], np.ndarray, str, None] = 'spectral',
        random_state: Optional[Union[int, np.random.RandomState]] = 0,
        a: Optional[float] = None,
        b: Optional[float] = None,
        method: Literal['umap', 'rapids'] = 'umap',
        key_added: Optional[str] = None,
        neighbors_key: str = 'neighbors',
        copy: bool = False,
) -> Optional[AnnData]:
    """
    Embed the neighborhood graph using UMAP (Uniform Manifold Approximation and Projection).

    Requires having run `neighbors()` first. This projects the data into a lower-dimensional
    space (usually 2D) that preserves the global and local topology.

    Parameters:
        adata: The annotated data matrix.
        min_dist: The effective minimum distance between embedded points.
                  Lower values result in clumpier embeddings.
        spread: The effective scale of embedded points.
        n_components: The number of dimensions of the embedding (default: 2).
        maxiter: The number of iterations (epochs) of the optimization.
        alpha: The initial learning rate for the embedding optimization.
        gamma: Weighting applied to negative samples in optimization.
        negative_sample_rate: Number of negative samples to use per positive sample.
        init_pos: How to initialize the low dimensional embedding ('paga', 'spectral', 'random', or key).
        random_state: Seed for reproducibility.
        a: More specific parameter controlling the embedding (auto-set if None).
        b: More specific parameter controlling the embedding (auto-set if None).
        method: Chosen implementation ('umap' or 'rapids' for GPU).
        key_added: If specified, stores embedding in `obsm[key_added]`. Default is 'X_umap'.
        neighbors_key: Key in `uns` where neighbor settings are found.
        copy: Whether to return a copy or modify inplace.

    Returns:
        Returns None if copy=False, else returns an AnnData object.
        Sets `adata.obsm['X_umap']` (coordinates) and `adata.uns['umap']` (parameters).
    """
    return sc.tl.umap(
        adata,
        min_dist=min_dist,
        spread=spread,
        n_components=n_components,
        maxiter=maxiter,
        alpha=alpha,
        gamma=gamma,
        negative_sample_rate=negative_sample_rate,
        init_pos=init_pos,
        random_state=random_state,
        a=a,
        b=b,
        method=method,
        key_added=key_added,
        neighbors_key=neighbors_key,
        copy=copy
    )
