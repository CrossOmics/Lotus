from typing import Collection

import numpy as np
import pandas as pd
from anndata import AnnData
import scanpy.preprocessing as pp_module

'''
implement domain knowledge or customized methods in this core file
'''


def run_preprocessing(
        adata: AnnData,
        *,
        n_pcs: int | None = None,
        target_sum: float = 1e4,
        n_top_genes: int | None = None,
        n_neighbors: int = 15,
        use_rep: str = "X_pca",
        save_raw: bool = True,
        raw_layer: str = "raw_counts",
        min_genes: int | None = None,
        min_cells: int | None = None,
        min_counts: int | None = None,
        max_counts: int | None = None,
        max_genes: int | None = None,
        pct_mt_max: float | None = None,
        hvg_flavor: str = "seurat_v3",
        batch_key: str | None = None,
        regress_out_keys: list[str] | None = None,
        use_combat: bool = False,
        qc_vars: Collection[str] = (),
) -> None:
    """
    Run a full preprocessing pipeline on an AnnData object.

    Steps included:
        1. Restore raw counts if preprocessing was run before.
        2. Save raw counts to a layer for reproducibility.
        3. Compute QC metrics.
        4. Filter low-quality cells/genes (counts, genes, MT%).
        5. Normalize total counts and apply log1p.
        6. Optionally regress out unwanted covariates.
        7. Optionally apply ComBat batch correction.
        8. Select highly variable genes and subset.
        9. Scale features.
       10. Run PCA.
       11. Build the neighbor graph.

    The function updates `adata` in place and populates:
           - QC metrics in `adata.obs` and `adata.var`
           - normalized/log1p data in `adata.X` and `adata.raw`
           - HVGs in `adata.var['highly_variable']`
           - PCA embeddings in `adata.obsm['X_pca']`
           - neighbor graph in `adata.uns['neighbors']` and `adata.obsp`
    """
    # 1. Restore raw counts if preprocessing was run before
    # This prevents errors when re-running preprocessing on already processed data
    if adata.raw is not None or (save_raw and raw_layer and raw_layer in adata.layers):
        if save_raw and raw_layer and raw_layer in adata.layers:
            # Restore from raw_layer (original count matrix)
            print("Restoring from raw counts layer for re-preprocessing...")
            adata.X = adata.layers[raw_layer].copy()
            # Clear preprocessing results that will be recomputed
            keys_to_clear_obsm = ['X_pca', 'X_latent']
            keys_to_clear_uns = ['neighbors', 'pca']
            keys_to_clear_var = ['highly_variable', 'means', 'dispersions', 'dispersions_norm',
                                 'highly_variable_nbatches', 'highly_variable_intersection']

            for k in keys_to_clear_obsm: adata.obsm.pop(k, None)
            for k in keys_to_clear_uns: adata.uns.pop(k, None)
            for k in keys_to_clear_var: adata.var.pop(k, None)

        elif adata.raw is not None:
            # Restore from adata.raw (normalized+log1p data)
            # We try to use this, but warn that it might not be raw counts
            print("Restoring from adata.raw for re-preprocessing...")
            print("Warning: Cannot fully restore from adata.raw. Using adata.raw.X as starting point.")
            adata.X = adata.raw.X.copy()
            # Clear preprocessing results
            for k in ['X_pca', 'X_latent']: adata.obsm.pop(k, None)
            for k in ['neighbors', 'pca']: adata.uns.pop(k, None)

    # 2. Save raw counts to a layer for reproducibility
    if save_raw and raw_layer:
        if raw_layer not in adata.layers:
            adata.layers[raw_layer] = adata.X.copy()

    # 3. Compute QC metrics
    calculate_qc_metrics(adata, qc_vars=qc_vars, inplace=True)

    # 4. Filter low-quality cells/genes
    # Filter cells based on counts and genes
    if min_counts is not None:
        filter_cells(adata, min_counts=min_counts, inplace=True)

    if min_genes is not None:
        filter_cells(adata, min_genes=min_genes, inplace=True)

    if max_counts is not None:
        filter_cells(adata, max_counts=max_counts, inplace=True)

    if max_genes is not None:
        filter_cells(adata, max_genes=max_genes, inplace=True)

    # Filter genes based on cell appearance
    if min_cells is not None:
        filter_genes(adata, min_cells=min_cells, inplace=True)

    # Filter based on mitochondrial percentage if requested
    if pct_mt_max is not None:
        pct_mt_key = None
        if 'pct_counts_mt' in adata.obs.columns:
            pct_mt_key = 'pct_counts_mt'
        elif 'pct_mt' in adata.obs.columns:
            pct_mt_key = 'pct_mt'

        if pct_mt_key is not None:
            n_before = adata.n_obs
            # Standard Scanpy way to filter obs in-place
            mask = adata.obs[pct_mt_key] < pct_mt_max
            adata._inplace_subset_obs(mask)
            if n_before != adata.n_obs:
                print(f"Filtered {n_before - adata.n_obs} cells with {pct_mt_key} >= {pct_mt_max}")

    # 5. Normalize total counts and apply log1p
    normalization(adata, target_sum=target_sum, inplace=True)
    log1p(adata)

    # Save the normalized, log-transformed data to adata.raw (Standard Scanpy Workflow)
    if save_raw:
        adata.raw = adata

    # 6. Optionally regress out unwanted covariates
    if regress_out_keys is not None:
        regress_out(adata, keys=regress_out_keys)

    # 7. Optionally apply ComBat batch correction
    if use_combat and batch_key is not None:
        if batch_key in adata.obs.columns:
            combat(adata, key=batch_key, inplace=True)
        else:
            print(f"Warning: batch_key '{batch_key}' not found in adata.obs. Skipping ComBat correction.")

    # 8. Select highly variable genes and subset
    # This reduces dimensionality before scaling and PCA
    if n_top_genes is None:
        n_top_genes = min(2000, adata.n_vars)

    highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor=hvg_flavor,
        subset=True,
        inplace=True
    )

    # 9. Scale features (Zero-center and unit variance)
    scale(adata, zero_center=True, max_value=10)

    # 10. Run PCA
    pca(adata, n_comps=n_pcs)

    # 11. Build the neighbor graph
    neighbors(adata, n_neighbors=n_neighbors, use_rep=use_rep)


def calculate_qc_metrics(
        adata: AnnData,
        *,
        expr_type: str = "counts",
        var_type: str = "genes",
        qc_vars: Collection[str] | str = (),
        percent_top: Collection[int] | None = (50, 100, 200, 500),
        layer: str | None = None,
        use_raw: bool = False,
        inplace: bool = False,
        log1p: bool = True,
        parallel: bool | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    """
    Calculates a number of qc metrics for an AnnData object.

    Returns:
        Depending on `inplace`, returns calculated metrics (as DataFrame) or updates
        adata's obs and var.
    """
    return pp_module.calculate_qc_metrics(
        adata=adata,
        expr_type=expr_type,
        var_type=var_type,
        qc_vars=qc_vars,
        percent_top=percent_top,
        layer=layer,
        use_raw=use_raw,
        inplace=inplace,
        log1p=log1p,
        parallel=parallel
    )


def filter_cells(
        data: AnnData,
        *,
        min_counts: int | None = None,
        min_genes: int | None = None,
        max_counts: int | None = None,
        max_genes: int | None = None,
        inplace: bool = True,
        copy: bool = False,
) -> AnnData | tuple[np.ndarray, np.ndarray] | None:
    """
    Filter cell outliers based on counts and numbers of genes expressed.
    For instance, only keep cells with at least `min_counts` counts or `min_genes`
    genes expressed.

    Parameters:
        data: The annotated data matrix of shape n_obs × n_vars.
        min_counts: Minimum number of counts required for a cell to pass filtering.
        min_genes: Minimum number of genes expressed required for a cell to pass filtering.
        max_counts: Maximum number of counts required for a cell to pass filtering.
        max_genes: Maximum number of genes expressed required for a cell to pass filtering.
        inplace: Perform computation inplace or return result.
        copy: Whether to modify specific copy of object.

    Returns:
        Depending on `inplace`, returns the filtered AnnData object, a tuple of masks,
        or None (updates data in place).
    """
    return pp_module.filter_cells(
        data,
        min_counts=min_counts,
        min_genes=min_genes,
        max_counts=max_counts,
        max_genes=max_genes,
        inplace=inplace,
        copy=copy
    )


def filter_genes(
        data: AnnData,
        *,
        min_counts: int | None = None,
        min_cells: int | None = None,
        max_counts: int | None = None,
        max_cells: int | None = None,
        inplace: bool = True,
        copy: bool = False,
) -> AnnData | tuple[np.ndarray, np.ndarray] | None:
    """
    Filter genes based on number of cells or counts.

    Keep genes that have at least `min_counts` counts or are expressed in at least
    `min_cells` cells or have at most `max_counts` counts or are expressed in at
    most `max_cells` cells.

    Note:
        Only provide one of the optional parameters `min_counts`, `min_cells`,
        `max_counts`, `max_cells` per call.

    Parameters:
        data: An annotated data matrix of shape n_obs × n_vars.
        min_counts: Minimum number of counts required for a gene to pass filtering.
        min_cells: Minimum number of cells expressed required for a gene to pass filtering.
        max_counts: Maximum number of counts required for a gene to pass filtering.
        max_cells: Maximum number of cells expressed required for a gene to pass filtering.
        inplace: Perform computation inplace or return result.
        copy: Whether to modify specific copy of object.

    Returns:
        Depending on `inplace`, returns the filtered AnnData object, a tuple of masks,
        or None (updates data in place).
    """
    return pp_module.filter_genes(
        data,
        min_counts=min_counts,
        min_cells=min_cells,
        max_counts=max_counts,
        max_cells=max_cells,
        inplace=inplace,
        copy=copy
    )


def highly_variable_genes(
        adata: AnnData,
        *,
        layer: str | None = None,
        n_top_genes: int | None = None,
        min_disp: float = 0.5,
        max_disp: float = float("inf"),
        min_mean: float = 0.0125,
        max_mean: float = 3,
        span: float = 0.3,
        n_bins: int = 20,
        flavor: str = "seurat",
        subset: bool = False,
        inplace: bool = True,
        batch_key: str | None = None,
        check_values: bool = True,
) -> pd.DataFrame | None:
    """
    Identify highly variable genes (HVGs) using several popular methods
    (Seurat, Cell Ranger, Seurat v3).

    This function computes gene-level variability statistics and marks genes
    that show strong biological variation across cells. HVGs are typically
    used for downstream steps such as PCA, neighbors, clustering, and UMAP.

    Different `flavor` options reproduce the behavior of common pipelines:

    - 'seurat' / 'cell_ranger': dispersion-based HVG selection.
      Requires log-transformed data.
    - 'seurat_v3' / 'seurat_v3_paper': variance-stabilizing transformation (VST).
      Requires raw count data and `scikit-misc`.

    If `batch_key` is provided, HVGs are selected within each batch and then
    merged to avoid batch-specific genes dominating the selection.

    Parameters
    ----------
    adata : AnnData
        Input data matrix (cells × genes).
    layer : str, optional
        Use `adata.layers[layer]` instead of `adata.X`.
    n_top_genes : int, optional
        Number of HVGs to keep. Required for Seurat v3 flavors.
    flavor : {'seurat', 'cell_ranger', 'seurat_v3', 'seurat_v3_paper'}
        Method used to compute HVGs.
    subset : bool
        If True, subset `adata` to HVGs in place.
    inplace : bool
        If True, store results in `adata.var`; otherwise return a DataFrame.
    batch_key : str, optional
        Column in `adata.obs` for batch-aware HVG selection.

    Returns
    -------
    If `inplace=False`, returns a DataFrame with HVG metrics.
    Otherwise, stores the following in `adata.var`:

    - 'highly_variable' : bool indicator of HVGs
    - 'means' : mean expression per gene
    - 'dispersions' / 'dispersions_norm' : for dispersion-based methods
    - 'variances' / 'variances_norm' : for Seurat v3 methods
    - 'highly_variable_rank' : rank of each gene (Seurat v3)
    - 'highly_variable_nbatches' : number of batches where gene is HVG
    - 'highly_variable_intersection' : genes HVG in all batches
    """

    return pp_module.highly_variable_genes(
        adata,
        layer=layer,
        n_top_genes=n_top_genes,
        min_disp=min_disp,
        max_disp=max_disp,
        min_mean=min_mean,
        max_mean=max_mean,
        span=span,
        n_bins=n_bins,
        flavor=flavor,
        subset=subset,
        inplace=inplace,
        batch_key=batch_key,
        check_values=check_values,
    )


def normalization(
        adata: AnnData,
        *,
        target_sum: float | None = None,
        exclude_highly_expressed: bool = False,
        max_fraction: float = 0.05,
        key_added: str | None = None,
        layer: str | None = None,
        inplace: bool = True,
        copy: bool = False,
) -> AnnData | dict[str, np.ndarray] | None:
    """
    Normalize counts per cell.

    Normalize each cell by total counts over all genes, so that every cell has the same
    total count after normalization. If choosing `target_sum=1e6`, this is CPM normalization.

    Parameters:
        adata: The annotated data matrix of shape n_obs × n_vars.
        target_sum: If None, after normalization, each observation (cell) has a total count
                    equal to the median of total counts for observations (cells) before normalization.
        exclude_highly_expressed: Exclude (very) highly expressed genes for the computation
                                  of the normalization factor (size factor) for each cell.
        max_fraction: If `exclude_highly_expressed=True`, consider cells as highly expressed
                      that have more counts than `max_fraction` of the original total counts
                      in at least one cell.
        key_added: Name of the field in `adata.obs` where the normalization factor is stored.
        layer: Layer to normalize instead of `X`. If None, `X` is normalized.
        inplace: Whether to update `adata` or return dictionary with normalized copies of
                 `adata.X` and `adata.layers`.
        copy: Whether to modify copied input object. Not compatible with inplace=False.

    Returns:
        Returns dictionary with normalized copies of `adata.X` and `adata.layers` or updates
        `adata` with normalized version of the original `adata.X` and `adata.layers`,
        depending on `inplace`.
    """
    return pp_module.normalize_total(
        adata,
        target_sum=target_sum,
        exclude_highly_expressed=exclude_highly_expressed,
        max_fraction=max_fraction,
        key_added=key_added,
        layer=layer,
        inplace=inplace,
        copy=copy
    )


def log1p(
        data: AnnData | np.ndarray,
        *,
        base: float | None = None,
        copy: bool = False,
        chunked: bool | None = None,
        chunk_size: int | None = None,
        layer: str | None = None,
        obsm: str | None = None,
) -> AnnData | np.ndarray | None:
    """
    Logarithmize the data matrix.

    Computes X = log(X + 1), where log denotes the natural logarithm unless a different
    base is given.

    Parameters:
        data: The (annotated) data matrix of shape n_obs × n_vars.
        base: Base of the logarithm. Natural logarithm is used by default.
        copy: If an `AnnData` is passed, determines whether a copy is returned.
        chunked: Process the data matrix in chunks, which will save memory. Applies only to `AnnData`.
        chunk_size: n_obs of the chunks to process the data in.
        layer: Entry of layers to transform.
        obsm: Entry of obsm to transform.

    Returns:
        Returns or updates `data`, depending on `copy`.
    """
    return pp_module.log1p(
        data,
        base=base,
        copy=copy,
        chunked=chunked,
        chunk_size=chunk_size,
        layer=layer,
        obsm=obsm
    )


def regress_out(
        adata: AnnData,
        keys: str | list[str],
        *,
        layer: str | None = None,
        n_jobs: int | None = None,
        copy: bool = False,
) -> AnnData | None:
    """
    Regress out (mostly) unwanted sources of variation.

    Uses simple linear regression. Note that this function tends to overcorrect in certain
    circumstances.

    Parameters:
        adata: The annotated data matrix.
        keys: Keys for observation annotation on which to regress on.
        layer: If provided, which element of layers to regress on.
        n_jobs: Number of jobs for parallel computation. None means using `scanpy.settings.n_jobs`.
        copy: Determines whether a copy of `adata` is returned.

    Returns:
        Returns `None` if `copy=False`, else returns an updated `AnnData` object.
        Sets `adata.X` or `adata.layers[layer]` with the corrected count data matrix.
    """
    return pp_module.regress_out(
        adata,
        keys,
        layer=layer,
        n_jobs=n_jobs,
        copy=copy
    )


def scale(
        data: AnnData | np.ndarray,
        *,
        zero_center: bool = True,
        max_value: float | None = None,
        copy: bool = False,
        layer: str | None = None,
        obsm: str | None = None,
        mask_obs: np.ndarray | str | None = None,
) -> AnnData | np.ndarray | None:
    """
    Scale data to unit variance and zero mean.

    Variables (genes) that do not display any variation (are constant across all observations)
    are retained and (for zero_center==True) set to 0 during this operation.

    Parameters:
        data: The (annotated) data matrix of shape n_obs × n_vars.
        zero_center: If False, omit zero-centering variables, which allows to handle sparse input efficiently.
        max_value: Clip (truncate) to this value after scaling. If None, do not clip.
        copy: Whether this function should be performed inplace. If an AnnData object is passed,
              this also determines if a copy is returned.
        layer: If provided, which element of layers to scale.
        obsm: If provided, which element of obsm to scale.
        mask_obs: Restrict both the derivation of scaling parameters and the scaling itself to a
                  certain set of observations.

    Returns:
        Returns `None` if `copy=False`, else returns an updated `AnnData` object.
        Sets `adata.X` (scaled matrix) and `adata.var` (mean, std, var before scaling).
    """
    return pp_module.scale(
        data,
        zero_center=zero_center,
        max_value=max_value,
        copy=copy,
        layer=layer,
        obsm=obsm,
        mask_obs=mask_obs
    )


def pca(
        data: AnnData | np.ndarray,
        n_comps: int | None = None,
        *,
        layer: str | None = None,
        zero_center: bool = True,
        svd_solver: str | None = None,
        chunked: bool = False,
        chunk_size: int | None = None,
        random_state: int | None = 0,
        return_info: bool = False,
        mask_var: np.ndarray | str | None = None,
        use_highly_variable: bool | None = None,
        dtype: str = 'float32',
        key_added: str | None = None,
        copy: bool = False,
) -> AnnData | np.ndarray | None:
    """
    Principal component analysis (PCA).

    Computes PCA coordinates, loadings and variance decomposition.

    Parameters:
        data: The (annotated) data matrix of shape n_obs × n_vars.
        n_comps: Number of principal components to compute. Defaults to 50.
        layer: If provided, which element of layers to use for PCA.
        zero_center: If True, compute PCA from covariance matrix. If False, perform a truncated SVD.
        svd_solver: SVD solver to use ('arpack', 'randomized', 'auto', 'covariance_eigh', etc.).
        chunked: If True, perform an incremental PCA on segments of `chunk_size`.
        chunk_size: Number of observations to include in each chunk.
        random_state: Change to use different initial states for the optimization.
        return_info: Only relevant when not passing an `AnnData`.
        mask_var: To run only on a certain set of genes given by a boolean array or a string.
        use_highly_variable: Whether to use highly variable genes only (Deprecated, use mask_var).
        dtype: Numpy data type string to which to convert the result.
        key_added: If specified, the embedding is stored as `obsm[key_added]`.
        copy: If an `AnnData` is passed, determines whether a copy is returned.

    Returns:
        If `data` is array-like and `return_info=False`, returns the PCA representation.
        Otherwise, returns `None` if `copy=False`, else an updated `AnnData` object.
        Updates `.obsm['X_pca']`, `.varm['PCs']`, and `.uns['pca']`.
    """
    # handle the deprecated _empty logic via default arguments or pass-through
    # we simply pass the arguments to the backend
    return pp_module.pca(
        data,
        n_comps=n_comps,
        layer=layer,
        zero_center=zero_center,
        svd_solver=svd_solver,
        chunked=chunked,
        chunk_size=chunk_size,
        random_state=random_state,
        return_info=return_info,
        mask_var=mask_var,
        use_highly_variable=use_highly_variable,
        dtype=dtype,
        key_added=key_added,
        copy=copy
    )


def sample(
        data: AnnData | np.ndarray,
        fraction: float | None = None,
        *,
        n: int | None = None,
        rng: int | None = None,
        copy: bool = False,
        replace: bool = False,
        axis: str | int = 'obs',
        p: str | np.ndarray | None = None,
) -> AnnData | np.ndarray | tuple | None:
    """
    Sample observations or variables with or without replacement.

    Parameters:
        data: The (annotated) data matrix of shape n_obs × n_vars.
        fraction: Sample to this fraction of the number of observations or variables.
        n: Sample to this number of observations or variables.
        rng: Random seed to change subsampling.
        copy: If an `AnnData` is passed, determines whether a copy is returned.
        replace: If True, samples are drawn with replacement.
        axis: Sample observations (axis 0) or variables (axis 1).
        p: Drawing probabilities (floats) or mask (bools).

    Returns:
        If `isinstance(data, AnnData)` and `copy=False`, this function returns `None`
        (updates inplace). Otherwise returns the subset.
    """
    return pp_module.sample(
        data,
        fraction=fraction,
        n=n,
        rng=rng,
        copy=copy,
        replace=replace,
        axis=axis,
        p=p
    )


from typing import Mapping, Any, Callable


def downsample_counts(
        adata: AnnData,
        counts_per_cell: int | Collection[int] | None = None,
        total_counts: int | None = None,
        *,
        random_state: int | np.random.RandomState | None = 0,
        replace: bool = False,
        copy: bool = False,
) -> AnnData | None:
    """
    Downsample counts from count matrix.

    If `counts_per_cell` is specified, each cell will be downsampled. If `total_counts`
    is specified, the expression matrix will be downsampled to contain at most
    `total_counts`.

    Parameters:
        adata: Annotated data matrix.
        counts_per_cell: Target total counts per cell. If a cell has more than
                         'counts_per_cell', it will be downsampled to this number.
        total_counts: Target total counts. If the count matrix has more than
                      `total_counts` it will be downsampled to have this number.
        random_state: Random seed for subsampling.
        replace: Whether to sample the counts with replacement.
        copy: Determines whether a copy of `adata` is returned.

    Returns:
        Returns `None` if `copy=False`, else returns an `AnnData` object.
        Updates `adata.X` with the downsampled counts matrix.
    """
    return pp_module.downsample_counts(
        adata,
        counts_per_cell=counts_per_cell,
        total_counts=total_counts,
        random_state=random_state,
        replace=replace,
        copy=copy
    )


def recipe_zheng17(
        adata: AnnData,
        *,
        n_top_genes: int = 1000,
        log: bool = True,
        plot: bool = False,
        copy: bool = False,
) -> AnnData | None:
    """
    Normalize and filter as of Zheng et al. [2017].

    Reproduces the preprocessing of Zheng et al. [2017] – the Cell Ranger R Kit of
    10x Genomics. Expects non-logarithmized data.

    Parameters:
        adata: Annotated data matrix.
        n_top_genes: Number of genes to keep.
        log: Take logarithm.
        plot: Show a plot of the gene dispersion vs. mean relation.
        copy: Return a copy of `adata` instead of updating it.

    Returns:
        Returns `None` if `copy=False`, else returns an `AnnData` object.
    """
    return pp_module.recipe_zheng17(
        adata,
        n_top_genes=n_top_genes,
        log=log,
        plot=plot,
        copy=copy
    )


def recipe_weinreb17(
        adata: AnnData,
        *,
        log: bool = True,
        mean_threshold: float = 0.01,
        cv_threshold: float = 2,
        n_pcs: int = 50,
        svd_solver: str = 'randomized',
        random_state: int = 0,
        copy: bool = False,
) -> AnnData | None:
    """
    Normalize and filter as of Weinreb et al. [2017].

    Expects non-logarithmized data.

    Parameters:
        adata: Annotated data matrix.
        log: Logarithmize data?
        mean_threshold: Mean expression threshold.
        cv_threshold: Coefficient of variation threshold.
        n_pcs: Number of Principal Components.
        svd_solver: SVD solver to use.
        random_state: Random seed.
        copy: Return a copy if true.

    Returns:
        Returns `None` if `copy=False`, else returns an `AnnData` object.
    """
    return pp_module.recipe_weinreb17(
        adata,
        log=log,
        mean_threshold=mean_threshold,
        cv_threshold=cv_threshold,
        n_pcs=n_pcs,
        svd_solver=svd_solver,
        random_state=random_state,
        copy=copy
    )


def recipe_seurat(
        adata: AnnData,
        *,
        log: bool = True,
        plot: bool = False,
        copy: bool = False,
) -> AnnData | None:
    """
    Normalize and filter as of Seurat [Satija et al., 2015].

    Expects non-logarithmized data.

    Parameters:
        adata: Annotated data matrix.
        log: Logarithmize data?
        plot: Show a plot of the gene dispersion vs. mean relation.
        copy: Return a copy if true.

    Returns:
        Returns `None` if `copy=False`, else returns an `AnnData` object.
    """
    return pp_module.recipe_seurat(
        adata,
        log=log,
        plot=plot,
        copy=copy
    )


def combat(
        adata: AnnData,
        key: str = 'batch',
        *,
        covariates: Collection[str] | None = None,
        inplace: bool = True,
) -> np.ndarray | None:
    """
    ComBat function for batch effect correction.

    Corrects for batch effects by fitting linear models.

    Parameters:
        adata: Annotated data matrix.
        key: Key to a categorical annotation from `obs` that will be used for batch effect removal.
        covariates: Additional covariates besides the batch variable.
        inplace: Whether to replace `adata.X` or to return the corrected data.

    Returns:
        Returns `numpy.ndarray` if `inplace=False`, else returns `None` and sets `adata.X`.
    """
    return pp_module.combat(
        adata,
        key=key,
        covariates=covariates,
        inplace=inplace
    )


def scrublet(
        adata: AnnData,
        adata_sim: AnnData | None = None,
        *,
        batch_key: str | None = None,
        sim_doublet_ratio: float = 2.0,
        expected_doublet_rate: float = 0.05,
        stdev_doublet_rate: float = 0.02,
        synthetic_doublet_umi_subsampling: float = 1.0,
        knn_dist_metric: str | Callable[[np.ndarray, np.ndarray], float] = 'euclidean',
        normalize_variance: bool = True,
        log_transform: bool = False,
        mean_center: bool = True,
        n_prin_comps: int = 30,
        use_approx_neighbors: bool | None = None,
        get_doublet_neighbor_parents: bool = False,
        n_neighbors: int | None = None,
        threshold: float | None = None,
        verbose: bool = True,
        copy: bool = False,
        random_state: int | np.random.RandomState | None = 0,
) -> AnnData | None:
    """
    Predict doublets using Scrublet.

    Predict cell doublets using a nearest-neighbor classifier of observed transcriptomes
    and simulated doublets.

    Parameters:
        adata: The annotated data matrix.
        adata_sim: Optional AnnData object generated by `scrublet_simulate_doublets`.
        batch_key: Optional `obs` column name discriminating between batches.
        sim_doublet_ratio: Number of doublets to simulate relative to the number of
                           observed transcriptomes.
        expected_doublet_rate: Estimated doublet rate for the experiment.
        stdev_doublet_rate: Uncertainty in the expected doublet rate.
        synthetic_doublet_umi_subsampling: Rate for sampling UMIs when creating synthetic doublets.
        knn_dist_metric: Distance metric used when finding nearest neighbors.
        normalize_variance: If True, normalize the data such that each gene has a variance of 1.
        log_transform: Whether to use `log1p()` to log-transform the data prior to PCA.
        mean_center: If True, center the data such that each gene has a mean of 0.
        n_prin_comps: Number of principal components used to embed the transcriptomes.
        use_approx_neighbors: Use approximate nearest neighbor method (annoy) for the KNN classifier.
        get_doublet_neighbor_parents: If True, return parent transcriptomes of doublet neighbors.
        n_neighbors: Number of neighbors used to construct the KNN graph.
        threshold: Doublet score threshold.
        verbose: If True, log progress updates.
        copy: If True, return a copy of the input `adata` with Scrublet results added.
        random_state: Initial state for doublet simulation and nearest neighbors.

    Returns:
        If `copy=True` it returns a new AnnData object, else adds fields to `adata`.
        Updates `.obs['doublet_score']`, `.obs['predicted_doublet']`, and `.uns['scrublet']`.
    """
    return pp_module.scrublet(
        adata,
        adata_sim=adata_sim,
        batch_key=batch_key,
        sim_doublet_ratio=sim_doublet_ratio,
        expected_doublet_rate=expected_doublet_rate,
        stdev_doublet_rate=stdev_doublet_rate,
        synthetic_doublet_umi_subsampling=synthetic_doublet_umi_subsampling,
        knn_dist_metric=knn_dist_metric,
        normalize_variance=normalize_variance,
        log_transform=log_transform,
        mean_center=mean_center,
        n_prin_comps=n_prin_comps,
        use_approx_neighbors=use_approx_neighbors,
        get_doublet_neighbor_parents=get_doublet_neighbor_parents,
        n_neighbors=n_neighbors,
        threshold=threshold,
        verbose=verbose,
        copy=copy,
        random_state=random_state
    )


def scrublet_simulate_doublets(
        adata: AnnData,
        *,
        layer: str | None = None,
        sim_doublet_ratio: float = 2.0,
        synthetic_doublet_umi_subsampling: float = 1.0,
        random_seed: int = 0,
) -> AnnData:
    """
    Simulate doublets by adding the counts of random observed transcriptome pairs.

    Parameters:
        adata: The annotated data matrix.
        layer: Layer of adata where raw values are stored, or 'X' if values are in .X.
        sim_doublet_ratio: Number of doublets to simulate relative to the number of
                           observed transcriptomes.
        synthetic_doublet_umi_subsampling: Rate for sampling UMIs when creating synthetic doublets.
        random_seed: Random seed.

    Returns:
        AnnData with simulated doublets in .X.
        Adds fields to `.obsm['scrublet']['doublet_parents']` and `.uns['scrublet']['parameters']`.
    """
    return pp_module.scrublet_simulate_doublets(
        adata,
        layer=layer,
        sim_doublet_ratio=sim_doublet_ratio,
        synthetic_doublet_umi_subsampling=synthetic_doublet_umi_subsampling,
        random_seed=random_seed
    )


def neighbors(
        adata: AnnData,
        n_neighbors: int = 15,
        n_pcs: int | None = None,
        *,
        use_rep: str | None = None,
        knn: bool = True,
        method: str = 'umap',
        transformer: Any | None = None,
        metric: str | Callable[[np.ndarray, np.ndarray], float] = 'euclidean',
        metric_kwds: Mapping[str, Any] = {},
        random_state: int | np.random.RandomState | None = 0,
        key_added: str | None = None,
        copy: bool = False,
) -> AnnData | None:
    """
    Compute the nearest neighbors distance matrix and a neighborhood graph of observations.

    Parameters:
        adata: Annotated data matrix.
        n_neighbors: The size of local neighborhood (in terms of number of neighboring data points).
        n_pcs: Use this many PCs.
        use_rep: Use the indicated representation. 'X' or any key for .obsm is valid.
        knn: If True, use a hard threshold to restrict the number of neighbors to `n_neighbors`.
        method: Use 'umap' or 'gauss' for computing connectivities.
        transformer: Approximate kNN search implementation.
        metric: A known metric’s name or a callable that returns a distance.
        metric_kwds: Options for the metric.
        random_state: A numpy random seed.
        key_added: If specified, the neighbors data is added to .uns[key_added].
        copy: Return a copy instead of writing to adata.

    Returns:
        Returns `None` if `copy=False`, else returns an `AnnData` object.
        Updates `.obsp['distances']`, `.obsp['connectivities']`, and `.uns['neighbors']`.
    """
    return pp_module.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep=use_rep,
        knn=knn,
        method=method,
        transformer=transformer,
        metric=metric,
        metric_kwds=metric_kwds,
        random_state=random_state,
        key_added=key_added,
        copy=copy
    )
