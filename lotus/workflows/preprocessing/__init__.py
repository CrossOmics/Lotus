from __future__ import annotations

import numpy as np
from anndata import AnnData

import lotus as lt


def input(adata: AnnData, *, save_raw: bool = True, raw_layer: str = "raw_counts") -> None:
    """
    Input step: Save raw count matrix
    
    Parameters:
        adata: AnnData object
        save_raw: Whether to save the raw count matrix
        raw_layer: Name of the layer to save raw count matrix
    """
    if save_raw:
        adata.layers[raw_layer] = adata.X.copy()


def qc(
    adata: AnnData,
    *,
    expr_type: str = "counts",
    var_type: str = "genes",
    qc_vars: tuple = (),
    percent_top: tuple[int, ...] | None = None,
    layer: str | None = None,
    use_raw: bool = False,
    inplace: bool = True,
    log1p: bool = True,
) -> None:
    """
    QC (Quality Control) step: Calculate QC metrics
    
    Parameters:
        adata: AnnData object
        expr_type: Expression type for QC calculation
        var_type: Variable type for QC calculation
        qc_vars: QC variables to calculate
        percent_top: Percent top genes to calculate. If None, auto-adjusts based on n_vars
        layer: Layer to use for QC calculation
        use_raw: Whether to use raw data
        inplace: Whether to modify adata in place
        log1p: Whether to apply log1p transformation
    """
    # Auto-adjust percent_top based on number of genes
    if percent_top is None:
        n_vars = adata.n_vars
        percent_top = tuple(min(p, n_vars) for p in (50, 100, 200, 500) if p <= n_vars)
        if not percent_top:
            percent_top = (min(50, n_vars),)
    
    lt.pp.calculate_qc_metrics(
        adata,
        expr_type=expr_type,
        var_type=var_type,
        qc_vars=qc_vars,
        percent_top=percent_top,
        layer=layer,
        use_raw=use_raw,
        inplace=inplace,
        log1p=log1p,
    )


def filtering(
    adata: AnnData,
    *,
    min_counts: int | None = None,
    min_genes: int | None = None,
    max_counts: int | None = None,
    max_genes: int | None = None,
    min_cells: int | None = None,
    inplace: bool = True,
) -> None:
    """
    Filtering step: Filter cells and genes
    
    Parameters:
        adata: AnnData object
        min_counts: Minimum number of counts per cell
        min_genes: Minimum number of genes per cell
        max_counts: Maximum number of counts per cell
        max_genes: Maximum number of genes per cell
        min_cells: Minimum number of cells expressing a gene
        inplace: Whether to modify adata in place
    """
    # Filter cells
    if min_counts is not None or min_genes is not None or max_counts is not None or max_genes is not None:
        lt.pp.filter_cells(
            adata,
            min_counts=min_counts,
            min_genes=min_genes,
            max_counts=max_counts,
            max_genes=max_genes,
            inplace=inplace,
        )
    
    # Filter genes
    if min_cells is not None:
        lt.pp.filter_genes(
            adata,
            min_cells=min_cells,
            inplace=inplace,
        )


def normalization(adata: AnnData, *, target_sum: float = 1e4) -> None:
    """
    Normalization step: Normalize to target sum
    
    Parameters:
        adata: AnnData object
        target_sum: Target sum for normalization
    """
    # Normalize total counts (this is safe to repeat)
    lt.pp.normalize_total(adata, target_sum=target_sum)
    
    # Apply log1p transformation
    # If data is already log-transformed, scanpy will warn but we'll catch it
    # and handle gracefully to prevent crashes
    import warnings
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:
            lt.pp.log1p(adata)
        except Exception as e:
            # If there's an error, check if it's related to log transformation
            if "log" in str(e).lower() or "invalid value" in str(e).lower():
                print(f"Warning: Skipping log1p transformation - data may already be log-transformed. Error: {e}")
            else:
                raise
        
        # Check for warnings about already log-transformed data
        if w:
            for warning in w:
                if "already log-transformed" in str(warning.message).lower():
                    print("Warning: Data appears to be already log-transformed. Skipping log1p transformation.")
                    # Try to reverse the log1p that was just applied
                    # Since we applied log1p to already log-transformed data, we need expm1 to reverse
                    try:
                        if hasattr(adata.X, 'data'):  # Sparse matrix
                            adata.X.data = np.expm1(adata.X.data)
                        else:  # Dense array
                            adata.X = np.expm1(adata.X)
                    except Exception:
                        # If we can't reverse, that's okay - the data might still be usable
                        pass
                    break


def hvg(adata: AnnData, *, n_top_genes: int | None = None) -> None:
    """
    HVG (Highly Variable Genes) step: Select highly variable genes
    
    Parameters:
        adata: AnnData object
        n_top_genes: Number of highly variable genes. If None, uses min(2000, adata.n_vars)
    """
    if n_top_genes is None:
        n_top_genes = min(2000, adata.n_vars)
    lt.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor="seurat",
        subset=True,
    )


def scaling(adata: AnnData, *, zero_center: bool = True, max_value: float = 10) -> None:
    """
    Scaling step: Standardize data
    
    Parameters:
        adata: AnnData object
        zero_center: Whether to zero center
        max_value: Maximum value for clipping
    """
    lt.pp.scale(adata, zero_center=zero_center, max_value=max_value)


def pca(adata: AnnData, *, n_pcs: int | None = None, svd_solver: str = "arpack") -> None:
    """
    PCA step: Principal component analysis
    
    Parameters:
        adata: AnnData object
        n_pcs: Number of principal components. If None, automatically determined by scanpy
               (default: min(50, n_obs, n_vars))
        svd_solver: SVD solver to use ('arpack', 'randomized', 'auto')
    """
    lt.tl.pca(adata, n_comps=n_pcs, svd_solver=svd_solver)
    
    # Plot PCA variance ratio (standard Scanpy workflow)
    # This visualizes the variance explained by each principal component
    # Use n_pcs if provided, otherwise use default 30 for plotting
    plot_n_pcs = n_pcs if n_pcs is not None else 30
    lt.pl.pca_variance_ratio(adata, n_pcs=plot_n_pcs, log=True, show=False)
    
    # Combine initial latent vectors with PCA results (if exists)
    latent_init = adata.obsm.pop("X_latent_init", None)
    if latent_init is not None:
        combined = np.hstack([adata.obsm["X_pca"], latent_init])
    else:
        combined = adata.obsm["X_pca"]
    adata.obsm["X_latent"] = combined.astype(np.float32, copy=False)


def neighbors(adata: AnnData, *, use_rep: str = "X_pca", n_neighbors: int = 15) -> None:
    """
    Neighbors step: Construct neighbor graph
    
    Parameters:
        adata: AnnData object
        use_rep: Representation to use for neighbor graph construction, default is "X_pca"
        n_neighbors: Number of neighbors
    """
    lt.pp.neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors)


def log1p(adata: AnnData, *, base: float | None = None, layer: str | None = None) -> None:
    """
    Log1p step: Log transform data (log(x + 1))
    
    Parameters:
        adata: AnnData object
        base: Base of logarithm. If None, uses natural logarithm
        layer: Layer to transform. If None, transforms X
    """
    lt.pp.log1p(adata, base=base, layer=layer)


def regress_out(
    adata: AnnData,
    keys: str | list[str],
    *,
    layer: str | None = None,
    n_jobs: int | None = None,
) -> None:
    """
    Regress out step: Regress out unwanted variables (e.g., cell cycle, batch effects)
    
    Parameters:
        adata: AnnData object
        keys: Keys of observations to regress out
        layer: Layer to use for regression. If None, uses X
        n_jobs: Number of jobs for parallel processing
    """
    lt.pp.regress_out(adata, keys, layer=layer, n_jobs=n_jobs)


def combat(
    adata: AnnData,
    key: str = "batch",
    *,
    covariates: list[str] | None = None,
) -> None:
    """
    ComBat step: Batch effect correction using ComBat
    
    Parameters:
        adata: AnnData object
        key: Key in adata.obs that contains batch information
        covariates: Additional covariates to preserve
    """
    lt.pp.combat(adata, key, covariates=covariates)


def scrublet(
    adata: AnnData,
    *,
    batch_key: str | None = None,
    sim_doublet_ratio: float = 2.0,
    expected_doublet_rate: float = 0.05,
    threshold: float | None = None,
    random_state: int = 0,
) -> None:
    """
    Scrublet step: Predict doublets using Scrublet
    
    Parameters:
        adata: AnnData object
        batch_key: Key in adata.obs that contains batch information
        sim_doublet_ratio: Ratio of simulated doublets to observed cells
        expected_doublet_rate: Expected doublet rate
        threshold: Doublet score threshold. If None, automatically determined
        random_state: Random seed
    """
    lt.pp.scrublet(
        adata,
        batch_key=batch_key,
        sim_doublet_ratio=sim_doublet_ratio,
        expected_doublet_rate=expected_doublet_rate,
        threshold=threshold,
        random_state=random_state,
    )


def scrublet_simulate_doublets(
    adata: AnnData,
    *,
    layer: str | None = None,
    sim_doublet_ratio: float = 2.0,
    random_seed: int = 0,
) -> AnnData:
    """
    Scrublet simulate doublets step: Simulate doublets by adding counts of random pairs
    
    Parameters:
        adata: AnnData object
        layer: Layer to use for simulation. If None, uses X
        sim_doublet_ratio: Ratio of simulated doublets to observed cells
        random_seed: Random seed
    
    Returns:
        AnnData object with simulated doublets
    """
    return lt.pp.scrublet_simulate_doublets(
        adata,
        layer=layer,
        sim_doublet_ratio=sim_doublet_ratio,
        random_seed=random_seed,
    )


def sample(
    adata: AnnData,
    *,
    fraction: float | None = None,
    n: int | None = None,
    axis: str = "obs",
    replace: bool = False,
    random_state: int = 0,
) -> AnnData:
    """
    Sample step: Sample observations or variables
    
    Parameters:
        adata: AnnData object
        fraction: Fraction of data to sample (0-1)
        n: Number of items to sample
        axis: Axis to sample along ('obs' or 'var')
        replace: Whether to sample with replacement
        random_state: Random seed
    
    Returns:
        Sampled AnnData object
    """
    return lt.pp.sample(
        adata,
        fraction=fraction,
        n=n,
        axis=axis,
        replace=replace,
        rng=random_state,
    )


def downsample_counts(
    adata: AnnData,
    *,
    counts_per_cell: int | list[int] | None = None,
    total_counts: int | None = None,
    random_state: int = 0,
) -> AnnData:
    """
    Downsample counts step: Downsample counts to a target number
    
    Parameters:
        adata: AnnData object
        counts_per_cell: Target counts per cell (can be a list for different targets per cell)
        total_counts: Target total counts
        random_state: Random seed
    
    Returns:
        AnnData object with downsampled counts
    """
    return lt.pp.downsample_counts(
        adata,
        counts_per_cell=counts_per_cell,
        total_counts=total_counts,
        random_state=random_state,
    )


def recipe_zheng17(
    adata: AnnData,
    *,
    n_top_genes: int = 1000,
    log: bool = True,
) -> None:
    """
    Recipe Zheng17 step: Normalize and filter as of Zheng et al. 2017
    
    Parameters:
        adata: AnnData object
        n_top_genes: Number of highly variable genes
        log: Whether to apply log transformation
    """
    lt.pp.recipe_zheng17(adata, n_top_genes=n_top_genes, log=log)


def recipe_weinreb17(
    adata: AnnData,
    *,
    log: bool = True,
    mean_threshold: float = 0.01,
    cv_threshold: int = 2,
    n_pcs: int = 50,
    random_state: int = 0,
) -> None:
    """
    Recipe Weinreb17 step: Normalize and filter as of Weinreb et al. 2017
    
    Parameters:
        adata: AnnData object
        log: Whether to apply log transformation
        mean_threshold: Mean expression threshold
        cv_threshold: Coefficient of variation threshold
        n_pcs: Number of principal components
        random_state: Random seed
    """
    lt.pp.recipe_weinreb17(
        adata,
        log=log,
        mean_threshold=mean_threshold,
        cv_threshold=cv_threshold,
        n_pcs=n_pcs,
        random_state=random_state,
    )


def recipe_seurat(
    adata: AnnData,
    *,
    log: bool = True,
) -> None:
    """
    Recipe Seurat step: Normalize and filter as of Seurat (Satija et al. 2015)
    
    Parameters:
        adata: AnnData object
        log: Whether to apply log transformation
    """
    lt.pp.recipe_seurat(adata, log=log)


def preprocess(
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
) -> None:
    """
    Complete preprocessing pipeline following standard Scanpy workflow:
    QC → Filtering → Normalization → Save Raw → HVG → Scaling → PCA → Neighbors
    
    Standard Scanpy preprocessing steps:
    1. Calculate QC metrics (n_counts, n_genes, pct_mt, etc.)
    2. Filter cells and genes
    3. Normalize total counts and apply log1p transformation
    4. Save raw data (normalized + log1p) to adata.raw
    5. Select highly variable genes (HVG) and subset to HVG only
    6. Scale data (zero-center and clip)
    7. Principal component analysis (PCA)
    8. Compute neighborhood graph
    
    Parameters:
        adata: AnnData object
        n_pcs: Number of principal components for PCA. If None, automatically determined by scanpy
               (default: min(50, n_obs, n_vars))
        target_sum: Target sum for normalization (default: 1e4)
        n_top_genes: Number of highly variable genes. If None, uses min(2000, adata.n_vars)
        n_neighbors: Number of neighbors for neighbor graph construction
        use_rep: Representation to use for neighbor graph construction, default is "X_pca"
        save_raw: Whether to save the normalized+log1p data to adata.raw (standard Scanpy practice)
        raw_layer: Name of the layer to save raw count matrix (if save_raw=True, also saves to layer)
        min_genes: Minimum number of genes per cell for filtering
        min_cells: Minimum number of cells expressing a gene for filtering
        min_counts: Minimum number of counts per cell for filtering
        max_counts: Maximum number of counts per cell for filtering
        max_genes: Maximum number of genes per cell for filtering
        pct_mt_max: Maximum percentage of mitochondrial genes (filters cells with pct_mt > pct_mt_max)
        hvg_flavor: Flavor for highly variable genes selection ('seurat_v3', 'seurat', 'cell_ranger')
        batch_key: Key in adata.obs for batch information (for batch correction)
        regress_out_keys: Keys in adata.obs to regress out (e.g., ['total_counts', 'pct_counts_mt'])
        use_combat: Whether to use ComBat for batch effect correction (requires batch_key)
    """
    # Check if data has already been preprocessed and restore from raw if needed
    # This prevents errors when re-running preprocessing on already processed data
    if adata.raw is not None or (save_raw and raw_layer and raw_layer in adata.layers):
        # Data has been preprocessed before, restore from raw counts
        if save_raw and raw_layer and raw_layer in adata.layers:
            # Restore from raw_layer (original count matrix)
            print("Restoring from raw counts layer for re-preprocessing...")
            adata.X = adata.layers[raw_layer].copy()
            # Clear preprocessing results that will be recomputed
            if 'X_pca' in adata.obsm:
                del adata.obsm['X_pca']
            if 'X_latent' in adata.obsm:
                del adata.obsm['X_latent']
            if 'neighbors' in adata.uns:
                del adata.uns['neighbors']
            if 'pca' in adata.uns:
                del adata.uns['pca']
            # Clear obs columns that will be recomputed
            if 'highly_variable' in adata.var:
                del adata.var['highly_variable']
            if 'highly_variable_nbatches' in adata.var:
                del adata.var['highly_variable_nbatches']
            if 'highly_variable_intersection' in adata.var:
                del adata.var['highly_variable_intersection']
            if 'means' in adata.var:
                del adata.var['means']
            if 'dispersions' in adata.var:
                del adata.var['dispersions']
            if 'dispersions_norm' in adata.var:
                del adata.var['dispersions_norm']
        elif adata.raw is not None:
            # Restore from adata.raw (normalized+log1p data)
            # We need to reverse the normalization and log1p to get back to counts
            print("Restoring from adata.raw for re-preprocessing...")
            # adata.raw contains normalized+log1p data, so we need to reverse:
            # 1. Reverse log1p: expm1
            # 2. Reverse normalization: multiply by original total counts
            # But we don't have the original total counts, so we'll use raw_layer if available
            if save_raw and raw_layer and raw_layer in adata.layers:
                adata.X = adata.layers[raw_layer].copy()
            else:
                # Try to reverse from adata.raw
                # This is approximate and may not be exact
                print("Warning: Cannot fully restore from adata.raw. Using adata.raw.X as starting point.")
                adata.X = adata.raw.X.copy()
                # Clear preprocessing results
                if 'X_pca' in adata.obsm:
                    del adata.obsm['X_pca']
                if 'X_latent' in adata.obsm:
                    del adata.obsm['X_latent']
                if 'neighbors' in adata.uns:
                    del adata.uns['neighbors']
                if 'pca' in adata.uns:
                    del adata.uns['pca']
    
    # Save raw count matrix at the beginning (before any transformations)
    if save_raw and raw_layer:
        if raw_layer not in adata.layers:
            adata.layers[raw_layer] = adata.X.copy()
    
    # Step 1: QC - Calculate quality control metrics
    qc(adata)
    
    # Step 2: Filtering - Filter cells and genes
    filtering(
        adata,
        min_genes=min_genes,
        min_cells=min_cells,
        min_counts=min_counts,
        max_counts=max_counts,
        max_genes=max_genes,
    )
    
    # Additional filtering: mitochondrial percentage
    if pct_mt_max is not None:
        pct_mt_key = None
        if 'pct_counts_mt' in adata.obs.columns:
            pct_mt_key = 'pct_counts_mt'
        elif 'pct_mt' in adata.obs.columns:
            pct_mt_key = 'pct_mt'
        
        if pct_mt_key is not None:
            n_before = adata.n_obs
            # Create boolean mask and filter
            mask = adata.obs[pct_mt_key] < pct_mt_max
            # Use AnnData's internal method for in-place subsetting
            # This is the standard way to filter observations in place
            adata._inplace_subset_obs(mask)
            if n_before != adata.n_obs:
                print(f"Filtered {n_before - adata.n_obs} cells with {pct_mt_key} >= {pct_mt_max}")
    
    # Step 3: Normalization - Normalize total counts and log1p transform
    normalization(adata, target_sum=target_sum)
    
    # Step 4: Save raw data (standard Scanpy practice: save normalized+log1p to adata.raw)
    # This is the standard Scanpy workflow: adata.raw contains normalized+log1p data
    if save_raw:
        adata.raw = adata
    
    # Optional: Regress out unwanted variables (e.g., total_counts, pct_counts_mt)
    if regress_out_keys is not None:
        regress_out(adata, regress_out_keys)
    
    # Optional: Batch effect correction using ComBat
    if use_combat and batch_key is not None:
        if batch_key in adata.obs.columns:
            combat(adata, key=batch_key)
        else:
            print(f"Warning: batch_key '{batch_key}' not found in adata.obs. Skipping ComBat correction.")
    
    # Step 5: HVG - Select highly variable genes and subset
    if n_top_genes is None:
        n_top_genes = min(2000, adata.n_vars)
    
    # Use the specified flavor for HVG selection
    if hvg_flavor == "seurat_v3":
        lt.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            flavor="seurat_v3",
            subset=True,
        )
    else:
        lt.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            flavor=hvg_flavor,
            subset=True,
        )
    
    # Step 6: Scaling - Zero-center and clip
    scaling(adata)
    
    # Step 7: PCA - Principal component analysis
    pca(adata, n_pcs=n_pcs)
    
    # Step 8: Neighbors - Compute neighborhood graph
    neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors)


__all__ = [
    "qc",
    "filtering",
    "normalization",
    "hvg",
    "scaling",
    "pca",
    "neighbors",
    "preprocess",
    "log1p",
    "regress_out",
    "combat",
    "scrublet",
    "scrublet_simulate_doublets",
    "sample",
    "downsample_counts",
    "recipe_zheng17",
    "recipe_weinreb17",
    "recipe_seurat",
]
