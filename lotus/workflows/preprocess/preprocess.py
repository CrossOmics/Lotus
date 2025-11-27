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
    lt.pp.normalize_total(adata, target_sum=target_sum)
    lt.pp.log1p(adata)


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


def pca(adata: AnnData, *, n_pcs: int = 20) -> None:
    """
    PCA step: Principal component analysis
    
    Parameters:
        adata: AnnData object
        n_pcs: Number of principal components
    """
    lt.tl.pca(adata, n_comps=n_pcs)
    
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


def preprocess(
    adata: AnnData,
    *,
    n_pcs: int = 20,
    target_sum: float = 1e4,
    n_top_genes: int | None = None,
    n_neighbors: int = 15,
    use_rep: str = "X_pca",
    save_raw: bool = True,
    raw_layer: str = "raw_counts",
    min_genes: int | None = None,
    min_cells: int | None = None,
) -> None:
    """
    Complete preprocessing pipeline: QC, filtering, normalization, HVG selection, scaling, PCA, and neighbors.
    
    This function performs a complete preprocessing pipeline including quality control, filtering, normalization, highly variable genes selection, scaling, principal component analysis, and neighbor graph construction. Raw count matrix is automatically saved at the beginning if save_raw=True.
    
    Parameters:
        adata: AnnData object
        n_pcs: Number of principal components for PCA
        target_sum: Target sum for normalization
        n_top_genes: Number of highly variable genes. If None, uses min(2000, adata.n_vars)
        n_neighbors: Number of neighbors for neighbor graph construction
        use_rep: Representation to use for neighbor graph construction, default is "X_pca"
        save_raw: Whether to save the raw count matrix
        raw_layer: Name of the layer to save raw count matrix
        min_genes: Minimum number of genes per cell for filtering
        min_cells: Minimum number of cells expressing a gene for filtering
    """
    # Save raw count matrix at the beginning
    if save_raw:
        adata.layers[raw_layer] = adata.X.copy()
    
    qc(adata)
    filtering(adata, min_genes=min_genes, min_cells=min_cells)
    normalization(adata, target_sum=target_sum)
    hvg(adata, n_top_genes=n_top_genes)
    scaling(adata)
    pca(adata, n_pcs=n_pcs)
    neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors)
