from __future__ import annotations

import numpy as np
from anndata import AnnData

import lotus as lt


def input(adata: AnnData, *, save_raw: bool = True, raw_layer: str = "raw_counts") -> None:
    """
    Save raw count matrix.
    
    Parameters:
        adata (AnnData): AnnData object
        save_raw (bool): Save raw count matrix. Default: True
        raw_layer (str): Layer name for raw counts. Default: "raw_counts"
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
    Calculate QC metrics.
    
    Parameters:
        adata (AnnData): AnnData object
        expr_type (str): Expression type for QC calculation. Default: "counts"
        var_type (str): Variable type for QC calculation. Default: "genes"
        qc_vars (tuple): QC variables to calculate. Default: ()
        percent_top (tuple[int, ...] | None): Percent top genes to calculate. Default: None (auto-adjust)
        layer (str | None): Layer to use for QC calculation. Default: None
        use_raw (bool): Use raw data. Default: False
        inplace (bool): Modify adata in place. Default: True
        log1p (bool): Apply log1p transformation. Default: True
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
    Filter cells and genes.
    
    Parameters:
        adata (AnnData): AnnData object
        min_counts (int | None): Minimum number of counts per cell. Default: None
        min_genes (int | None): Minimum number of genes per cell. Default: None
        max_counts (int | None): Maximum number of counts per cell. Default: None
        max_genes (int | None): Maximum number of genes per cell. Default: None
        min_cells (int | None): Minimum number of cells expressing a gene. Default: None
        inplace (bool): Modify adata in place. Default: True
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
    Normalize to target sum.
    
    Parameters:
        adata (AnnData): AnnData object
        target_sum (float): Target sum for normalization. Default: 1e4
    """
    lt.pp.normalize_total(adata, target_sum=target_sum)
    lt.pp.log1p(adata)


def hvg(adata: AnnData, *, n_top_genes: int | None = None) -> None:
    """
    Select highly variable genes.
    
    Parameters:
        adata (AnnData): AnnData object
        n_top_genes (int | None): Number of highly variable genes. Default: None (min(2000, n_vars))
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
    Standardize data.
    
    Parameters:
        adata (AnnData): AnnData object
        zero_center (bool): Zero center the data. Default: True
        max_value (float): Maximum value for clipping. Default: 10
    """
    lt.pp.scale(adata, zero_center=zero_center, max_value=max_value)


def pca(adata: AnnData, *, n_pcs: int = 20) -> None:
    """
    Principal component analysis.
    
    Parameters:
        adata (AnnData): AnnData object
        n_pcs (int): Number of principal components. Default: 20
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
    Construct neighbor graph.
    
    Parameters:
        adata (AnnData): AnnData object
        use_rep (str): Representation to use for neighbor graph. Default: "X_pca"
        n_neighbors (int): Number of neighbors. Default: 15
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
    
    Parameters:
        adata (AnnData): AnnData object
        n_pcs (int): Number of principal components. Default: 20
        target_sum (float): Target sum for normalization. Default: 1e4
        n_top_genes (int | None): Number of highly variable genes. Default: None (min(2000, n_vars))
        n_neighbors (int): Number of neighbors. Default: 15
        use_rep (str): Representation for neighbor graph. Default: "X_pca"
        save_raw (bool): Save raw count matrix. Default: True
        raw_layer (str): Layer name for raw counts. Default: "raw_counts"
        min_genes (int | None): Minimum genes per cell. Default: None
        min_cells (int | None): Minimum cells per gene. Default: None
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
