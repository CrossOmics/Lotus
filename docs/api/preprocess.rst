Preprocessing Module
=====================

The ``lotus.workflows.preprocess`` module provides data preprocessing related functions.

Overview
--------

Preprocessing is the first step in single-cell data analysis, encompassing quality control (QC) to calculate quality metrics, filtering to remove low-quality cells and genes, normalization to eliminate technical variation, highly variable genes (HVG) selection for subsequent analysis, scaling to standardize data, principal component analysis (PCA) for dimensionality reduction, and neighbor graph construction to prepare for clustering and visualization.

Complete Preprocessing Pipeline
--------------------------------

.. autofunction:: lotus.workflows.preprocess.preprocess.preprocess

   **Usage Example:**

   .. code-block:: python

      from lotus.workflows import preprocess
      import scanpy as sc
      
      # Load data
      adata = sc.datasets.pbmc3k()
      
      # Complete preprocessing pipeline
      preprocess(
          adata,
          n_pcs=20,
          target_sum=1e4,
          n_top_genes=2000,
          n_neighbors=15,
      )

Step-by-Step Preprocessing Functions
-------------------------------------

Quality Control
~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.preprocess.preprocess.qc

   Calculate quality control metrics to help identify low-quality cells.

   **Biological Background:**
   
   In single-cell sequencing, low-quality cells may be caused by cell death or damage, doublets (multiple cells), or sequencing technical issues. By calculating metrics such as total counts and number of expressed genes, these low-quality cells can be identified.

Filtering
~~~~~~~~~

.. autofunction:: lotus.workflows.preprocess.preprocess.filtering

   Filter low-quality cells and genes based on QC metrics.

   **Biological Background:**
   
   Filtering is an important step to remove technical noise by filtering low-quality cells such as dead cells and doublets, and by filtering low-expression genes that are not expressed in any cells.

Normalization
~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.preprocess.preprocess.normalization

   Normalize data to eliminate the effects of sequencing depth.

   **Biological Background:**
   
   Different cells may have different sequencing depths. Normalization eliminates sequencing depth differences and makes expression levels comparable between different cells, usually normalizing to 10,000 counts per cell.

Highly Variable Genes Selection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.preprocess.preprocess.hvg

   Select highly variable genes for subsequent analysis.

   **Biological Background:**
   
   Highly Variable Genes (HVG) are genes with large expression variation between different cells, usually including cell type-specific genes. They are used for dimensionality reduction and clustering analysis.

Scaling
~~~~~~~

.. autofunction:: lotus.workflows.preprocess.preprocess.scaling

   Scale data (z-score standardization).

   **Biological Background:**
   
   Scaling puts expression levels of different genes on the same scale and avoids high-expression genes dominating the analysis results. It is usually used before PCA.

Principal Component Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.preprocess.preprocess.pca

   Perform principal component analysis for dimensionality reduction.

   **Biological Background:**
   
   PCA (Principal Component Analysis) projects high-dimensional data into low-dimensional space while preserving the main variation in the data. It is used for subsequent clustering and visualization.

Neighbor Graph Construction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.preprocess.preprocess.neighbors

   Build neighbor graph between cells.

   **Biological Background:**
   
   The neighbor graph represents similarity between cells and is used for clustering algorithms such as Leiden and dimensionality reduction algorithms like UMAP.
