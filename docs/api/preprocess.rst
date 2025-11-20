Preprocessing Module
=====================

The ``lotus.workflows.preprocess`` module provides data preprocessing related functions.

Overview
--------

Preprocessing is the first step in single-cell data analysis, including:

- **Quality Control (QC)**: Calculate quality metrics
- **Filtering**: Remove low-quality cells and genes
- **Normalization**: Eliminate technical variation
- **Highly Variable Genes (HVG) Selection**: Select genes for subsequent analysis
- **Scaling**: Standardize data
- **Principal Component Analysis (PCA)**: Dimensionality reduction
- **Neighbor Graph Construction**: Prepare for clustering and visualization

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
   
   In single-cell sequencing, low-quality cells may be caused by:
   - Cell death or damage
   - Doublets (multiple cells)
   - Sequencing technical issues
   
   By calculating metrics such as total counts and number of expressed genes, these low-quality cells can be identified.

Filtering
~~~~~~~~~

.. autofunction:: lotus.workflows.preprocess.preprocess.filtering

   Filter low-quality cells and genes based on QC metrics.

   **Biological Background:**
   
   Filtering is an important step to remove technical noise:
   - Filter low-quality cells: Remove dead cells, doublets, etc.
   - Filter low-expression genes: Remove genes that are not expressed in any cells

Normalization
~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.preprocess.preprocess.normalization

   Normalize data to eliminate the effects of sequencing depth.

   **Biological Background:**
   
   Different cells may have different sequencing depths. Normalization can:
   - Eliminate sequencing depth differences
   - Make expression levels comparable between different cells
   - Usually normalize to 10,000 counts per cell

Highly Variable Genes Selection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.preprocess.preprocess.hvg

   Select highly variable genes for subsequent analysis.

   **Biological Background:**
   
   Highly Variable Genes (HVG) are:
   - Genes with large expression variation between different cells
   - Usually include cell type-specific genes
   - Used for dimensionality reduction and clustering analysis

Scaling
~~~~~~~

.. autofunction:: lotus.workflows.preprocess.preprocess.scaling

   Scale data (z-score standardization).

   **Biological Background:**
   
   Scaling can:
   - Put expression levels of different genes on the same scale
   - Avoid high-expression genes dominating the analysis results
   - Usually used before PCA

Principal Component Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.preprocess.preprocess.pca

   Perform principal component analysis for dimensionality reduction.

   **Biological Background:**
   
   PCA (Principal Component Analysis):
   - Projects high-dimensional data into low-dimensional space
   - Preserves the main variation in the data
   - Used for subsequent clustering and visualization

Neighbor Graph Construction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.preprocess.preprocess.neighbors

   Build neighbor graph between cells.

   **Biological Background:**
   
   Neighbor graph:
   - Represents similarity between cells
   - Used for clustering algorithms (such as Leiden)
   - Used for dimensionality reduction algorithms like UMAP
