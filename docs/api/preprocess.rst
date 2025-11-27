Preprocessing Module
=====================

The ``lotus.workflows.preprocess`` module provides data preprocessing related functions.

Overview
--------

Preprocessing is the first step in single-cell data analysis, encompassing quality control (QC) to calculate quality metrics, filtering to remove low-quality cells and genes, normalization to eliminate technical variation, highly variable genes (HVG) selection for subsequent analysis, scaling to standardize data, principal component analysis (PCA) for dimensionality reduction, and neighbor graph construction to prepare for clustering and visualization.

**Citation:** Preprocessing functions are based on `scanpy <https://scanpy.readthedocs.io/>`_. Please cite:

* **scanpy**: F. Alexander Wolf, Philipp Angerer, and Fabian J. Theis. *SCANPY: large-scale single-cell gene expression data analysis.* Genome Biology 19, 15 (2018). `doi:10.1186/s13059-017-1382-0 <https://doi.org/10.1186/s13059-017-1382-0>`_

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

Filtering
~~~~~~~~~

.. autofunction:: lotus.workflows.preprocess.preprocess.filtering

   Filter low-quality cells and genes based on QC metrics.

Normalization
~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.preprocess.preprocess.normalization

   Normalize data to eliminate the effects of sequencing depth.

Highly Variable Genes Selection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.preprocess.preprocess.hvg

   Select highly variable genes for subsequent analysis.

Scaling
~~~~~~~

.. autofunction:: lotus.workflows.preprocess.preprocess.scaling

   Scale data (z-score standardization).

Principal Component Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.preprocess.preprocess.pca

   Perform principal component analysis for dimensionality reduction.

Neighbor Graph Construction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.preprocess.preprocess.neighbors

   Build neighbor graph between cells.
