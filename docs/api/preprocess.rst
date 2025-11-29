Preprocessing Module
=====================

The ``lotus.workflows.preprocessing`` module provides data preprocessing related functions.

Overview
--------

Preprocessing is the first step in single-cell data analysis, encompassing quality control (QC) to calculate quality metrics, filtering to remove low-quality cells and genes, normalization to eliminate technical variation, highly variable genes (HVG) selection for subsequent analysis, scaling to standardize data, principal component analysis (PCA) for dimensionality reduction, and neighbor graph construction to prepare for clustering and visualization.

**Citation:** Preprocessing functions are based on `scanpy <https://scanpy.readthedocs.io/>`_. Please cite:

* **scanpy**: F. Alexander Wolf, Philipp Angerer, and Fabian J. Theis. *SCANPY: large-scale single-cell gene expression data analysis.* Genome Biology 19, 15 (2018). `doi:10.1186/s13059-017-1382-0 <https://doi.org/10.1186/s13059-017-1382-0>`_

Complete Preprocessing Pipeline
--------------------------------

.. autofunction:: lotus.workflows.preprocessing.preprocess

   **Usage Example:**

   .. code-block:: python

      from lotus.workflows import preprocessing
      import scanpy as sc
      
      # Load data
      adata = sc.datasets.pbmc3k()
      
      # Complete preprocessing pipeline
      preprocessing.preprocess(
          adata,
          n_pcs=20,
          target_sum=1e4,
          n_top_genes=2000,
          n_neighbors=15,
      )

Step-by-Step Preprocessing Functions
-------------------------------------

Input
~~~~~

.. autofunction:: lotus.workflows.preprocessing.input

Quality Control
~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.preprocessing.qc

Filtering
~~~~~~~~~

.. autofunction:: lotus.workflows.preprocessing.filtering

Normalization
~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.preprocessing.normalization

Highly Variable Genes Selection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.preprocessing.hvg

Scaling
~~~~~~~

.. autofunction:: lotus.workflows.preprocessing.scaling

Principal Component Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.preprocessing.pca

Neighbor Graph Construction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.preprocessing.neighbors

Additional Preprocessing Functions
----------------------------------

Log Transformation
~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.preprocessing.log1p

Regress Out Variables
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.preprocessing.regress_out

Batch Effect Correction
~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.preprocessing.combat

Doublet Detection
~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.preprocessing.scrublet

.. autofunction:: lotus.workflows.preprocessing.scrublet_simulate_doublets

Sampling and Downsampling
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.preprocessing.sample

.. autofunction:: lotus.workflows.preprocessing.downsample_counts

Preprocessing Recipes
~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.preprocessing.recipe_zheng17

.. autofunction:: lotus.workflows.preprocessing.recipe_weinreb17

.. autofunction:: lotus.workflows.preprocessing.recipe_seurat
