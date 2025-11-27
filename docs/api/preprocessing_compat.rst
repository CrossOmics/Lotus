Preprocessing Module (lotus.pp)
================================

The ``lotus.pp`` module provides scanpy-compatible preprocessing functions. This module is a complete wrapper around ``scanpy.pp``, providing all scanpy preprocessing functions with identical interfaces.

Overview
--------

``lotus.pp`` is a direct wrapper around ``scanpy.pp``. All functions work exactly the same way as in scanpy, making Lotus fully compatible with scanpy preprocessing workflows.

**Key Features:**

The module provides complete scanpy.pp compatibility with the same function signatures and behavior. It can be used as a drop-in replacement for scanpy.pp, and all scanpy preprocessing functions are available.

Compatibility
-------------

You can use ``lotus.pp`` exactly like ``scanpy.pp``:

.. code-block:: python

   import lotus as lt
   
   # These work exactly like scanpy.pp
   lt.pp.calculate_qc_metrics(adata)
   lt.pp.filter_cells(adata, min_genes=200)
   lt.pp.normalize_total(adata, target_sum=1e4)
   lt.pp.log1p(adata)
   lt.pp.highly_variable_genes(adata, n_top_genes=2000)
   lt.pp.scale(adata)
   lt.pp.pca(adata, n_comps=50)
   lt.pp.neighbors(adata, n_neighbors=15)

Main Functions
--------------

Quality Control
~~~~~~~~~~~~~~~

.. autofunction:: lotus.methods.scanpy.preprocessing.calculate_qc_metrics

   Calculate quality control metrics for cells and genes.

Filtering
~~~~~~~~~

.. autofunction:: lotus.methods.scanpy.preprocessing.filter_cells

   Filter cells based on QC metrics.

.. autofunction:: lotus.methods.scanpy.preprocessing.filter_genes

   Filter genes based on expression.

Normalization
~~~~~~~~~~~~~

.. autofunction:: lotus.methods.scanpy.preprocessing.normalize_total

   Normalize counts per cell.

.. autofunction:: lotus.methods.scanpy.preprocessing.log1p

   Logarithmize the data matrix.

Highly Variable Genes
~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.methods.scanpy.preprocessing.highly_variable_genes

   Identify highly variable genes.

Scaling
~~~~~~~

.. autofunction:: lotus.methods.scanpy.preprocessing.scale

   Scale data to unit variance and zero mean.

Dimensionality Reduction
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.methods.scanpy.preprocessing.pca

   Principal component analysis.

Neighbor Graph
~~~~~~~~~~~~~~

.. autofunction:: lotus.methods.scanpy.preprocessing.neighbors

   Compute neighborhood graph of cells.

Batch Correction
~~~~~~~~~~~~~~~~

.. autofunction:: lotus.methods.scanpy.preprocessing.combat

   ComBat batch correction for removing batch effects.

Doublet Detection
~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.methods.scanpy.preprocessing.scrublet

   Detect doublets using Scrublet.

Other Functions
~~~~~~~~~~~~~~~~

.. autofunction:: lotus.methods.scanpy.preprocessing.regress_out

   Regress out unwanted sources of variation.

.. autofunction:: lotus.methods.scanpy.preprocessing.downsample_counts

   Downsample counts to a target number.

Usage Examples
--------------

Complete Preprocessing Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import lotus as lt
   
   # Calculate QC metrics
   lt.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
   
   # Filter cells and genes
   lt.pp.filter_cells(adata, min_genes=200)
   lt.pp.filter_genes(adata, min_cells=3)
   
   # Normalize and log transform
   lt.pp.normalize_total(adata, target_sum=1e4)
   lt.pp.log1p(adata)
   
   # Find highly variable genes
   lt.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
   
   # Keep only highly variable genes
   adata.raw = adata
   adata = adata[:, adata.var.highly_variable]
   
   # Scale and PCA
   lt.pp.scale(adata, max_value=10)
   lt.pp.pca(adata, n_comps=50)
   
   # Compute neighbors
   lt.pp.neighbors(adata, n_neighbors=15, n_pcs=50)

Batch Correction
~~~~~~~~~~~~~~~~

.. code-block:: python

   import lotus as lt
   
   # Apply ComBat batch correction
   lt.pp.combat(adata, key='batch')

Doublet Detection
~~~~~~~~~~~~~~~~~

.. code-block:: python

   import lotus as lt
   
   # Detect doublets
   lt.pp.scrublet(adata)

Complete Function List
----------------------

All functions from ``scanpy.pp`` are available in ``lotus.pp``:

- ``calculate_qc_metrics()`` - Calculate QC metrics
- ``filter_cells()`` - Filter cells
- ``filter_genes()`` - Filter genes
- ``filter_genes_dispersion()`` - Filter genes by dispersion
- ``normalize_total()`` - Normalize counts
- ``normalize_per_cell()`` - Normalize per cell (deprecated)
- ``log1p()`` - Log transform
- ``sqrt()`` - Square root transform
- ``highly_variable_genes()`` - Find HVG
- ``scale()`` - Scale data
- ``pca()`` - Principal component analysis
- ``neighbors()`` - Compute neighbor graph
- ``combat()`` - Batch correction
- ``regress_out()`` - Regress out variables
- ``scrublet()`` - Doublet detection
- ``scrublet_simulate_doublets()`` - Simulate doublets
- ``downsample_counts()`` - Downsample counts
- ``sample()`` - Sample cells/genes
- ``subsample()`` - Subsample data
- ``recipe_seurat()`` - Seurat preprocessing recipe
- ``recipe_zheng17()`` - Zheng et al. 2017 preprocessing recipe
- ``recipe_weinreb17()`` - Weinreb et al. 2017 preprocessing recipe

For detailed documentation of each function, please refer to the `scanpy official documentation <https://scanpy.readthedocs.io/en/stable/api/scanpy.pp.html>`__.

Note
----

All functions in ``lotus.pp`` are direct wrappers around ``scanpy.pp`` functions. They have identical signatures, parameters, and behavior. You can use them as drop-in replacements for scanpy.pp functions.

