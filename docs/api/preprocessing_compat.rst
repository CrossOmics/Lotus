Preprocessing Module (lotus.pp)
================================

The ``lotus.pp`` module provides preprocessing functions for single-cell data analysis.

Overview
--------

``lotus.pp`` serves as a building block API that provides direct access to all preprocessing functions. This module is used internally by ``lotus.workflows`` to implement high-level preprocessing functionality, and can also be used directly when you need fine-grained control or want to access advanced features not covered by the workflows module, such as batch correction and doublet detection.

**Key Features:**

The module provides complete preprocessing functionality with standard function signatures and behavior. All preprocessing functions are available. For standard preprocessing workflows, we recommend using ``lotus.workflows.preprocess()`` which provides a complete preprocessing pipeline with optimized defaults. For advanced features like batch correction and doublet detection, use the building blocks API directly.

Compatibility
-------------

Usage example:

.. code-block:: python

   import lotus as lt
   
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

   ⚙️ Calculate quality control metrics for cells and genes (also available as ``lotus.workflows.qc()``).

Filtering
~~~~~~~~~

.. autofunction:: lotus.methods.scanpy.preprocessing.filter_cells

   ⚙️ Filter cells based on QC metrics (also available as ``lotus.workflows.filtering()``).

.. autofunction:: lotus.methods.scanpy.preprocessing.filter_genes

   ⚙️ Filter genes based on expression (also available as ``lotus.workflows.filtering()``).

Normalization
~~~~~~~~~~~~~

.. autofunction:: lotus.methods.scanpy.preprocessing.normalize_total

   ⚙️ Normalize counts per cell (also available as ``lotus.workflows.normalization()``).

.. autofunction:: lotus.methods.scanpy.preprocessing.log1p

   Logarithmize the data matrix.

Highly Variable Genes
~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.methods.scanpy.preprocessing.highly_variable_genes

   ⚙️ Identify highly variable genes (also available as ``lotus.workflows.hvg()``).

Scaling
~~~~~~~

.. autofunction:: lotus.methods.scanpy.preprocessing.scale

   ⚙️ Scale data to unit variance and zero mean (also available as ``lotus.workflows.scaling()``).

Dimensionality Reduction
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.methods.scanpy.preprocessing.pca

   ⚙️ Principal component analysis (also available as ``lotus.workflows.pca()``).

Neighbor Graph
~~~~~~~~~~~~~~

.. autofunction:: lotus.methods.scanpy.preprocessing.neighbors

   ⚙️ Compute neighborhood graph of cells (also available as ``lotus.workflows.neighbors()``).

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

All preprocessing functions are available in ``lotus.pp``. Functions marked with ⚙️ are also available as high-level wrappers in ``lotus.workflows``, which we recommend for standard preprocessing workflows. For advanced features like batch correction and doublet detection, use the building blocks API directly.

**Quality Control Functions:**

The ``calculate_qc_metrics()`` function calculates quality control metrics for cells and genes (also available as ``lotus.workflows.qc()``).

**Filtering Functions:**

The ``filter_cells()`` function filters cells based on QC metrics (also available as ``lotus.workflows.filtering()``). The ``filter_genes()`` function filters genes based on expression (also available as ``lotus.workflows.filtering()``). The ``filter_genes_dispersion()`` function filters genes by dispersion.

**Normalization Functions:**

The ``normalize_total()`` function normalizes counts per cell (also available as ``lotus.workflows.normalization()``). The ``normalize_per_cell()`` function normalizes per cell but is deprecated. The ``log1p()`` function applies log transformation to the data matrix.

**Highly Variable Genes Functions:**

The ``highly_variable_genes()`` function identifies highly variable genes (also available as ``lotus.workflows.hvg()``).

**Scaling Functions:**

The ``scale()`` function scales data to unit variance and zero mean (also available as ``lotus.workflows.scaling()``).

**Dimensionality Reduction Functions:**

The ``pca()`` function performs principal component analysis (also available as ``lotus.workflows.pca()``). The ``neighbors()`` function computes neighborhood graph of cells (also available as ``lotus.workflows.neighbors()``).

**Batch Correction Functions:**

The ``combat()`` function performs ComBat batch correction for removing batch effects. This advanced feature is not covered by the workflows module and should be used directly from the building blocks API. The ``regress_out()`` function regresses out unwanted sources of variation.

**Doublet Detection Functions:**

The ``scrublet()`` function detects doublets using Scrublet. This advanced feature is not covered by the workflows module and should be used directly from the building blocks API. The ``scrublet_simulate_doublets()`` function simulates doublets for training.

**Other Functions:**

The ``downsample_counts()`` function downsamples counts to a target number. The ``sample()`` function samples cells or genes. The ``subsample()`` function subsamples data. The ``sqrt()`` function applies square root transformation. The ``recipe_seurat()`` function provides Seurat preprocessing recipe. The ``recipe_zheng17()`` function provides Zheng et al. 2017 preprocessing recipe. The ``recipe_weinreb17()`` function provides Weinreb et al. 2017 preprocessing recipe.

For detailed documentation of each function, please refer to the `scanpy official documentation <https://scanpy.readthedocs.io/en/stable/api/scanpy.pp.html>`__.

