Preprocessing Workflow
======================

Preprocessing is the first step in single-cell data analysis, including quality control, data filtering, normalization, and dimensionality reduction.

Overview
--------

The goals of preprocessing are:
- Identify and remove low-quality cells (such as dead cells, doublets)
- Filter low-expression or highly variable genes
- Normalize data to eliminate technical variation
- Select highly variable genes for subsequent analysis
- Perform dimensionality reduction (PCA) and build neighbor graph

Complete Preprocessing Pipeline
-------------------------------

The simplest way is to use the ``preprocess()`` function, which automatically executes all preprocessing steps:

.. code-block:: python

    from lotus.workflows import preprocess
    
    preprocess(
        adata,
        n_pcs=20,              # Number of PCA principal components
        target_sum=1e4,        # Normalization target sum
        n_top_genes=2000,       # Number of highly variable genes
        n_neighbors=15,        # Number of neighbors for neighbor graph
        save_raw=True,          # Save raw count matrix
    )

Step-by-Step Preprocessing
--------------------------

If you need finer control, you can call each step separately:

Quality Control (QC)
~~~~~~~~~~~~~~~~~~~~

Calculate quality control metrics to help identify low-quality cells:

.. code-block:: python

    from lotus.workflows import qc
    
    qc(adata)
    
    # View QC metrics
    print(adata.obs.columns)  # Will include total_counts, n_genes_by_counts, etc.

Filtering
~~~~~~~~~

Filter low-quality cells and genes based on QC metrics:

.. code-block:: python

    from lotus.workflows import filtering
    
    filtering(
        adata,
        min_genes=200,      # Minimum number of genes expressed per cell
        min_cells=3,       # Minimum number of cells expressing each gene
    )

Normalization
~~~~~~~~~~~~~

Normalize data to eliminate the effects of sequencing depth:

.. code-block:: python

    from lotus.workflows import normalization
    
    normalization(adata, target_sum=1e4)  # Normalize to 10,000 counts per cell

Highly Variable Genes (HVG)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Select highly variable genes for subsequent analysis:

.. code-block:: python

    from lotus.workflows import hvg
    
    hvg(adata, n_top_genes=2000)  # Select top 2000 highly variable genes

Scaling
~~~~~~~

Scale data (z-score standardization):

.. code-block:: python

    from lotus.workflows import scaling
    
    scaling(adata, zero_center=True, max_value=10)

Principal Component Analysis (PCA)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Perform principal component analysis for dimensionality reduction:

.. code-block:: python

    from lotus.workflows import pca
    
    pca(adata, n_pcs=20)  # Compute 20 principal components

Building Neighbor Graph
~~~~~~~~~~~~~~~~~~~~~~~

Build neighbor graph between cells for subsequent clustering and UMAP:

.. code-block:: python

    from lotus.workflows import neighbors
    
    neighbors(adata, use_rep="X_pca", n_neighbors=15)

Parameter Description
---------------------

- ``n_pcs``: Number of PCA principal components, typically 20-50
- ``target_sum``: Normalization target sum, typically 10,000
- ``n_top_genes``: Number of highly variable genes, typically 2000-3000
- ``n_neighbors``: Number of neighbors for neighbor graph, typically 15-30
- ``min_genes``: Threshold for filtering low-quality cells
- ``min_cells``: Threshold for filtering low-expression genes

Frequently Asked Questions
---------------------------

**Q: How to choose appropriate filtering thresholds?**
A: Usually determined by viewing the distribution plots of QC metrics. You can use scanpy's plotting functions to visualize these metrics.

**Q: How many principal components should be selected?**
A: Usually 20-50 principal components are sufficient. This can be determined by viewing the variance explained ratio of PCA.

**Q: How to choose the number of highly variable genes?**
A: Usually 2000-3000 highly variable genes. If the dataset is very large, this can be increased to 4000-5000.

Next Steps
----------

After completing preprocessing, you can:
- Perform clustering analysis: :doc:`clustering`
- Perform visualization: :doc:`visualization`
