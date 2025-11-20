Clustering Analysis
===================

Clustering analysis is used to identify cell types and subpopulations, and is one of the core steps in single-cell data analysis.

Overview
--------

The goals of clustering analysis are:
- Group cells based on gene expression patterns
- Identify different cell types or states
- Provide basis for subsequent differential expression analysis

Lotus uses the cplearn algorithm for clustering, which can:
- Auto-detect the best data representation (X_latent, X_pca, or X)
- Handle large-scale datasets
- Generate stable clustering results

Basic Usage
-----------

.. code-block:: python

    from lotus.workflows import clustering
    
    model = clustering(
        adata,
        use_rep="X_latent",        # Data representation to use (optional, auto-detected)
        key_added="cplearn_labels", # Key name for storing cluster labels
        cluster_resolution=1.2,    # Clustering resolution
        stable_core_frac=0.25,     # Stable core fraction
        stable_ng_num=8,           # Number of neighbors for stable core
    )
    
    # View clustering results
    print(adata.obs["cplearn_labels"].value_counts())

Parameter Description
---------------------

- ``use_rep``: Data representation to use
  - ``"X_latent"``: Latent representation (recommended if available)
  - ``"X_pca"``: PCA representation
  - ``"X"``: Raw data
  - ``None``: Auto-detect (default)

- ``key_added``: Key name for storing cluster labels in ``adata.obs``, default is ``"cplearn_labels"``

- ``cluster_resolution``: Clustering resolution, controls clustering granularity
  - Smaller values (0.5-1.0): Produce fewer, larger clusters
  - Larger values (1.5-2.0): Produce more, smaller clusters
  - Default: 1.2

- ``stable_core_frac``: Fraction of stable core, default 0.25

- ``stable_ng_num``: Number of neighbors for stable core, default 8

- ``fine_grained``: Whether to use fine-grained clustering, default False

- ``propagate``: Whether to propagate labels, default True

Adjusting Clustering Resolution
-------------------------------

Clustering resolution is a key parameter that controls clustering granularity:

.. code-block:: python

    # Coarser clustering (fewer clusters)
    model_coarse = clustering(adata, cluster_resolution=0.8)
    
    # Medium clustering (default)
    model_medium = clustering(adata, cluster_resolution=1.2)
    
    # Finer clustering (more clusters)
    model_fine = clustering(adata, cluster_resolution=1.8)

Viewing Clustering Results
---------------------------

After clustering is complete, results are stored in ``adata.obs``:

.. code-block:: python

    # View cluster labels
    print(adata.obs["cplearn_labels"])
    
    # Count cells per cluster
    print(adata.obs["cplearn_labels"].value_counts())
    
    # View number of clusters
    print(f"Number of clusters: {adata.obs['cplearn_labels'].nunique()}")

Compatibility with scanpy
-------------------------

Lotus clustering results are fully compatible with scanpy:

.. code-block:: python

    # Cluster labels are automatically converted to categorical type
    print(type(adata.obs["cplearn_labels"]))  # <class 'pandas.core.arrays.categorical.Categorical'>
    
    # Can be used for scanpy visualization
    import scanpy as sc
    sc.pl.umap(adata, color="cplearn_labels")

Frequently Asked Questions
---------------------------

**Q: How to choose appropriate clustering resolution?**
A: Usually start with 1.0-1.5, adjust based on biological knowledge and clustering results. You can cluster with different resolutions and compare through visualization.

**Q: What if clustering results are unstable?**
A: You can adjust the ``stable_core_frac`` and ``stable_ng_num`` parameters, or set ``fine_grained=True``.

**Q: How to know if clustering is reasonable?**
A: View cluster distribution in reduced-dimensional space through UMAP visualization, and validate biological differences between clusters through differential expression analysis.

Next Steps
----------

After completing clustering, you can:
- Perform differential expression analysis: :doc:`deg`
- Perform visualization: :doc:`visualization`
