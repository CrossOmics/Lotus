Core Selection
==============

Core selection is used to identify core cells and compute core map embedding, which can improve clustering stability and provide additional insights.

Overview
--------

Core selection serves two purposes in the Lotus workflow:

1. **Before Clustering** (optional): Identify core cells that can help improve clustering stability
2. **After Clustering** (for cplearn workflow): Compute core map embedding using the clustering model

The core selection step helps:
- Identify stable, high-quality cells
- Improve clustering results
- Provide additional dimensionality reduction representation

Core Selection Before Clustering
----------------------------------

In some workflows, core selection can be performed before clustering to identify core cells:

.. code-block:: python

    from lotus.workflows import core_selection
    import scanpy as sc
    
    # First, ensure neighbors graph exists
    sc.pp.neighbors(adata, n_neighbors=15)
    
    # Core selection before clustering (if supported by your workflow)
    # Note: This typically requires additional setup depending on your analysis needs

Core Selection After Clustering (cplearn workflow)
---------------------------------------------------

For the cplearn clustering workflow, core selection is typically performed after clustering:

.. code-block:: python

    from lotus.workflows import clustering, core_selection
    
    # First perform clustering
    model = clustering(
        adata,
        use_rep="X_latent",
        key_added="cplearn_labels",
        cluster_resolution=1.2,
    )
    
    # Then compute core map embedding
    core_selection(
        adata,
        model=model,
        use_rep="X_latent",
        key_added="X_cplearn_coremap",
    )
    
    # View results
    print(adata.obsm["X_cplearn_coremap"].shape)

Parameter Description
---------------------

- ``model``: CorespectModel object from clustering step (required for post-clustering core selection)
- ``use_rep``: Data representation to use
  - ``"X_latent"``: Latent representation (recommended)
  - ``"X_pca"``: PCA representation
  - ``"X"``: Raw data
  - ``None``: Auto-detect (default)

- ``key_added``: Key name for embedding results in ``adata.obsm``, default is ``"X_cplearn_coremap"``

Compatibility with scanpy
-------------------------

Core selection results are fully compatible with scanpy:

.. code-block:: python

    # Core map embedding can be used with scanpy visualization
    import scanpy as sc
    sc.pl.umap(adata, color="cplearn_labels", obsm="X_cplearn_coremap")

Next Steps
----------

After core selection, you can:
- Perform visualization: :doc:`visualization`
- Perform differential expression analysis: :doc:`deg`

