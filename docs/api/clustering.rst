Clustering Module
=================

The ``lotus.workflows.clustering`` module provides cell clustering analysis functionality.

Overview
--------

Clustering analysis is used to group cells based on gene expression patterns and identify different cell types or states.

Main Functions
--------------

.. autofunction:: lotus.workflows.clustering.clustering.clustering

   **Biological Background:**
   
   Cell clustering is a core step in single-cell analysis:
   - Group cells based on gene expression patterns
   - Identify different cell types or states
   - Provide basis for subsequent differential expression analysis
   
   **Lotus supports multiple clustering methods:**
   
   - **cplearn** (default): Lotus cplearn clustering algorithm - auto-detect the best data representation, handle large-scale datasets, generate stable clustering results
   - **leiden**: Scanpy Leiden algorithm - fast and scalable community detection
   - **louvain**: Scanpy Louvain algorithm - classic community detection method
   
   You can switch between methods using the ``method`` parameter. All methods output cluster labels in scanpy-compatible format.

   **Usage Examples:**

   .. code-block:: python

      from lotus.workflows import clustering
      
      # Option 1: Use cplearn (default)
      model = clustering(
          adata,
          method="cplearn",  # or omit for default
          use_rep="X_latent",
          key_added="cplearn_labels",
          cluster_resolution=1.2,
      )
      
      # Option 2: Use scanpy Leiden algorithm
      clustering(
          adata,
          method="leiden",
          cluster_resolution=0.5,
          key_added="leiden",  # auto-set if None
      )
      
      # Option 3: Use scanpy Louvain algorithm
      clustering(
          adata,
          method="louvain",
          cluster_resolution=0.5,
          key_added="louvain",  # auto-set if None
      )
      
      # View clustering results (works for all methods)
      print(adata.obs["cplearn_labels"].value_counts())
      print(adata.obs["leiden"].value_counts())
      print(adata.obs["louvain"].value_counts())

   **Parameter Description:**

   - ``method``: Clustering method to use. Options: ``"cplearn"`` (default), ``"leiden"``, ``"louvain"``
   
   - ``use_rep``: Data representation to use (cplearn only)
     - ``"X_latent"``: Latent representation (recommended)
     - ``"X_pca"``: PCA representation
     - ``"X"``: Raw data
     - ``None``: Auto-detect (default)
   
   - ``key_added``: Key name for cluster labels in ``adata.obs``
     - If ``None``, uses method-specific default: ``"cplearn_labels"``, ``"leiden"``, or ``"louvain"``
   
   - ``cluster_resolution``: Clustering resolution (applies to all methods)
     - Smaller values (0.5-1.0): Fewer, larger clusters
     - Larger values (1.5-2.0): More, smaller clusters
     - Default: 1.2
   
   - ``stable_core_frac``: Stable core fraction (cplearn only), default 0.25
   
   - ``stable_ng_num``: Number of neighbors for stable core (cplearn only), default 8
   
   - ``fine_grained``: Whether to use fine-grained clustering (cplearn only), default False
   
   - ``propagate``: Whether to propagate labels (cplearn only), default True
   
   - ``random_state``: Random seed for reproducibility (scanpy methods only), default 0
   
   - ``neighbors_key``: Key in ``adata.uns`` where neighbors are stored (scanpy methods only)
   
   - ``obsp``: Key in ``adata.obsp`` to use as adjacency matrix (scanpy methods only)
   
   - ``print_summary``: Whether to print cluster summary, default True

.. autofunction:: lotus.workflows.clustering.clustering.run_clustering

   This is an alias for ``clustering()`` for backward compatibility.

Helper Functions
----------------

.. autofunction:: lotus.workflows.clustering.clustering.summarize_clusters

