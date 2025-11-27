Clustering Module
=================

The ``lotus.workflows.clustering`` module provides cell clustering analysis functionality.

Overview
--------

Clustering analysis is used to group cells based on gene expression patterns and identify different cell types or states.

Main Functions
--------------

.. autofunction:: lotus.workflows.clustering.clustering.clustering

   **Lotus supports multiple clustering methods:**
   
   The **cplearn** method (default) is Lotus cplearn clustering algorithm that auto-detects the best data representation, handles large-scale datasets, and generates stable clustering results. The **leiden** method uses the Scanpy Leiden algorithm for fast and scalable community detection. The **louvain** method uses the Scanpy Louvain algorithm, a classic community detection method. You can switch between methods using the ``method`` parameter. All methods output cluster labels in scanpy-compatible format.

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

   The ``method`` parameter specifies the clustering method to use, with options ``"cplearn"`` (default), ``"leiden"``, or ``"louvain"``. The ``use_rep`` parameter specifies the data representation to use (cplearn only), with options ``"X_latent"`` (latent representation, recommended), ``"X_pca"`` (PCA representation), ``"X"`` (raw data), or ``None`` (auto-detect, default). The ``key_added`` parameter specifies the key name for cluster labels in ``adata.obs``. If ``None``, it uses method-specific defaults: ``"cplearn_labels"``, ``"leiden"``, or ``"louvain"``. The ``cluster_resolution`` parameter controls clustering resolution (applies to all methods), where smaller values (0.5-1.0) produce fewer and larger clusters, larger values (1.5-2.0) produce more and smaller clusters, with default 1.2. The ``stable_core_frac`` parameter specifies stable core fraction (cplearn only), default 0.25. The ``stable_ng_num`` parameter specifies the number of neighbors for stable core (cplearn only), default 8. The ``fine_grained`` parameter controls whether to use fine-grained clustering (cplearn only), default False. The ``propagate`` parameter controls whether to propagate labels (cplearn only), default True. The ``random_state`` parameter sets the random seed for reproducibility (scanpy methods only), default 0. The ``neighbors_key`` parameter specifies the key in ``adata.uns`` where neighbors are stored (scanpy methods only). The ``obsp`` parameter specifies the key in ``adata.obsp`` to use as adjacency matrix (scanpy methods only). The ``print_summary`` parameter controls whether to print cluster summary, default True.

Helper Functions
----------------

.. autofunction:: lotus.workflows.clustering.clustering.summarize_clusters

