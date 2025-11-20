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
   
   - **Lotus cplearn** (default): Auto-detect the best data representation, handle large-scale datasets, generate stable clustering results
   - **scanpy** (alternative): Use proven scanpy algorithms (Leiden or Louvain), fully compatible with Lotus workflow
   
   You can freely switch between methods based on your needs.

   **Usage Example:**

   .. code-block:: python

      from lotus.workflows import clustering
      
      # Basic usage
      model = clustering(
          adata,
          use_rep="X_latent",
          key_added="cplearn_labels",
          cluster_resolution=1.2,
      )
      
      # View clustering results
      print(adata.obs["cplearn_labels"].value_counts())

   **Switching to scanpy:**
   
   If you prefer to use scanpy's clustering algorithms, you can easily switch:
   
   .. code-block:: python
   
      import scanpy as sc
      
      # Option 1: Use scanpy Leiden algorithm
      sc.tl.leiden(adata, resolution=0.5, key_added="leiden")
      
      # Option 2: Use scanpy Louvain algorithm
      sc.tl.louvain(adata, resolution=0.5, key_added="louvain")
      
      # Use the cluster key in subsequent Lotus functions
      # e.g., marker_genes(adata, cluster_key="leiden")

   **Parameter Description:**

   - ``use_rep``: Data representation to use
     - ``"X_latent"``: Latent representation (recommended)
     - ``"X_pca"``: PCA representation
     - ``"X"``: Raw data
     - ``None``: Auto-detect
   
   - ``cluster_resolution``: Clustering resolution
     - Smaller values (0.5-1.0): Fewer, larger clusters
     - Larger values (1.5-2.0): More, smaller clusters
     - Default: 1.2

.. autofunction:: lotus.workflows.clustering.clustering.run_clustering

   This is an alias for ``clustering()`` for backward compatibility.

Helper Functions
----------------

.. autofunction:: lotus.workflows.clustering.clustering.summarize_clusters

   Summarize cluster label statistics.
