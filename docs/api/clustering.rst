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
   
   Lotus uses the cplearn algorithm for clustering, which can:
   - Auto-detect the best data representation
   - Handle large-scale datasets
   - Generate stable clustering results

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
