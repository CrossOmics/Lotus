Clustering Module
=================

The ``lotus.workflows.clustering`` module provides cell clustering analysis functionality.

Overview
--------

Clustering analysis is used to group cells based on gene expression patterns and identify different cell types or states.

Main Functions
--------------

.. autofunction:: lotus.workflows.clustering.clustering.clustering

   **Usage Examples:**

   .. code-block:: python

      from lotus.workflows import clustering
      
      # Use Leiden algorithm (default)
      clustering(adata, method="leiden", cluster_resolution=0.5)
      
      # Use Louvain algorithm
      clustering(adata, method="louvain", cluster_resolution=0.5)
      
      # View clustering results
      print(adata.obs["leiden"].value_counts())

Helper Functions
----------------

.. autofunction:: lotus.workflows.clustering.clustering.summarize_clusters

