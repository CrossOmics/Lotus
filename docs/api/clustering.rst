Clustering Module
=================

The ``lotus.workflows.clustering`` module provides cell clustering analysis functionality.

Overview
--------

Clustering analysis is used to group cells based on gene expression patterns and identify different cell types or states.

**Citation:** 

* **Leiden Algorithm** (default): V. A. Traag, L. Waltman, and N. J. van Eck. *From Louvain to Leiden: guaranteeing well-connected communities.* Scientific Reports 9, 5233 (2019). `doi:10.1038/s41598-019-41695-z <https://doi.org/10.1038/s41598-019-41695-z>`_

* **Louvain Algorithm**: V. D. Blondel, J.-L. Guillaume, R. Lambiotte, and E. Lefebvre. *Fast unfolding of communities in large networks.* Journal of Statistical Mechanics: Theory and Experiment 2008, P10008 (2008). `doi:10.1088/1742-5468/2008/10/P10008 <https://doi.org/10.1088/1742-5468/2008/10/P10008>`_

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

