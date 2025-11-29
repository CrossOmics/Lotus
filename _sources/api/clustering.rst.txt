Clustering Module
=================

The ``lotus.workflows.clustering`` module provides cell clustering analysis functionality.

Overview
--------

Clustering analysis is used to group cells based on gene expression patterns and identify different cell types or states.

**Citation:** 

* **Leiden Algorithm** (default): V. A. Traag, L. Waltman, and N. J. van Eck. *From Louvain to Leiden: guaranteeing well-connected communities.* Scientific Reports 9, 5233 (2019). `doi:10.1038/s41598-019-41695-z <https://doi.org/10.1038/s41598-019-41695-z>`_

* **Louvain Algorithm**: V. D. Blondel, J.-L. Guillaume, R. Lambiotte, and E. Lefebvre. *Fast unfolding of communities in large networks.* Journal of Statistical Mechanics: Theory and Experiment 2008, P10008 (2008). `doi:10.1088/1742-5468/2008/10/P10008 <https://doi.org/10.1088/1742-5468/2008/10/P10008>`_

* **CoreSPECT (cplearn)**: If you use cplearn clustering, please cite the `cplearn <https://github.com/csmukherjee/cplearn>`_ package and the CoreSPECT paper (see :doc:`../citations`).

Main Functions
--------------

Unified Clustering Function
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.clustering.cluster

   **Usage Examples:**

   .. code-block:: python

      from lotus.workflows import cluster
      
      # Use Leiden algorithm (default)
      cluster(adata, method="leiden", cluster_resolution=0.5)
      
      # Use Louvain algorithm
      cluster(adata, method="louvain", cluster_resolution=0.5)
      
      # View clustering results
      print(adata.obs["leiden"].value_counts())

Individual Clustering Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Leiden Clustering
^^^^^^^^^^^^^^^^^

.. autofunction:: lotus.workflows.clustering.leiden

   **Usage Example:**

   .. code-block:: python

      from lotus.workflows import leiden
      
      # Leiden clustering
      leiden(adata, resolution=0.5, key_added="leiden")
      
      # View clustering results
      print(adata.obs["leiden"].value_counts())

Louvain Clustering
^^^^^^^^^^^^^^^^^^

.. autofunction:: lotus.workflows.clustering.louvain

   **Usage Example:**

   .. code-block:: python

      from lotus.workflows import louvain
      
      # Louvain clustering
      louvain(adata, resolution=0.5, key_added="louvain")
      
      # View clustering results
      print(adata.obs["louvain"].value_counts())

Cplearn Clustering (CoreSPECT)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: lotus.workflows.clustering.cplearn_cluster

   **Usage Example:**

   .. code-block:: python

      from lotus.workflows import cplearn_cluster
      
      # Cplearn clustering (CoreSPECT)
      model = cplearn_cluster(
          adata,
          use_rep="X_pca",
          key_added="cplearn",
          cluster_resolution=1.2,
      )
      
      # View clustering results
      print(adata.obs["cplearn"].value_counts())
      
      # The model can be used for core analysis
      from lotus.workflows import core_analyze
      core_analyze(adata, model=model)

Helper Functions
----------------

.. autofunction:: lotus.workflows.clustering.summarize_clusters

   Summarize cluster label statistics.

