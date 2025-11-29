Core Analysis Module
=====================

The ``lotus.workflows.core_analysis`` module provides core map embedding functionality.

Overview
--------

Core analysis is used to compute core map embedding, a special dimensionality reduction representation that can be used for further analysis and visualization.

**Note on workflow order:** Core analysis is performed before clustering to identify core cells and prepare for stable clustering.

**Citation:** If you use core analysis features, please cite the `cplearn <https://github.com/csmukherjee/cplearn>`_ package and the following papers:

* **CoreSPECT**: Chandra Sekhar Mukherjee, Joonyoung Bae, and Jiapeng Zhang. *CoreSPECT: Enhancing Clustering Algorithms via an Interplay of Density and Geometry.* `arXiv:2507.08243 <https://arxiv.org/abs/2507.08243>`_

* **Balanced Ranking**: Chandra Sekhar Mukherjee and Jiapeng Zhang. *Balanced Ranking with Relative Centrality: A Multi-Core Periphery Perspective.* ICLR 2025. 
Main Functions
--------------

Core Analysis
~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.core_analysis.core_analyze

   **Usage Example:**

   .. code-block:: python

      from lotus.workflows import core_analyze
      from lotus.methods.cplearn.external import cplearn
      
      # First perform clustering to get a model
      model = cplearn.corespect(adata, use_rep="X_pca", key_added="cplearn")
      
      # Then compute core map embedding (after clustering)
      # Note: This requires a neighbors graph from preprocessing
      core_analyze(
          adata,
          model=model,
          use_rep="X_pca",
          key_added="X_cplearn_coremap",
      )
      
      # View results
      print(adata.obsm["X_cplearn_coremap"].shape)
      print(adata.obs["X_cplearn_coremap_is_core"].head())
