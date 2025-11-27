Core Analysis Module
=====================

The ``lotus.workflows.core_analysis`` module provides core map embedding functionality.

Overview
--------

Core analysis is used to compute core map embedding, a special dimensionality reduction representation that can be used for further analysis and visualization.

**Note on workflow order:** Core analysis is performed before clustering to identify core cells and prepare for stable clustering. It is compatible with both cplearn and scanpy clustering workflows.
Main Functions
--------------

Core Analysis
~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.core_analysis.core_analysis.core_analysis

   **Biological Background:**
   
   Core map embedding is a special embedding computed before clustering to identify core cells. It can highlight core features of cells, is used to prepare data for stable clustering, and helps identify stable cell populations. This function should be used before clustering analysis and requires a neighbor graph from preprocessing.

   **Usage Example:**

   .. code-block:: python

      from lotus.workflows import core_analysis
      from lotus.methods.cplearn.external import cplearn
      
      # First compute core map embedding (before clustering)
      # Note: This requires a neighbors graph from preprocessing
      model = cplearn.corespect(adata, use_rep="X_latent", key_added="cplearn")
      core_analysis(
          adata,
          model=model,
          use_rep="X_latent",
          key_added="X_cplearn_coremap",
      )
      
      # Then perform clustering using the core analysis results
      # (clustering can use the core map embedding for better results)
      
      # View results
      print(adata.obsm["X_cplearn_coremap"].shape)

.. autofunction:: lotus.workflows.core_analysis.core_analysis.compute_coremap_embedding

   This is an alias for ``core_analysis()`` for backward compatibility.
