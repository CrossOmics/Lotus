Core Selection Module
======================

The ``lotus.workflows.core_selection`` module provides core map embedding functionality.

Overview
--------

Core selection is used to compute core map embedding, a special dimensionality reduction representation that can be used for further analysis and visualization.

Main Functions
--------------

Core Selection
~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.core_selection.core_selection.core_selection

   **Biological Background:**
   
   Core map embedding:
   - A special embedding computed based on clustering results
   - Can highlight core features of cells
   - Used for further analysis and visualization
   
   Note: This function needs to be used after clustering analysis and requires a neighbor graph.

   **Usage Example:**

   .. code-block:: python

      from lotus.workflows import clustering, core_selection
      
      # First perform clustering
      model = clustering(adata, key_added="cplearn_labels")
      
      # Then compute core map embedding
      core_selection(
          adata,
          model=model,
          use_rep="X_latent",
          key_added="X_cplearn_coremap",
      )
      
      # View results
      print(adata.obsm["X_cplearn_coremap"].shape)

.. autofunction:: lotus.workflows.core_selection.core_selection.compute_coremap_embedding

   This is an alias for ``core_selection()`` for backward compatibility.
