Visualization Module
====================

The ``lotus.workflows.visualization`` module provides result visualization functionality.

Overview
--------

Visualization is an important way to display analysis results and help understand the structure of the data and clustering results.

Main Functions
--------------

UMAP Visualization
~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.visualization.visualization.umap

   **Biological Background:**
   
   UMAP (Uniform Manifold Approximation and Projection) projects high-dimensional data into 2D or 3D space while preserving both local and global structure of the data. It is used to visualize cell distribution in expression space. UMAP plots can help understand relationships between cell types, identify cell subpopulations, and validate clustering results.

   **Usage Example:**

   .. code-block:: python

      from lotus.workflows import umap
      
      # Basic usage (with Lotus cplearn clustering)
      umap(
          adata,
          cluster_key="cplearn_labels",  # Can also use "leiden" or "louvain" if using scanpy
          truth_key="truth",
          output_dir="./results",
          save="_clusters.png",
      )
      
      # If using scanpy clustering, use the corresponding cluster key:
      # cluster_key="leiden"  # if using sc.tl.leiden()
      # cluster_key="louvain"  # if using sc.tl.louvain()

Marker Gene Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.visualization.visualization.render_visualizations

   **Biological Background:**
   
   Marker gene visualization includes UMAP clustering plots that show cell distribution in reduced-dimensional space, and violin plots that show expression distribution of marker genes across different clusters. These visualizations can help validate marker gene specificity, understand characteristics of different cell types, and prepare publication-quality figures.

   **Usage Example:**

   .. code-block:: python

      from lotus.workflows import render_visualizations
      
      # Assume you have already found marker genes
      marker_genes = ["CD3D", "CD79A", "MS4A1", "CD14", "FCGR3A"]
      
      render_visualizations(
          adata,
          marker_genes=marker_genes,
          output_dir="./results",
          cluster_key="cplearn_labels",  # Can also use "leiden" or "louvain" if using scanpy
      )

CoreMap Visualization
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.visualization.visualization.coremap

   **Biological Background:**
   
   CoreMap visualization visualizes the core map embedding computed by core_analysis() and provides interactive Plotly plots with layer sliders. It highlights core cell states and cell relationships, using cplearn's own visualization function for optimal display. This visualization is particularly useful for understanding core-periphery structure in cell populations, exploring cell state transitions, and interactive exploration of core map embeddings.
   
   **Usage Example:**
   
   .. code-block:: python
   
      from lotus.workflows import coremap, core_analysis
      from lotus.methods.cplearn.external import cplearn
      
      # First perform core analysis (before clustering)
      model = cplearn.corespect(adata, use_rep="X_pca", key_added="cplearn")
      core_analysis(adata, model=model, key_added="X_cplearn_coremap")
      
      # Then perform clustering using the core analysis results
      
      # Then visualize the core map
      coremap(
          adata,
          coremap_key="X_cplearn_coremap",
          cluster_key="cplearn",
          model=model,
          output_dir="./results",
          save="_coremap.html",
      )
   
   **Note:** The output is an interactive HTML file (Plotly format) that can be opened in a web browser. The visualization includes layer sliders to explore different core layers.
