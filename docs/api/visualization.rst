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
   
   UMAP (Uniform Manifold Approximation and Projection):
   - Projects high-dimensional data into 2D or 3D space
   - Preserves local and global structure of the data
   - Used to visualize cell distribution in expression space
   
   UMAP plots can help:
   - Understand relationships between cell types
   - Identify cell subpopulations
   - Validate clustering results

   **Usage Example:**

   .. code-block:: python

      from lotus.workflows import umap
      
      # Basic usage
      umap(
          adata,
          cluster_key="cplearn_labels",
          truth_key="truth",
          output_dir="./results",
          save="_clusters.png",
      )

Marker Gene Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.visualization.visualization.render_visualizations

   **Biological Background:**
   
   Marker gene visualization includes:
   - UMAP clustering plot: Shows cell distribution in reduced-dimensional space
   - Violin plot: Shows expression distribution of marker genes across different clusters
   
   These visualizations can help:
   - Validate marker gene specificity
   - Understand characteristics of different cell types
   - Prepare publication-quality figures

   **Usage Example:**

   .. code-block:: python

      from lotus.workflows import render_visualizations
      
      # Assume you have already found marker genes
      marker_genes = ["CD3D", "CD79A", "MS4A1", "CD14", "FCGR3A"]
      
      render_visualizations(
          adata,
          marker_genes=marker_genes,
          output_dir="./results",
          cluster_key="cplearn_labels",
      )
