Visualization Module
====================

The ``lotus.workflows.visualization`` module provides result visualization functionality.

Overview
--------

Visualization is an important way to display analysis results and help understand the structure of the data and clustering results.

**Citation:**

* **UMAP**: L. McInnes, J. Healy, and J. Melville. *UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction.* `arXiv:1802.03426 <https://arxiv.org/abs/1802.03426>`_

* **t-SNE**: L. van der Maaten and G. Hinton. *Visualizing Data using t-SNE.* Journal of Machine Learning Research 9, 2579-2605 (2008).

* **Diffusion Maps**: R. R. Coifman and S. Lafon. *Diffusion maps.* Applied and Computational Harmonic Analysis 21, 5-30 (2006).

* **CoreMap**: If you use CoreMap visualization, please cite the `cplearn <https://github.com/csmukherjee/cplearn>`_ package and the CoreSPECT paper (see :doc:`../citations`).

Main Functions
--------------

UMAP Visualization
~~~~~~~~~~~~~~~~~~~
Unified Visualization Function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.visualization.visualization.visualization

   **Usage Example:**

   .. code-block:: python

      from lotus.workflows import visualization
      
      # Use different visualization methods
      visualization(adata, method="umap", cluster_key="leiden")
      visualization(adata, method="tsne", perplexity=50)
      visualization(adata, method="diffmap", n_comps=20)
      visualization(adata, method="draw_graph", layout="fa")
      visualization(adata, method="coremap", model=model)


.. autofunction:: lotus.workflows.visualization.visualization.umap

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


t-SNE Visualization
~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.visualization.visualization.tsne

   **Usage Example:**

   .. code-block:: python

      from lotus.workflows import tsne
      
      tsne(
          adata,
          cluster_key="leiden",
          truth_key="truth",
          output_dir="./results",
          save="_tsne.png",
          perplexity=50,
      )

Diffusion Map Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.visualization.visualization.diffmap

   **Usage Example:**

   .. code-block:: python

      from lotus.workflows import diffmap
      
      diffmap(
          adata,
          cluster_key="leiden",
          truth_key="truth",
          output_dir="./results",
          save="_diffmap.png",
          n_comps=20,
      )

Force-Directed Graph Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.visualization.visualization.draw_graph

   **Usage Example:**

   .. code-block:: python

      from lotus.workflows import draw_graph
      
      draw_graph(
          adata,
          cluster_key="leiden",
          truth_key="truth",
          output_dir="./results",
          save="_draw_graph.png",
          layout="fa",
      )

Marker Gene Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.visualization.visualization.render_visualizations

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
