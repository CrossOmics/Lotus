Visualization
=============

Visualization is an important way to display analysis results and help understand the structure of the data and clustering results.

Overview
--------

Lotus provides two main visualization functions:
- **UMAP Visualization**: Display cell distribution in reduced-dimensional space
- **Marker Gene Visualization**: Display expression patterns of marker genes

UMAP Visualization
------------------

UMAP (Uniform Manifold Approximation and Projection) is a dimensionality reduction technique used to project high-dimensional data into 2D or 3D space for visualization.

Basic Usage
~~~~~~~~~~~

.. code-block:: python

    from lotus.workflows import umap
    
    umap(
        adata,
        cluster_key="cplearn_labels",  # Cluster key for coloring (can be "cplearn_labels", "leiden", or "louvain")
        truth_key="truth",            # Truth labels (optional)
        output_dir="./results",       # Output directory
        save="_clusters.png",         # Save filename
        show=False,                   # Whether to show
    )
    
    # Note: If you used scanpy for clustering, use the corresponding cluster key:
    # cluster_key="leiden"  # if using sc.tl.leiden()
    # cluster_key="louvain"  # if using sc.tl.louvain()

Parameter Description
~~~~~~~~~~~~~~~~~~~~

- ``cluster_key``: Cluster key for coloring, default ``"cplearn_labels"``

- ``truth_key``: Truth label key (if available), for comparison, default ``"truth"``

- ``output_dir``: Output directory, if ``None``, uses current figdir setting

- ``save``: Save filename, if ``False``, don't save

- ``show``: Whether to show the figure, default False

- ``compute_umap``: Whether to compute UMAP (if not already computed), default True

- ``min_dist``: UMAP minimum distance parameter, default 0.5

- ``spread``: UMAP spread parameter, default 1.0

- ``n_components``: Number of UMAP components, default 2

Auto-compute UMAP
~~~~~~~~~~~~~~~~~

If ``adata.obsm["X_umap"]`` doesn't exist, the function will automatically compute it:

.. code-block:: python

    # Auto-compute UMAP (if not already computed)
    umap(adata, compute_umap=True)

    # Use already computed UMAP
    umap(adata, compute_umap=False)

Marker Gene Visualization
-------------------------

Use ``render_visualizations()`` to generate both UMAP plots and marker gene violin plots:

.. code-block:: python

    from lotus.workflows import render_visualizations
    
    # Assume you have already found marker genes
    marker_genes_list = ["Gene1", "Gene2", "Gene3", "Gene4", "Gene5"]
    
    render_visualizations(
        adata,
        marker_genes=marker_genes_list,
        output_dir="./results",
        cluster_key="cplearn_labels",
        truth_key="truth",
    )

This generates:
- UMAP clustering plot (``umap_clusters.png``)
- Marker gene violin plot (``violin_markers.png``)

Using scanpy for More Visualizations
--------------------------------------

Since Lotus is fully compatible with scanpy, you can use all of scanpy's visualization functions:

.. code-block:: python

    import scanpy as sc
    
    # UMAP plot
    sc.pl.umap(adata, color="cplearn_labels")
    
    # Marker gene expression plot
    sc.pl.umap(adata, color=["Gene1", "Gene2", "Gene3"])
    
    # Violin plot
    sc.pl.violin(adata, keys=["Gene1"], groupby="cplearn_labels")
    
    # Heatmap
    sc.pl.heatmap(adata, var_names=marker_genes, groupby="cplearn_labels")

Adjusting UMAP Parameters
--------------------------

If UMAP plot quality is not ideal, you can adjust parameters:

.. code-block:: python

    # Tighter clusters (smaller min_dist)
    umap(adata, min_dist=0.1, spread=0.5)
    
    # More dispersed clusters (larger min_dist)
    umap(adata, min_dist=1.0, spread=2.0)

Next Steps
----------

After completing visualization, you can:
- Analyze clustering results
- Perform further biological interpretation
- Prepare publication-quality figures
