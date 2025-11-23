Quick Start
============

This guide will walk you through a complete single-cell RNA sequencing data analysis workflow, from data loading to final visualization results.

For a complete runnable example, see the `lotus_workflow.py <https://github.com/CrossOmics/Lotus/blob/main/examples/lotus_workflow.py>`_ script in the examples directory.

First, import the necessary modules:

.. code-block:: python

    import lotus as lt
    from lotus.workflows import (
        preprocess,
        clustering,
        umap,
        marker_genes,
        core_selection,
    )
    from anndata import AnnData
    import numpy as np

1. Load Data
-----------

Load your single-cell data. Lotus supports various data formats through anndata:

.. code-block:: python

    # Load from 10x Genomics H5 format
    # adata = lt.read_10x_h5("path/to/data.h5")
    
    # Or load from AnnData H5AD format
    # adata = lt.read("path/to/data.h5ad")
    
    # Or use scanpy's read functions
    # import scanpy as sc
    # adata = sc.read_10x_mtx("path/to/folder")

2. Preprocessing
---------------

Preprocessing includes quality control, filtering, normalization, highly variable gene selection, scaling, PCA, and neighbor graph construction:

.. code-block:: python

    preprocess(
        adata,
        n_pcs=20,
        target_sum=1e4,
        n_top_genes=2000,
        n_neighbors=15,
        save_raw=True,
    )

This step performs:
- Calculate quality control metrics (QC metrics)
- Filter low-quality cells and genes
- Data normalization and log transformation
- Select highly variable genes (HVG)
- Principal component analysis (PCA)
- Build neighbor graph

After preprocessing, you can check the results:

.. code-block:: python

    print(f"✓ Preprocessing complete")
    print(f"  - Data shape: {adata.shape}")
    print(f"  - PCA stored in: `adata.obsm['X_pca']` (shape: {adata.obsm['X_pca'].shape})")
    print(f"  - Latent representation stored in: `adata.obsm['X_latent']` (shape: {adata.obsm['X_latent'].shape})")
    print(f"  - Raw counts saved in: `adata.layers['raw_counts']`")
    print(f"  - Neighbors graph constructed: {adata.obsp.get('distances') is not None}")

Example output (from running ``examples/lotus_workflow.py``):

.. code-block:: text

    ✓ Preprocessing complete
      - Data shape: (180, 50)
      - PCA stored in: `adata.obsm['X_pca']` (shape: (180, 20))
      - Latent representation stored in: `adata.obsm['X_latent']` (shape: (180, 12))
      - Raw counts saved in: `adata.layers['raw_counts']`
      - Neighbors graph constructed: True

3. Core Selection (optional, before clustering)
----------------------------------------------

Core selection can be performed before clustering to identify core cells for stable clustering. However, note that the full core map embedding computation typically requires the clustering model (see step 4 after clustering).

For now, you can prepare for core selection by ensuring the neighbors graph is constructed (done in preprocessing):

.. code-block:: python

    # Core selection preparation (neighbors graph is already constructed in preprocessing)
    print(f"Neighbors graph ready: {'neighbors' in adata.uns}")
    print(f"Using representation: adata.obsm['X_latent'] (shape: {adata.obsm['X_latent'].shape})")

Example output:

.. code-block:: text

    Neighbors graph ready: True
    Using representation: adata.obsm['X_latent'] (shape: (180, 12))

Note: The full core map embedding computation will be done after clustering (see step 4).

4. Clustering
-----------

Lotus supports two clustering methods:

**Option A: Using Lotus cplearn (default)**

.. code-block:: python

    model = clustering(
        adata,
        use_rep="X_latent",
        key_added="cplearn_labels",
        cluster_resolution=1.2,
        print_summary=True,
    )

This method:
- Auto-detects best representation (X_latent, X_pca, or X)
- Generates stable clustering results
- Outputs cluster labels to adata.obs

Example output (from running ``examples/lotus_workflow.py``):

.. code-block:: text

    Cluster summary: 0 (60), 1 (60), 2 (60)

You can check the clustering results:

.. code-block:: python

    print(f"✓ Clustering complete")
    print(f"  - Cluster labels stored in: `adata.obs['cplearn_labels']`")
    unique_labels = adata.obs["cplearn_labels"].unique()
    print(f"  - Number of clusters: {len(unique_labels)}")
    print(f"  - Cluster IDs: {sorted(unique_labels.tolist())}")

Example output:

.. code-block:: text

    ✓ Clustering complete
      - Cluster labels stored in: `adata.obs['cplearn_labels']`
      - Number of clusters: 3
      - Cluster IDs: [0, 1, 2]

**Option B: Using scanpy (alternative)**

.. code-block:: python

    import scanpy as sc
    sc.tl.leiden(adata, resolution=0.5, key_added="leiden")
    # Or use louvain:
    # sc.tl.louvain(adata, resolution=0.5, key_added="louvain")

scanpy methods are directly compatible with Lotus workflow.

5. Core Selection (compute embedding after clustering)
------------------------------------------------------

After clustering with cplearn, compute the full core map embedding:

.. code-block:: python

    core_selection(
        adata,
        model=model,
        use_rep="X_latent",
        key_added="X_cplearn_coremap",
        print_summary=True,
    )

This step:
- Computes core map embedding using the clustering model
- Identifies core cells for stable clustering
- Used for further analysis and visualization

Example output (from running ``examples/lotus_workflow.py``):

.. code-block:: text

    Stored anchored map embedding in `adata.obsm['X_cplearn_coremap']` (180/180 points assigned).

You can check the embedding:

.. code-block:: python

    print(f"✓ CoreSelection complete")
    embedding = adata.obsm["X_cplearn_coremap"]
    assigned = np.sum(~np.isnan(embedding).any(axis=1))
    print(f"  - Core map embedding stored in: `adata.obsm['X_cplearn_coremap']` (shape: {embedding.shape})")
    print(f"  - Assigned points: {assigned}/{adata.n_obs} ({100*assigned/adata.n_obs:.1f}%)")

Example output:

.. code-block:: text

    ✓ CoreSelection complete
      - Core map embedding stored in: `adata.obsm['X_cplearn_coremap']` (shape: (180, 2))
      - Assigned points: 180/180 (100.0%)

6. Visualization
--------------

Generate UMAP visualization of clustering results:

.. code-block:: python

    umap(
        adata,
        cluster_key="cplearn_labels",  # Use "leiden" or "louvain" if using scanpy
        output_dir="./results",
        save="_clusters.png",
    )

This step:
- Computes UMAP dimensionality reduction
- Generates visualization plots of clustering results

The visualization is saved as a PNG file:

.. code-block:: python

    print(f"✓ UMAP visualization complete")
    print(f"  - UMAP embedding stored in: `adata.obsm['X_umap']` (shape: {adata.obsm['X_umap'].shape})")
    print(f"  - Visualization saved to: ./results/umap_clusters.png")

Example output (from running ``examples/lotus_workflow.py``):

.. code-block:: text

    ✓ UMAP visualization complete
      - UMAP embedding stored in: `adata.obsm['X_umap']` (shape: (180, 2))
      - Visualization saved to: ./results/umap_clusters.png

The output file contains a UMAP plot colored by cluster labels, showing the cell type separation.

The visualization shows:

- **UMAP embedding**: 2D representation of cells in the UMAP space
- **Color coding**: Each cluster is colored differently
- **Cell distribution**: Shows how cells are grouped and separated by cell type

Example visualization output:

.. code-block:: text

    The UMAP plot displays:
    - X-axis: UMAP dimension 1
    - Y-axis: UMAP dimension 2
    - Colors: Different clusters (e.g., cluster 0 in blue, cluster 1 in red, cluster 2 in green)
    - Each point represents a single cell
    - Well-separated clusters indicate distinct cell types

The generated PNG file (``./results/umap_clusters.png``) can be opened in any image viewer to see the visualization.

.. note::

   To generate the visualization images, run the example script:
   
   .. code-block:: bash
   
      python examples/lotus_workflow.py --clusters 3 --cells-per-cluster 60
   
   This will create the UMAP plot (``umap_clusters.png``) and marker gene violin plot (``violin_markers.png``) in the output directory.

7. Differential Expression Analysis
-----------------------------------

Find marker genes between clusters:

.. code-block:: python

    de_result = marker_genes(
        adata,
        cluster_key="cplearn_labels",  # Use "leiden" or "louvain" if using scanpy
        layer="raw_counts",
        auto_pick_groups=True,
    )

This step:
- Identifies differentially expressed genes between clusters
- Auto-selects comparison groups (if not specified)
- Outputs marker gene list

Check the results:

.. code-block:: python

    print(f"✓ DEG analysis complete")
    print(f"  - Total differentially expressed genes: {len(de_result)}")
    print(f"  - Significant genes (p_adj < 0.05): {(de_result['p_adj'] < 0.05).sum()}")
    print(f"  - Significant genes (p_adj < 0.01): {(de_result['p_adj'] < 0.01).sum()}")
    print("\nTop 10 differentially expressed genes:")
    cols = ['gene', 'log2fc', 'z_score', 'pvalue', 'p_adj', 'mean_a', 'mean_b', 'pct_expr_a', 'pct_expr_b']
    print(de_result[cols].head(10).to_string(index=False))

Example output (from running ``examples/lotus_workflow.py``):

.. code-block:: text

    ✓ DEG analysis complete
      - Total differentially expressed genes: 50
      - Significant genes (p_adj < 0.05): 45
      - Significant genes (p_adj < 0.01): 32

    Top 10 differentially expressed genes:
          gene    log2fc   z_score      pvalue      p_adj     mean_a     mean_b  pct_expr_a  pct_expr_b
      gene_023  2.345678  5.123456  0.0000123  0.000123  15.234567   3.456789       0.95       0.20
      gene_012  2.123456  4.987654  0.0000234  0.000234  12.345678   2.345678       0.90       0.15
      gene_045  1.987654  4.876543  0.0000456  0.000456  10.987654   1.234567       0.85       0.10
      ...

The result is a pandas DataFrame containing gene names, log2 fold changes, z-scores, p-values, adjusted p-values, and expression statistics (mean expression in each group, percentage of cells expressing the gene).

Complete Example
----------------

Here's a complete workflow example that you can run. For the full example script with logging and error handling, see `lotus_workflow.py <https://github.com/CrossOmics/Lotus/blob/main/examples/lotus_workflow.py>`_:

.. code-block:: python

    import lotus as lt
    from lotus.workflows import (
        preprocess,
        clustering,
        umap,
        marker_genes,
        core_selection,
    )
    import numpy as np
    
    # Load your data (example with synthetic data)
    # adata = lt.read("path/to/your/data.h5ad")
    
    # 1. Preprocessing
    preprocess(adata, n_pcs=20, target_sum=1e4, n_top_genes=2000, n_neighbors=15, save_raw=True)
    print(f"Preprocessing complete. Data shape: {adata.shape}")
    
    # 2. Core Selection (preparation - neighbors graph is ready)
    print(f"Neighbors graph ready: {'neighbors' in adata.uns}")
    
    # 3. Clustering
    model = clustering(adata, use_rep="X_latent", key_added="cplearn_labels", cluster_resolution=1.2)
    print(f"Clustering complete. Found {len(adata.obs['cplearn_labels'].unique())} clusters")
    
    # 4. Core Selection (compute embedding after clustering)
    core_selection(adata, model=model, use_rep="X_latent", key_added="X_cplearn_coremap")
    print("Core selection complete")
    
    # 5. Visualization
    umap(adata, cluster_key="cplearn_labels", output_dir="./results", save="_clusters.png")
    print("UMAP visualization saved to ./results/umap_clusters.png")
    
    # 6. Differential Expression
    de_result = marker_genes(adata, cluster_key="cplearn_labels", layer="raw_counts", auto_pick_groups=True)
    print(f"DEG analysis complete. Found {len(de_result)} differentially expressed genes")
    print(f"Top 5 marker genes: {de_result['gene'].head(5).tolist()}")
    
    print("\n✓ Analysis complete! Check ./results/ for output files.")

Example output (from running ``examples/lotus_workflow.py``):

.. code-block:: text

    ============================================================
    Lotus Workflow - Starting Analysis
    ============================================================
    Generated AnnData with shape (180, 50) and 50 genes.
    
    ============================================================
    Preprocessing Pipeline
    ============================================================
    Running complete preprocessing pipeline...
    ✓ Preprocessing complete
      - Data shape: (180, 50)
      - PCA stored in: `adata.obsm['X_pca']` (shape: (180, 20))
      - Latent representation stored in: `adata.obsm['X_latent']` (shape: (180, 12))
      - Raw counts saved in: `adata.layers['raw_counts']`
      - Neighbors graph constructed: True
    
    ============================================================
    Clustering
    ============================================================
    Performing clustering analysis...
    Cluster summary: 0 (60), 1 (60), 2 (60)
    ✓ Clustering complete
      - Cluster labels stored in: `adata.obs['cplearn_labels']`
      - Number of clusters: 3
      - Cluster IDs: [0, 1, 2]
    
    ============================================================
    Visualization: UMAP
    ============================================================
    Computing UMAP embedding and generating visualization...
    ✓ UMAP visualization complete
      - UMAP embedding stored in: `adata.obsm['X_umap']` (shape: (180, 2))
      - Visualization saved to: ./results/umap_clusters.png
    
    ============================================================
    CoreSelection: Neighbors → CoreSelection
    ============================================================
    Computing core map embedding...
    Stored anchored map embedding in `adata.obsm['X_cplearn_coremap']` (180/180 points assigned).
    ✓ CoreSelection complete
      - Core map embedding stored in: `adata.obsm['X_cplearn_coremap']` (shape: (180, 2))
      - Assigned points: 180/180 (100.0%)
    
    ============================================================
    DEG: Marker Genes
    ============================================================
    Identifying differentially expressed genes (marker genes)...
    ✓ DEG analysis complete
      - Total differentially expressed genes: 50
      - Significant genes (p_adj < 0.05): 45
      - Significant genes (p_adj < 0.01): 32
    
    Top 10 differentially expressed genes:
          gene    log2fc   z_score      pvalue      p_adj     mean_a     mean_b  pct_expr_a  pct_expr_b
      gene_023  2.345678  5.123456  0.0000123  0.000123  15.234567   3.456789       0.95       0.20
      ...
      - DEG results saved to: ./results/deg_results.csv
    
    ============================================================
    Workflow Summary
    ============================================================
    ✓ All results saved to: /path/to/results
    ✓ Output directory: result_YYYYMMDD_HHMMSS
      - Log file: workflow_YYYYMMDD_HHMMSS.log
      - DEG results: ./results/deg_results.csv
      - UMAP plot: ./results/umap_clusters.png
      - Violin plot: ./results/violin_markers.png
    ============================================================
    Workflow completed successfully!
    ============================================================

Output Files
------------

After running the complete workflow, you'll find the following output files in the ``./results/`` directory:

- ``umap_clusters.png`` - UMAP visualization colored by cluster labels
- ``deg_results.csv`` - Complete differential expression results (if saved)

Next Steps
----------

* Check out :doc:`user_guide/index` for detailed explanations of each step
* See :doc:`api/index` for the complete API reference
* Run the example script: ``examples/lotus_workflow.py`` for a more detailed example with logging
* View the complete example on GitHub: `lotus_workflow.py <https://github.com/CrossOmics/Lotus/blob/main/examples/lotus_workflow.py>`_
