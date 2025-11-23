Quick Start
============

This guide will walk you through a complete single-cell RNA sequencing data analysis workflow, from data loading to final visualization results.

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

    print(f"Data shape: {adata.shape}")
    print(f"PCA stored in: adata.obsm['X_pca'] (shape: {adata.obsm['X_pca'].shape})")
    print(f"Latent representation: adata.obsm['X_latent'] (shape: {adata.obsm['X_latent'].shape})")
    print(f"Raw counts saved in: adata.layers['raw_counts']")
    print(f"Neighbors graph constructed: {'neighbors' in adata.uns}")

Example output:

.. code-block:: text

    Data shape: (180, 2000)
    PCA stored in: adata.obsm['X_pca'] (shape: (180, 20))
    Latent representation: adata.obsm['X_latent'] (shape: (180, 12))
    Raw counts saved in: adata.layers['raw_counts']
    Neighbors graph constructed: True

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

Example output:

.. code-block:: text

    Cluster summary: 0 (60), 1 (60), 2 (60)

You can check the clustering results:

.. code-block:: python

    print(f"Cluster labels stored in: adata.obs['cplearn_labels']")
    unique_labels = adata.obs["cplearn_labels"].unique()
    print(f"Number of clusters: {len(unique_labels)}")
    print(f"Cluster IDs: {sorted(unique_labels.tolist())}")

Example output:

.. code-block:: text

    Cluster labels stored in: adata.obs['cplearn_labels']
    Number of clusters: 3
    Cluster IDs: [0, 1, 2]

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

Example output:

.. code-block:: text

    Stored anchored map embedding in `adata.obsm['X_cplearn_coremap']` (180/180 points assigned).

You can check the embedding:

.. code-block:: python

    embedding = adata.obsm["X_cplearn_coremap"]
    assigned = np.sum(~np.isnan(embedding).any(axis=1))
    print(f"Core map embedding shape: {embedding.shape}")
    print(f"Assigned points: {assigned}/{adata.n_obs} ({100*assigned/adata.n_obs:.1f}%)")

Example output:

.. code-block:: text

    Core map embedding shape: (180, 2)
    Assigned points: 180/180 (100.0%)

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

    print(f"UMAP embedding stored in: adata.obsm['X_umap'] (shape: {adata.obsm['X_umap'].shape})")
    print(f"Visualization saved to: ./results/umap_clusters.png")

Example output:

.. code-block:: text

    UMAP embedding stored in: adata.obsm['X_umap'] (shape: (180, 2))
    Visualization saved to: ./results/umap_clusters.png

The output file contains a UMAP plot colored by cluster labels, showing the cell type separation.

The visualization shows:

- **UMAP embedding**: 2D representation of cells in the UMAP space
- **Color coding**: Each cluster is colored differently
- **Cell distribution**: Shows how cells are grouped and separated by cell type

Example visualization (saved as PNG):

.. code-block:: text

    The UMAP plot displays:
    - X-axis: UMAP dimension 1
    - Y-axis: UMAP dimension 2
    - Colors: Different clusters (e.g., cluster 0 in blue, cluster 1 in red, cluster 2 in green)
    - Each point represents a single cell
    - Well-separated clusters indicate distinct cell types

The generated PNG file (``./results/umap_clusters.png``) can be opened in any image viewer to see the visualization.

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

    print(f"Total differentially expressed genes: {len(de_result)}")
    print(f"Significant genes (p_adj < 0.05): {(de_result['p_adj'] < 0.05).sum()}")
    print(f"Significant genes (p_adj < 0.01): {(de_result['p_adj'] < 0.01).sum()}")
    print("\nTop 10 differentially expressed genes:")
    print(de_result[['gene', 'log2fc', 'p_adj', 'mean_a', 'mean_b']].head(10))

Example output:

.. code-block:: text

    Total differentially expressed genes: 50
    Significant genes (p_adj < 0.05): 45
    Significant genes (p_adj < 0.01): 32

    Top 10 differentially expressed genes:
          gene    log2fc      p_adj     mean_a     mean_b
    0  gene_023  2.345678  0.000123  15.234567   3.456789
    1  gene_012  2.123456  0.000234  12.345678   2.345678
    2  gene_045  1.987654  0.000456  10.987654   1.234567
    ...

The result is a pandas DataFrame containing gene names, log2 fold changes, p-values, and expression statistics.

Complete Example
----------------

Here's a complete workflow example that you can run:

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

Example output:

.. code-block:: text

    Preprocessing complete. Data shape: (180, 2000)
    Neighbors graph ready: True
    Cluster summary: 0 (60), 1 (60), 2 (60)
    Clustering complete. Found 3 clusters
    Stored anchored map embedding in `adata.obsm['X_cplearn_coremap']` (180/180 points assigned).
    Core selection complete
    UMAP visualization saved to ./results/umap_clusters.png
    DEG analysis complete. Found 50 differentially expressed genes
    Top 5 marker genes: ['gene_023', 'gene_012', 'gene_045', 'gene_067', 'gene_089']

    ✓ Analysis complete! Check ./results/ for output files.

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
