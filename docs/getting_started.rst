Quick Start
============

This guide will walk you through a complete single-cell RNA sequencing data analysis workflow, from data loading to final visualization results.

All example outputs shown below are **real outputs** from running ``examples/lotus_workflow.py`` with default parameters (``--clusters 3 --cells-per-cluster 60``).

For a complete runnable example, see the `lotus_workflow.py <https://github.com/CrossOmics/Lotus/blob/main/examples/lotus_workflow.py>`_ script in the examples directory.

First, import the necessary modules:

.. code-block:: python

    import lotus as lt
    from lotus.workflows import (
        preprocess,
        clustering,
        umap,
        marker_genes,
        core_analysis,
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

    2025-11-23 00:13:39,358 - INFO - ✓ Preprocessing complete
    2025-11-23 00:13:39,358 - INFO -   - Data shape: (180, 50)
    2025-11-23 00:13:39,358 - INFO -   - PCA stored in: `adata.obsm['X_pca']` (shape: (180, 20))
    2025-11-23 00:13:39,358 - INFO -   - Latent representation stored in: `adata.obsm['X_latent']` (shape: (180, 32))
    2025-11-23 00:13:39,358 - INFO -   - Raw counts saved in: `adata.layers['raw_counts']`
    2025-11-23 00:13:39,358 - INFO -   - Neighbors graph constructed: True

3. Core Analysis (after clustering)
------------------------------------

Core analysis is performed after clustering to compute core map embedding. The neighbors graph is already constructed in the preprocessing step:

.. code-block:: python

    # Core analysis preparation (neighbors graph is already constructed in preprocessing)
    print(f"Neighbors graph ready: {'neighbors' in adata.uns}")
    print(f"Using representation: adata.obsm['X_latent'] (shape: {adata.obsm['X_latent'].shape})")

Example output:

.. code-block:: text

    Neighbors graph ready: True
    Using representation: adata.obsm['X_latent'] (shape: (180, 32))

This step prepares the data for stable clustering by ensuring the neighbors graph is ready. The actual core map embedding computation will be performed after clustering (as shown in the complete workflow output below).

4. Clustering
-----------

Lotus supports multiple clustering methods through a unified interface:

**Option A: Using Lotus cplearn (default)**

.. code-block:: python

    model = clustering(
        adata,
        method="cplearn",  # or omit for default
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

    Initiating FlowRank with r= 10
    Original average degree: 20.0
    Cluster summary: 0: 60, 1: 60, 2: 60
    2025-11-23 00:13:47,577 - INFO - ✓ Clustering complete
    2025-11-23 00:13:47,577 - INFO -   - Cluster labels stored in: `adata.obs['cplearn_labels']`
    2025-11-23 00:13:47,577 - INFO -   - Number of clusters: 3
    2025-11-23 00:13:47,577 - INFO -   - Cluster IDs: [0, 1, 2]

**Option B: Using scanpy Leiden algorithm**

.. code-block:: python

    clustering(
        adata,
        method="leiden",
        cluster_resolution=0.5,
        key_added="leiden",  # auto-set if None
    )

**Option C: Using scanpy Louvain algorithm**

.. code-block:: python

    clustering(
        adata,
        method="louvain",
        cluster_resolution=0.5,
        key_added="louvain",  # auto-set if None
    )

All methods output cluster labels in scanpy-compatible format and can be used with subsequent Lotus functions.

5. Visualization
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

    WARNING: saving figure to file figures/umap_clusters.png
    2025-11-23 00:13:48,497 - INFO - ✓ UMAP visualization complete
    2025-11-23 00:13:48,497 - INFO -   - UMAP embedding stored in: `adata.obsm['X_umap']` (shape: (180, 2))
    2025-11-23 00:13:48,497 - INFO -   - Visualization saved to: /tmp/lotus_real_output/umap_clusters.png

The output file contains a UMAP plot colored by cluster labels, showing the cell type separation.

The visualization shows:

- **UMAP embedding**: 2D representation of cells in the UMAP space
- **Color coding**: Each cluster is colored differently
- **Cell distribution**: Shows how cells are grouped and separated by cell type

Example UMAP visualization:

.. figure:: _static/umap_clusters_example.png
   :alt: UMAP visualization colored by cluster labels
   :width: 600px
   :align: center

   UMAP plot showing cell clusters. Each color represents a different cell type/cluster identified by the clustering algorithm. The plot displays cells in 2D UMAP space, with well-separated clusters indicating distinct cell types.
   
   *To generate this image, run: ``python examples/lotus_workflow.py --clusters 3 --cells-per-cluster 60``, then copy ``result_*/umap_clusters.png`` to ``docs/_static/umap_clusters_example.png``*

The UMAP plot shows:
- **X-axis**: UMAP dimension 1
- **Y-axis**: UMAP dimension 2  
- **Colors**: Different clusters (e.g., cluster 0 in blue, cluster 1 in red, cluster 2 in green)
- **Points**: Each point represents a single cell
- **Separation**: Well-separated clusters indicate distinct cell types

.. note::

   To generate the visualization images yourself, run the example script:
   
   .. code-block:: bash
   
      python examples/lotus_workflow.py --clusters 3 --cells-per-cluster 60
   
   This will create the UMAP plot (``umap_clusters.png``) and marker gene violin plot (``violin_markers.png``) in the output directory.
   
   The generated PNG file (``./results/umap_clusters.png``) can be opened in any image viewer to see the visualization.

6. Differential Expression Analysis
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

    Comparing groups {0} vs {1}
    2025-11-23 00:13:51,088 - INFO - ✓ DEG analysis complete
    2025-11-23 00:13:51,088 - INFO -   - Total differentially expressed genes: 50
    2025-11-23 00:13:51,088 - INFO -   - Significant genes (p_adj < 0.05): 41
    2025-11-23 00:13:51,088 - INFO -   - Significant genes (p_adj < 0.01): 39
    2025-11-23 00:13:51,088 - INFO - 
    Top 10 differentially expressed genes:
    2025-11-23 00:13:51,091 - INFO - 
        gene   log2fc   z_score       pvalue        p_adj   mean_a   mean_b  pct_expr_a  pct_expr_b
    gene_002 2.187279 18.400315 2.021035e-20 3.368391e-19 6.666667 0.683333    1.000000    0.500000
    gene_027 2.102686 19.505699 4.567841e-21 1.170854e-19 7.733333 1.033333    1.000000    0.600000
    gene_021 2.078502 16.554285 4.683415e-21 1.170854e-19 9.700000 1.533333    1.000000    0.733333
    gene_029 1.939962 12.955106 7.978928e-18 4.432738e-17 4.883333 0.533333    0.966667    0.366667
    gene_032 1.911840 13.144602 6.643614e-18 4.152259e-17 5.083333 0.616667    0.966667    0.433333
    gene_031 1.902821 15.913483 3.142358e-20 3.927948e-19 7.850000 1.366667    1.000000    0.750000
    gene_045 1.866375 14.216678 6.965233e-19 4.975166e-18 6.900000 1.166667    1.000000    0.716667
    gene_004 1.847075 14.939477 2.073976e-19 1.728313e-18 4.216667 0.450000    0.983333    0.350000
    gene_023 1.603167 10.278753 1.255096e-15 4.827291e-15 4.316667 0.750000    0.983333    0.500000
    gene_044 1.584963 12.506953 1.768928e-17 8.040582e-17 5.950000 1.316667    1.000000    0.666667
    2025-11-23 00:13:51,093 - INFO -   - DEG results saved to: /tmp/lotus_real_output/deg_results.csv

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
        core_analysis,
    )
    import numpy as np
    
    # Load your data (example with synthetic data)
    # adata = lt.read("path/to/your/data.h5ad")
    
    # 1. Preprocessing
    preprocess(adata, n_pcs=20, target_sum=1e4, n_top_genes=2000, n_neighbors=15, save_raw=True)
    print(f"Preprocessing complete. Data shape: {adata.shape}")
    
    # 2. Core Analysis (after clustering)
    print(f"Neighbors graph ready: {'neighbors' in adata.uns}")
    print("Core analysis preparation complete")
    
    # 3. Clustering
    # Using cplearn (default)
    model = clustering(adata, method="cplearn", use_rep="X_latent", key_added="cplearn_labels", cluster_resolution=1.2)
    # Or use scanpy methods:
    # clustering(adata, method="leiden", cluster_resolution=0.5)
    # clustering(adata, method="louvain", cluster_resolution=0.5)
    print(f"Clustering complete. Found {len(adata.obs['cplearn_labels'].unique())} clusters")
    
    # 4. Visualization
    umap(adata, cluster_key="cplearn_labels", output_dir="./results", save="_clusters.png")
    print("UMAP visualization saved to ./results/umap_clusters.png")
    
    # 5. Differential Expression
    de_result = marker_genes(adata, cluster_key="cplearn_labels", layer="raw_counts", auto_pick_groups=True)
    print(f"DEG analysis complete. Found {len(de_result)} differentially expressed genes")
    print(f"Top 5 marker genes: {de_result['gene'].head(5).tolist()}")
    
    print("\n✓ Analysis complete! Check ./results/ for output files.")

Example output (from running ``examples/lotus_workflow.py``):

.. code-block:: text

    2025-11-23 00:13:35,955 - INFO - ============================================================
    2025-11-23 00:13:35,955 - INFO - Lotus Workflow - Starting Analysis
    2025-11-23 00:13:35,955 - INFO - ============================================================
    2025-11-23 00:13:35,957 - INFO - Generated AnnData with shape (180, 50) and 50 genes.
    2025-11-23 00:13:35,957 - INFO - 
    ============================================================
    2025-11-23 00:13:35,957 - INFO - Preprocessing Pipeline
    2025-11-23 00:13:35,957 - INFO - ============================================================
    2025-11-23 00:13:35,958 - INFO - Running complete preprocessing pipeline...
    2025-11-23 00:13:39,358 - INFO - ✓ Preprocessing complete
    2025-11-23 00:13:39,358 - INFO -   - Data shape: (180, 50)
    2025-11-23 00:13:39,358 - INFO -   - PCA stored in: `adata.obsm['X_pca']` (shape: (180, 20))
    2025-11-23 00:13:39,358 - INFO -   - Latent representation stored in: `adata.obsm['X_latent']` (shape: (180, 32))
    2025-11-23 00:13:39,358 - INFO -   - Raw counts saved in: `adata.layers['raw_counts']`
    2025-11-23 00:13:39,358 - INFO -   - Neighbors graph constructed: True
    2025-11-23 00:13:39,358 - INFO - 
    ============================================================
    2025-11-23 00:13:39,358 - INFO - Clustering
    2025-11-23 00:13:39,358 - INFO - ============================================================
    2025-11-23 00:13:39,358 - INFO - Performing clustering analysis...
    Initiating FlowRank with r= 10
    Original average degree: 20.0
    Cluster summary: 0: 60, 1: 60, 2: 60
    2025-11-23 00:13:47,577 - INFO - ✓ Clustering complete
    2025-11-23 00:13:47,577 - INFO -   - Cluster labels stored in: `adata.obs['cplearn_labels']`
    2025-11-23 00:13:47,577 - INFO -   - Number of clusters: 3
    2025-11-23 00:13:47,577 - INFO -   - Cluster IDs: [0, 1, 2]
    2025-11-23 00:13:47,577 - INFO - 
    ============================================================
    2025-11-23 00:13:47,577 - INFO - Visualization: UMAP
    2025-11-23 00:13:47,577 - INFO - ============================================================
    2025-11-23 00:13:47,577 - INFO - Computing UMAP embedding and generating visualization...
    WARNING: saving figure to file figures/umap_clusters.png
    2025-11-23 00:13:48,497 - INFO - ✓ UMAP visualization complete
    2025-11-23 00:13:48,497 - INFO -   - UMAP embedding stored in: `adata.obsm['X_umap']` (shape: (180, 2))
    2025-11-23 00:13:48,497 - INFO -   - Visualization saved to: /tmp/lotus_real_output/umap_clusters.png
    2025-11-23 00:13:48,497 - INFO - 
    ============================================================
    2025-11-23 00:13:48,497 - INFO - CoreAnalysis: Neighbors → CoreAnalysis
    2025-11-23 00:13:48,497 - INFO - ============================================================
    2025-11-23 00:13:48,497 - INFO - Computing core map embedding...
    Total number of clusters= 3
    [60, 52, 48]
    Fitting GMM anchors: 100%|██████████| 3/3 [00:00<00:00, 10.31it/s]
    GMM time=0.292 seconds (corrected)
    Shape of embedding after round 0 is (160, 32)
    Shape of embedding after round 1 is (165, 32)
    Shape of embedding after round 2 is (170, 32)
    Shape of embedding after round 3 is (174, 32)
    Shape of embedding after round 4 is (177, 32)
    Shape of embedding after round 5 is (180, 32)
    Stored anchored map embedding in `adata.obsm['X_cplearn_coremap']` (180/180 points assigned).
    2025-11-23 00:13:51,072 - INFO - ✓ CoreAnalysis complete
    2025-11-23 00:13:51,072 - INFO -   - Core map embedding stored in: `adata.obsm['X_cplearn_coremap']` (shape: (180, 32))
    2025-11-23 00:13:51,072 - INFO -   - Assigned points: 180/180 (100.0%)
    2025-11-23 00:13:51,072 - INFO - 
    ============================================================
    2025-11-23 00:13:51,072 - INFO - DEG: Marker Genes
    2025-11-23 00:13:51,072 - INFO - ============================================================
    2025-11-23 00:13:51,072 - INFO - Identifying differentially expressed genes (marker genes)...
    Comparing groups {0} vs {1}
    2025-11-23 00:13:51,088 - INFO - ✓ DEG analysis complete
    2025-11-23 00:13:51,088 - INFO -   - Total differentially expressed genes: 50
    2025-11-23 00:13:51,088 - INFO -   - Significant genes (p_adj < 0.05): 41
    2025-11-23 00:13:51,088 - INFO -   - Significant genes (p_adj < 0.01): 39
    2025-11-23 00:13:51,088 - INFO - 
    Top 10 differentially expressed genes:
    2025-11-23 00:13:51,091 - INFO - 
        gene   log2fc   z_score       pvalue        p_adj   mean_a   mean_b  pct_expr_a  pct_expr_b
    gene_002 2.187279 18.400315 2.021035e-20 3.368391e-19 6.666667 0.683333    1.000000    0.500000
    gene_027 2.102686 19.505699 4.567841e-21 1.170854e-19 7.733333 1.033333    1.000000    0.600000
    gene_021 2.078502 16.554285 4.683415e-21 1.170854e-19 9.700000 1.533333    1.000000    0.733333
    gene_029 1.939962 12.955106 7.978928e-18 4.432738e-17 4.883333 0.533333    0.966667    0.366667
    gene_032 1.911840 13.144602 6.643614e-18 4.152259e-17 5.083333 0.616667    0.966667    0.433333
    gene_031 1.902821 15.913483 3.142358e-20 3.927948e-19 7.850000 1.366667    1.000000    0.750000
    gene_045 1.866375 14.216678 6.965233e-19 4.975166e-18 6.900000 1.166667    1.000000    0.716667
    gene_004 1.847075 14.939477 2.073976e-19 1.728313e-18 4.216667 0.450000    0.983333    0.350000
    gene_023 1.603167 10.278753 1.255096e-15 4.827291e-15 4.316667 0.750000    0.983333    0.500000
    gene_044 1.584963 12.506953 1.768928e-17 8.040582e-17 5.950000 1.316667    1.000000    0.666667
    2025-11-23 00:13:51,093 - INFO -   - DEG results saved to: /tmp/lotus_real_output/deg_results.csv
    2025-11-23 00:13:51,093 - INFO - 
    ============================================================
    2025-11-23 00:13:51,093 - INFO - Visualization: Rendering all visualizations
    2025-11-23 00:13:51,093 - INFO - ============================================================
    WARNING: saving figure to file figures/umap_clusters.png
    WARNING: saving figure to file figures/violin_markers.png
    2025-11-23 00:13:51,487 - INFO -   - Moved violin_markers.png to output directory
    2025-11-23 00:13:51,487 - INFO -   - Moved umap_clusters.png to output directory
    2025-11-23 00:13:51,487 - INFO - ✓ All visualizations complete
    2025-11-23 00:13:51,487 - INFO -   - UMAP clusters plot: /tmp/lotus_real_output/umap_clusters.png
    2025-11-23 00:13:51,487 - INFO -   - Marker genes violin plot: /tmp/lotus_real_output/violin_markers.png
    2025-11-23 00:13:51,487 - INFO -   - Top 5 marker genes visualized: gene_002, gene_027, gene_021, gene_029, gene_032
    2025-11-23 00:13:51,487 - INFO - 
    ============================================================
    2025-11-23 00:13:51,487 - INFO - Workflow Summary
    2025-11-23 00:13:51,487 - INFO - ============================================================
    2025-11-23 00:13:51,487 - INFO - ✓ All results saved to: /private/tmp/lotus_real_output
    2025-11-23 00:13:51,487 - INFO - ✓ Output directory: lotus_real_output
    2025-11-23 00:13:51,487 - INFO -   - Log file: /tmp/lotus_real_output/workflow_20251123_001335.log
    2025-11-23 00:13:51,487 - INFO -   - DEG results: /tmp/lotus_real_output/deg_results.csv
    2025-11-23 00:13:51,487 - INFO -   - UMAP plot: /tmp/lotus_real_output/umap_clusters.png
    2025-11-23 00:13:51,487 - INFO -   - Violin plot: /tmp/lotus_real_output/violin_markers.png
    2025-11-23 00:13:51,487 - INFO - ============================================================
    2025-11-23 00:13:51,487 - INFO - Workflow completed successfully!
    2025-11-23 00:13:51,487 - INFO - ============================================================

Output Files
------------

After running the complete workflow, you'll find the following output files in the ``./results/`` directory:

- ``umap_clusters.png`` - UMAP visualization colored by cluster labels
- ``deg_results.csv`` - Complete differential expression results (if saved)

Next Steps
----------

* See :doc:`api/index` for the complete API reference
* Run the example script: ``examples/lotus_workflow.py`` for a more detailed example with logging
* View the complete example on GitHub: `lotus_workflow.py <https://github.com/CrossOmics/Lotus/blob/main/examples/lotus_workflow.py>`_
