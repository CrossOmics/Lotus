Quick Start
============

This guide will walk you through a complete single-cell RNA sequencing data analysis workflow, from data loading to final visualization results.

Complete Workflow
-----------------

Here is a typical Lotus analysis workflow:

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
    import pandas as pd
    import numpy as np

    # 1. Load data (assuming you already have an AnnData object)
    # adata = lt.read_10x_h5("path/to/data.h5")
    # Or load from other formats
    # adata = lt.read("path/to/data.h5ad")

    # 2. Preprocessing: QC → Filtering → Normalization → HVG → Scaling → PCA → Neighbors
    preprocess(
        adata,
        n_pcs=20,
        target_sum=1e4,
        n_top_genes=2000,
        n_neighbors=15,
        save_raw=True,
    )

    # 3. Core selection (optional, before clustering)
    # Note: Core selection can be performed before clustering to identify core cells
    # For cplearn workflow, this step is typically done after clustering (see below)
    
    # 4. Clustering analysis
    # Option A: Using Lotus cplearn (default)
    model = clustering(
        adata,
        use_rep="X_latent",
        key_added="cplearn_labels",
        cluster_resolution=1.2,
    )
    
    # Option B: Using scanpy (alternative)
    # import scanpy as sc
    # sc.tl.leiden(adata, resolution=0.5, key_added="leiden")
    # # Or use louvain: sc.tl.louvain(adata, resolution=0.5, key_added="louvain")

    # 5. Core selection (for cplearn workflow, after clustering)
    # This step requires the model from cplearn clustering
    core_selection(
        adata,
        model=model,
        use_rep="X_latent",
        key_added="X_cplearn_coremap",
    )

    # 6. UMAP visualization
    umap(
        adata,
        cluster_key="cplearn_labels",  # Use "leiden" or "louvain" if using scanpy
        output_dir="./results",
        save="_clusters.png",
    )

    # 7. Differential expression analysis (finding marker genes)
    de_result = marker_genes(
        adata,
        cluster_key="cplearn_labels",  # Use "leiden" or "louvain" if using scanpy
        layer="raw_counts",
        auto_pick_groups=True,
    )

    print("Analysis complete!")

Workflow Overview
----------------

1. **Preprocessing**
   - Calculate quality control metrics (QC metrics)
   - Filter low-quality cells and genes
   - Data normalization and log transformation
   - Select highly variable genes (HVG)
   - Principal component analysis (PCA)
   - Build neighbor graph

2. **Core Selection** (optional, typically before clustering for cplearn workflow)
   - Identify core cells for stable clustering
   - Compute core map embedding
   - Used to improve clustering stability

3. **Clustering**
   - **Lotus cplearn** (default): Perform cell clustering using cplearn algorithm
     - Auto-detect best representation (X_latent, X_pca, or X)
     - Generate stable clustering results
   - **scanpy** (alternative): Use scanpy's leiden or louvain algorithms
     - Directly compatible with Lotus workflow
     - Simply use ``sc.tl.leiden()`` or ``sc.tl.louvain()`` instead
   - Output cluster labels to adata.obs

4. **Core Selection** (for cplearn workflow, after clustering)
   - Compute core map embedding using clustering model
   - For further analysis and visualization

5. **Visualization**
   - Compute UMAP dimensionality reduction
   - Generate visualization plots of clustering results

6. **Differential Expression Analysis (DEG)**
   - Identify differentially expressed genes between clusters
   - Auto-select comparison groups (if not specified)
   - Output marker gene list

Next Steps
----------

* Check out :doc:`user_guide/index` for detailed explanations of each step
* See :doc:`api/index` for the complete API reference
* Run the example script: ``examples/lotus_workflow.py``
