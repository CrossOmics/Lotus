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

    # 3. Clustering analysis
    model = clustering(
        adata,
        use_rep="X_latent",
        key_added="cplearn_labels",
        cluster_resolution=1.2,
    )

    # 4. UMAP visualization
    umap(
        adata,
        cluster_key="cplearn_labels",
        output_dir="./results",
        save="_clusters.png",
    )

    # 5. Differential expression analysis (finding marker genes)
    de_result = marker_genes(
        adata,
        cluster_key="cplearn_labels",
        layer="raw_counts",
        auto_pick_groups=True,
    )

    # 6. Core selection (optional)
    core_selection(
        adata,
        model=model,
        use_rep="X_latent",
        key_added="X_cplearn_coremap",
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

2. **Clustering**
   - Perform cell clustering using cplearn algorithm
   - Auto-detect best representation (X_latent, X_pca, or X)
   - Output cluster labels to adata.obs

3. **Visualization**
   - Compute UMAP dimensionality reduction
   - Generate visualization plots of clustering results

4. **Differential Expression Analysis (DEG)**
   - Identify differentially expressed genes between clusters
   - Auto-select comparison groups (if not specified)
   - Output marker gene list

5. **Core Selection**
   - Compute core map embedding
   - For further analysis and visualization

Next Steps
----------

* Check out :doc:`user_guide/index` for detailed explanations of each step
* See :doc:`api/index` for the complete API reference
* Run the example script: ``examples/lotus_workflow.py``
