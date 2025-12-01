Standard Scanpy Workflow
=========================

This section demonstrates the standard single-cell analysis workflow using Lotus's high-level wrappers, which are compatible with scanpy's API.

The workflow follows these steps:

1. **Preprocess** - Quality control, filtering, normalization, PCA, and neighbor graph construction
2. **Clustering** - Identify cell clusters using Leiden or Louvain
3. **Visualization** - Generate UMAP visualization
4. **DEG Analysis** - Find marker genes between clusters

Import Required Modules
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import lotus as lt
    from lotus.workflows.preprocessing import preprocess
    from lotus.workflows.clustering import leiden  # or louvain
    from lotus.workflows.visualization import umap
    from lotus.workflows.deg_analysis import rank_genes_groups  # or marker_genes
    from anndata import AnnData

Step 1: Load Data
~~~~~~~~~~~~~~~~~

Load your single-cell data. Lotus supports various data formats through anndata:

.. code-block:: python

    # Load from 10x Genomics H5 format
    # adata = lt.read_10x_h5("path/to/data.h5")
    
    # Or load from AnnData H5AD format
    # adata = lt.read("path/to/data.h5ad")
    
    # Or use scanpy's read functions
    # import scanpy as sc
    # adata = sc.read_10x_mtx("path/to/folder")

Step 2: Preprocessing
~~~~~~~~~~~~~~~~~~~~~~

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
- Quality control metrics (QC metrics) calculation
- Filters low-quality cells and genes
- Normalizes data and applies log transformation
- Selects highly variable genes (HVG)
- Performs principal component analysis (PCA)
- Builds the neighbor graph

After preprocessing, check the results:

.. code-block:: python

    print(f" Preprocessing complete")
    print(f"  - Data shape: {adata.shape}")
    print(f"  - PCA stored in: `adata.obsm['X_pca']`")
    print(f"  - Neighbors graph constructed: {'neighbors' in adata.uns}")

Step 3: Clustering
~~~~~~~~~~~~~~~~~~

Lotus supports multiple clustering methods. For the standard scanpy workflow, use Leiden or Louvain:

**Option A: Leiden Clustering (Recommended)**

.. code-block:: python

    leiden(
        adata,
        resolution=0.5,
        key_added="leiden",
    )
    
    print(f"Found {adata.obs['leiden'].nunique()} clusters")

**Option B: Louvain Clustering**

.. code-block:: python

    from lotus.workflows.clustering import louvain
    
    louvain(
        adata,
        resolution=0.5,
        key_added="louvain",
    )
    
    print(f"Found {adata.obs['louvain'].nunique()} clusters")

Step 4: Visualization
~~~~~~~~~~~~~~~~~~~~~

Generate UMAP visualization of clustering results:

.. code-block:: python

    umap(
        adata,
        cluster_key="leiden",  # or "louvain"
        output_dir="./results",
        save="_clusters.png",
    )

This step computes UMAP dimensionality reduction and generates visualization plots. The UMAP embedding is stored in `adata.obsm['X_umap']`.

Step 5: Differential Expression Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Find marker genes between clusters using scanpy's method:

.. code-block:: python

    from lotus.workflows.deg_analysis import rank_genes_groups
    
    rank_genes_groups(
        adata,
        groupby="leiden",  # or "louvain"
        method="wilcoxon",
        key_added="rank_genes_groups",
    )
    
    # View top marker genes for each cluster
    import pandas as pd
    result = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)
    print(result)

Complete Standard Workflow Example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example uses the demo dataset included in the repository:

.. code-block:: python

    import lotus as lt
    from lotus.workflows.preprocessing import preprocess
    from lotus.workflows.clustering import leiden
    from lotus.workflows.visualization import umap
    from lotus.workflows.deg_analysis import rank_genes_groups
    
    # 1. Load data (using demo dataset from repository)
    adata = lt.read("data/demo_data.h5ad")
    # Or load your own data:
    # adata = lt.read("path/to/your/data.h5ad")
    
    # 2. Preprocessing
    preprocess(adata, n_pcs=20, n_top_genes=2000, n_neighbors=15, save_raw=True)
    
    # 3. Clustering (Leiden)
    leiden(adata, resolution=0.5, key_added="leiden")
    
    # 4. Visualization
    umap(adata, cluster_key="leiden", output_dir="./results", save="_clusters.png")
    
    # 5. Differential Expression Analysis
    rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
    
    print(" Standard workflow complete!")

Example output visualization:

The visualization below was generated using the demo dataset (`data/demo_data.h5ad`) with the standard workflow:

.. figure:: _static/examples/umap_standard_workflow.png
   :alt: UMAP visualization from standard workflow
   :width: 600px
   :align: center

   UMAP visualization colored by Leiden clusters from the standard workflow example, generated using `data/demo_data.h5ad`.

