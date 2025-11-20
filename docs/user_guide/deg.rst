Differential Expression Analysis
=================================

Differential expression analysis is used to identify differentially expressed genes between different cell populations, which are often called "marker genes".

Overview
--------

The goals of differential expression analysis are:
- Identify differentially expressed genes between different clusters
- Find marker genes for each cell type
- Understand biological characteristics of different cell populations

Lotus differential expression analysis features:
- Auto-detect cluster keys and data layers to use
- Support automatic selection of comparison groups
- Output results in scanpy-compatible format
- Works with cluster labels from both Lotus cplearn and scanpy (Leiden/Louvain)

Basic Usage
-----------

.. code-block:: python

    from lotus.workflows import marker_genes
    
    # Automatically select first two clusters for comparison
    de_result = marker_genes(
        adata,
        cluster_key="cplearn_labels",  # Cluster key (optional, auto-detected)
        layer="raw_counts",            # Data layer to use (optional, auto-detected)
        auto_pick_groups=True,         # Automatically select comparison groups
    )
    
    # View results
    print(de_result.head())

Manually Specifying Comparison Groups
-------------------------------------

If you want to compare specific clusters:

.. code-block:: python

    de_result = marker_genes(
        adata,
        cluster_key="cplearn_labels",
        groups_a={0},  # First group: cluster 0
        groups_b={1}, # Second group: cluster 1
        layer="raw_counts",
    )

Parameter Description
---------------------

- ``cluster_key``: Key name for cluster labels in ``adata.obs``
  - If ``None``, auto-detects: ``"cplearn_labels"`` > ``"leiden"`` > ``"louvain"``

- ``groups_a``: Set of cluster labels for first comparison group, e.g., ``{0, 1}``

- ``groups_b``: Set of cluster labels for second comparison group, e.g., ``{2, 3}``

- ``layer``: Data layer to use
  - If ``None``, auto-detects: ``"raw_counts"`` > ``"raw"`` > ``"counts"`` > ``None`` (use X)

- ``auto_pick_groups``: Whether to automatically select first two non-negative clusters as comparison groups, default True

- ``min_detect_pct``: Minimum detection percentage, default 0.0

- ``min_cells_per_group``: Minimum number of cells per group, default 5

Result Interpretation
---------------------

The differential expression analysis result is a DataFrame containing the following columns:

- ``gene``: Gene name
- ``log2fc``: log2 fold change
- ``z_score``: Z score
- ``pvalue``: p value
- ``p_adj``: Adjusted p value (FDR)
- ``mean_a``: Mean expression in group A
- ``mean_b``: Mean expression in group B
- ``pct_expr_a``: Expression percentage in group A
- ``pct_expr_b``: Expression percentage in group B

Filtering Significant Genes
----------------------------

.. code-block:: python

    # Filter significantly differentially expressed genes (p_adj < 0.05, |log2fc| > 1)
    significant = de_result[
        (de_result["p_adj"] < 0.05) & 
        (abs(de_result["log2fc"]) > 1)
    ]
    
    print(f"Number of significantly differentially expressed genes: {len(significant)}")
    
    # View top 10 upregulated genes
    top_up = significant.nlargest(10, "log2fc")
    print(top_up[["gene", "log2fc", "p_adj"]])

Saving Results
---------------

.. code-block:: python

    de_result.to_csv("deg_results.csv", index=False)

Compatibility with scanpy
-------------------------

Result format is compatible with scanpy's ``rank_genes_groups`` and can be used for subsequent analysis.

Frequently Asked Questions
---------------------------

**Q: Which data layer should be used for analysis?**
A: Usually use raw count data (``"raw_counts"`` or ``"raw"``), because differential expression analysis requires unnormalized data.

**Q: How to interpret log2fc?**
A: log2fc > 0 means upregulated in group A, log2fc < 0 means downregulated in group A. |log2fc| > 1 is usually considered a significant fold change.

**Q: What's the difference between p_adj and pvalue?**
A: p_adj is the p value after multiple testing correction (FDR), which is more reliable. Usually use p_adj < 0.05 as the significance threshold.

**Q: How to find marker genes for each cluster?**
A: You can compare each cluster with all other clusters separately, or use scanpy's ``sc.tl.rank_genes_groups()`` function.

Next Steps
----------

After completing differential expression analysis, you can:
- Perform visualization: :doc:`visualization`
- Use marker genes for further analysis
