Differential Expression Analysis Module
========================================

The ``lotus.workflows.deg`` module provides differential expression analysis (DEG) and marker gene identification functionality.

Overview
--------

Differential expression analysis is used to identify differentially expressed genes between different cell populations, which are often called "marker genes".

Main Functions
--------------

Marker Gene Identification
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.deg.deg.marker_genes

   **Biological Background:**
   
   Marker genes are:
   - Genes highly expressed in specific cell types
   - Used to identify and distinguish different cell types
   - Help understand cell function and characteristics
   
   Differential expression analysis identifies by comparing gene expression between different cell populations:
   - Upregulated genes: Highly expressed in the target population
   - Downregulated genes: Lowly expressed in the target population
   - Significance: Determined through statistical tests

   **Usage Example:**

   .. code-block:: python

      from lotus.workflows import marker_genes
      
      # Auto-select comparison groups
      de_result = marker_genes(
          adata,
          cluster_key="cplearn_labels",
          layer="raw_counts",
          auto_pick_groups=True,
      )
      
      # View results
      print(de_result.head())
      
      # Filter significant genes
      significant = de_result[
          (de_result["p_adj"] < 0.05) & 
          (abs(de_result["log2fc"]) > 1)
      ]

   **Result Interpretation:**

   The result DataFrame contains the following columns:
   - ``gene``: Gene name
   - ``log2fc``: log2 fold change (> 0 means upregulated, < 0 means downregulated)
   - ``z_score``: Z score
   - ``pvalue``: p value
   - ``p_adj``: Adjusted p value (FDR, more reliable)
   - ``mean_a``, ``mean_b``: Mean expression in the two groups
   - ``pct_expr_a``, ``pct_expr_b``: Expression percentage

.. autofunction:: lotus.workflows.deg.deg.run_differential_expression

   This is an alias for ``marker_genes()`` for backward compatibility.

Helper Functions
----------------

.. autofunction:: lotus.workflows.deg.deg.pick_groups

   Automatically select two comparison groups (first two non-negative clusters).
