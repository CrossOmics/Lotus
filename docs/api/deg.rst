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

   **Compatibility:**
   
   This function works with cluster labels from both Lotus cplearn (``"cplearn_labels"``) and scanpy (``"leiden"``, ``"louvain"``), and auto-detects cluster keys if not specified.

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

   The result DataFrame contains the following columns: ``gene`` (gene name), ``log2fc`` (log2 fold change, where > 0 means upregulated and < 0 means downregulated), ``z_score`` (Z score), ``pvalue`` (p value), ``p_adj`` (adjusted p value using FDR, more reliable), ``mean_a`` and ``mean_b`` (mean expression in the two groups), and ``pct_expr_a`` and ``pct_expr_b`` (expression percentage).

.. autofunction:: lotus.workflows.deg.deg.run_differential_expression

   This is an alias for ``marker_genes()`` for backward compatibility.

Helper Functions
----------------

.. autofunction:: lotus.workflows.deg.deg.pick_groups

   Automatically select two comparison groups (first two non-negative clusters).
