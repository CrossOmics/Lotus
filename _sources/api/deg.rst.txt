Differential Expression Analysis Module
========================================

The ``lotus.workflows.deg_analysis`` module provides differential expression analysis (DEG) and marker gene identification functionality.

Overview
--------

Differential expression analysis is used to identify differentially expressed genes between different cell populations, which are often called "marker genes". The module provides both cplearn-based and scanpy-based DEG analysis functions.

**Citation:** 

* **Wilcoxon Rank-Sum Test**: F. Wilcoxon. *Individual comparisons by ranking methods.* Biometrics Bulletin 1, 80-83 (1945).

Main Functions
--------------

Cplearn-Based Marker Gene Identification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.workflows.deg_analysis.marker_genes
   
   **Compatibility:**
   
   This function works with cluster labels from both Lotus cplearn (``"cplearn"``) and scanpy (``"leiden"``, ``"louvain"``), and auto-detects cluster keys if not specified.

   **Usage Example:**

   .. code-block:: python

      from lotus.workflows import marker_genes
      
      # Auto-select comparison groups
      de_result = marker_genes(
          adata,
          cluster_key="cplearn",
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

Scanpy-Based Differential Expression Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Rank Genes Groups
^^^^^^^^^^^^^^^^^

.. autofunction:: lotus.workflows.deg_analysis.rank_genes_groups

   ⚙️ Rank genes for characterizing groups using scanpy's differential expression analysis. Supports multiple statistical methods (wilcoxon, t-test, logreg, etc.) and can compare multiple groups simultaneously.

   **Usage Example:**

   .. code-block:: python

      from lotus.workflows import rank_genes_groups
      
      # Rank genes for all groups
      rank_genes_groups(
          adata,
          groupby="leiden",
          method="wilcoxon",
          n_genes=100,
      )
      
      # Access results
      print(adata.uns["rank_genes_groups"]["names"].head())
      
      # Compare specific groups
      rank_genes_groups(
          adata,
          groupby="leiden",
          groups=["0", "1"],
          reference="2",
          method="t-test",
      )

Filter Rank Genes Groups
^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: lotus.workflows.deg_analysis.filter_rank_genes_groups

   ⚙️ Filter ranked genes groups based on expression criteria to identify high-quality marker genes.

   **Usage Example:**

   .. code-block:: python

      from lotus.workflows import rank_genes_groups, filter_rank_genes_groups
      
      # First rank genes
      rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
      
      # Then filter based on expression criteria
      filter_rank_genes_groups(
          adata,
          min_in_group_fraction=0.25,
          min_fold_change=2,
          max_out_group_fraction=0.5,
      )
      
      # Access filtered results
      print(adata.uns["rank_genes_groups_filtered"]["names"].head())

Marker Gene Overlap
^^^^^^^^^^^^^^^^^^^^

.. autofunction:: lotus.workflows.deg_analysis.marker_gene_overlap

   ⚙️ Analyze marker gene overlap between data and reference markers to assess similarity and validate cell type annotations.

   **Usage Example:**

   .. code-block:: python

      from lotus.workflows import rank_genes_groups, marker_gene_overlap
      
      # First rank genes
      rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
      
      # Define reference markers
      reference_markers = {
          "T cells": {"CD3D", "CD3E", "CD3G"},
          "B cells": {"CD79A", "CD79B", "MS4A1"},
          "Monocytes": {"CD14", "FCGR3A", "LYZ"},
      }
      
      # Analyze overlap
      overlap_result = marker_gene_overlap(
          adata,
          reference_markers,
          top_n_markers=50,
      )
      
      print(overlap_result)

Helper Functions
----------------

.. autofunction:: lotus.workflows.deg_analysis.pick_groups

   Automatically select two comparison groups (first two non-negative clusters).
