User Guide
===========

This user guide provides detailed explanations of Lotus's various functional modules, helping you understand the purpose of each step and how to use them.

.. toctree::
   :maxdepth: 2

   preprocessing
   core_selection
   clustering
   deg
   visualization

Overview
--------

Lotus provides a complete single-cell RNA sequencing data analysis workflow, organized according to standard analysis steps:

1. **Preprocessing**: Data quality control, filtering, normalization, and dimensionality reduction
2. **Core Selection** (optional): Identify core cells for stable clustering (typically before clustering)
3. **Clustering**: Identify cell types and subpopulations
   - Supports multiple methods: Lotus cplearn (default) or scanpy (Leiden/Louvain)
   - Easy to switch between methods
4. **Differential Expression Analysis (DEG)**: Find marker genes
5. **Visualization**: Display analysis results

Each module can be used independently or combined into a complete analysis workflow.
