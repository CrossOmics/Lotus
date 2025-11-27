Workflows Module
================

The ``lotus.workflows`` module provides complete single-cell analysis workflow functions, organized according to standard analysis steps.

Overview
--------

The Workflows module includes the following submodules:

- **Preprocessing**: Data quality control, filtering, normalization, and dimensionality reduction
- **CoreAnalysis**: Core map embedding (typically after clustering)
- **Clustering**: Cell clustering analysis
  - Supports multiple methods: Lotus cplearn (default) or scanpy (Leiden/Louvain)
  - Easy to switch between methods
- **Visualization**: Result visualization
- **DEG (Differential Expression Analysis)**: Marker gene identification

All Functions
-------------

.. automodule:: lotus.workflows
   :members:
   :undoc-members:
   :show-inheritance:

Preprocessing Functions
-----------------------

.. automodule:: lotus.workflows.preprocess.preprocess
   :members:
   :undoc-members:
   :show-inheritance:

Core Analysis Functions
------------------------

.. automodule:: lotus.workflows.core_analysis.core_analysis
   :members:
   :undoc-members:
   :show-inheritance:

Clustering Functions
--------------------

.. automodule:: lotus.workflows.clustering.clustering
   :members:
   :undoc-members:
   :show-inheritance:

Visualization Functions
-----------------------

.. automodule:: lotus.workflows.visualization.visualization
   :members:
   :undoc-members:
   :show-inheritance:

Differential Expression Analysis Functions
------------------------------------------

.. automodule:: lotus.workflows.deg.deg
   :members:
   :undoc-members:
   :show-inheritance:
