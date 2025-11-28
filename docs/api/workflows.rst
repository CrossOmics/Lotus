Workflows Module
================

The ``lotus.workflows`` module provides complete single-cell analysis workflow functions, organized according to standard analysis steps.

Overview
--------

The Workflows module includes several submodules organized according to standard single-cell analysis steps. The **Preprocessing** submodule provides data quality control, filtering, normalization, and dimensionality reduction functionality. The **CoreAnalysis** submodule performs core map embedding before clustering to identify core cells. The **Clustering** submodule provides cell clustering analysis with support for multiple methods, including Lotus cplearn (default) and scanpy methods (Leiden/Louvain), with easy switching between methods. The **Visualization** submodule generates result visualizations, and the **DEG (Differential Expression Analysis)** submodule identifies marker genes.

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
