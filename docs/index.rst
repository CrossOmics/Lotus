Lotus Documentation
===================

Welcome to the Lotus documentation! Lotus is a Python package for single-cell RNA sequencing (scRNA-seq) data analysis, providing a complete analysis pipeline from data preprocessing to clustering, differential expression analysis, and visualization.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   getting_started
   user_guide/index
   api/index

About Lotus
-----------

Lotus is a Python package compatible with various single-cell analysis tools (such as scanpy), providing a unified and easy-to-use API. It organizes functions according to standard single-cell analysis workflows, making it easy for biologists to perform data analysis.

Key Features
------------

* **Complete Analysis Pipeline**: Full pipeline from data preprocessing to visualization
* **Method Flexibility**: Support multiple clustering methods (cplearn or scanpy), easy to switch between them
* **scanpy Compatible**: Fully compatible with the scanpy ecosystem
* **Easy to Use**: Clear API design suitable for biologists
* **Automated Processing**: Smart parameter defaults and auto-detection features

Acknowledgments
---------------

Lotus is built on top of excellent open-source packages in the single-cell analysis ecosystem:

* **scanpy** - Single-cell analysis in Python
* **anndata** - Annotated data objects for single-cell omics
* **cplearn** - Core-periphery learning for single-cell analysis
* **pandas** - Data manipulation and analysis
* **numpy** - Numerical computing
* **scikit-learn** - Machine learning tools
* **umap-learn** - Uniform Manifold Approximation and Projection

We thank all the developers and contributors of these packages for their excellent work.

Quick Start
-----------

If you're new to Lotus, we recommend starting with :doc:`getting_started`, which includes a complete workflow example.

.. note::

    Lotus is designed to make single-cell data analysis easy for biologists. If you encounter any issues, please refer to the relevant API documentation or user guide.
