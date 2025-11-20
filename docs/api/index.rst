API Reference
=============

Complete API reference documentation for Lotus. All functions and classes are organized by module.

.. toctree::
   :maxdepth: 2

   workflows
   preprocess
   core_selection
   clustering
   visualization
   deg

Overview
--------

Lotus API is mainly divided into two parts:

1. **Workflows Module** (Main API): Provides high-level, workflow-oriented functions
2. **Compatibility API**: Low-level functions compatible with scanpy

For most users, we recommend using functions from the Workflows module, which provide a cleaner interface and smart defaults.

Workflows Module
----------------

The Workflows module is the core of Lotus, providing a complete single-cell analysis workflow:

- :doc:`workflows` - Overview of all workflows functions
- :doc:`preprocess` - Preprocessing related functions
- :doc:`core_selection` - Core selection functions (typically before clustering)
- :doc:`clustering` - Clustering analysis functions
- :doc:`visualization` - Visualization functions
- :doc:`deg` - Differential expression analysis functions

Compatibility API
-----------------

Lotus also provides scanpy-compatible APIs in the main ``lotus`` module. These functions are wrappers around scanpy and are used exactly the same way as scanpy. For detailed API documentation, please refer to the `scanpy official documentation <https://scanpy.readthedocs.io/>`__.
