Tools Module (lotus.tl)
========================

The ``lotus.tl`` module provides scanpy-compatible tools and analysis functions. This module is a complete wrapper around ``scanpy.tl``, providing all scanpy tools functions with identical interfaces.

Overview
--------

``lotus.tl`` serves as a building block API that provides direct access to all scanpy tools functions. All functions work exactly the same way as in scanpy, making Lotus fully compatible with scanpy workflows. This module is used internally by ``lotus.workflows`` to implement high-level functionality, and can also be used directly when you need fine-grained control or want to access advanced features not covered by the workflows module.

**Key Features:**

The module provides complete scanpy.tl compatibility with the same function signatures and behavior. It can be used as a drop-in replacement for scanpy.tl, and all scanpy tools functions are available. For standard analysis workflows, we recommend using the high-level functions in ``lotus.workflows`` instead, which provide optimized defaults and a streamlined interface.

Compatibility
-------------

You can use ``lotus.tl`` exactly like ``scanpy.tl``:

.. code-block:: python

   import lotus as lt
   
   # These work exactly like scanpy.tl
   lt.tl.leiden(adata, resolution=0.5)
   lt.tl.umap(adata)
   lt.tl.paga(adata, groups='leiden')
   lt.tl.rank_genes_groups(adata, groupby='leiden')

Main Functions
--------------

Clustering
~~~~~~~~~~

.. autofunction:: lotus.methods.scanpy.tools.leiden

   Leiden clustering algorithm (default clustering method in scanpy).

.. autofunction:: lotus.methods.scanpy.tools.louvain

   Louvain clustering algorithm.

Dimensionality Reduction
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.methods.scanpy.tools.umap

   UMAP dimensionality reduction.

.. autofunction:: lotus.methods.scanpy.tools.tsne

   t-SNE dimensionality reduction.

.. autofunction:: lotus.methods.scanpy.tools.diffmap

   Diffusion map dimensionality reduction.

.. autofunction:: lotus.methods.scanpy.tools.draw_graph

   Force-directed graph layout (e.g., Fruchterman-Reingold).

Trajectory Inference
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.methods.scanpy.tools.paga

   Partition-based graph abstraction (PAGA) for trajectory inference.

.. autofunction:: lotus.methods.scanpy.tools.dpt

   Diffusion pseudotime (DPT) analysis.

Differential Expression
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.methods.scanpy.tools.rank_genes_groups

   Rank genes for characterizing groups (differential expression analysis).

Gene Scoring
~~~~~~~~~~~~

.. autofunction:: lotus.methods.scanpy.tools.score_genes

   Score cells based on gene expression.

.. autofunction:: lotus.methods.scanpy.tools.score_genes_cell_cycle

   Score cells based on cell cycle phase genes.

Other Functions
~~~~~~~~~~~~~~~~

.. autofunction:: lotus.methods.scanpy.tools.dendrogram

   Compute hierarchical clustering dendrogram.

.. autofunction:: lotus.methods.scanpy.tools.ingest

   Map labels from reference dataset to query dataset.

Usage Examples
--------------

Clustering
~~~~~~~~~~

.. code-block:: python

   import lotus as lt
   
   # Leiden clustering
   lt.tl.leiden(adata, resolution=0.5, key_added='leiden')
   
   # Louvain clustering
   lt.tl.louvain(adata, resolution=0.5, key_added='louvain')

Dimensionality Reduction
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import lotus as lt
   
   # UMAP
   lt.tl.umap(adata, n_components=2, min_dist=0.5)
   
   # t-SNE
   lt.tl.tsne(adata, n_pcs=40, perplexity=30)

Trajectory Inference
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import lotus as lt
   
   # PAGA
   lt.tl.paga(adata, groups='leiden')
   
   # Diffusion pseudotime
   lt.tl.dpt(adata, n_dcs=10)

Differential Expression
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import lotus as lt
   
   # Rank genes by groups
   lt.tl.rank_genes_groups(
       adata,
       groupby='leiden',
       method='wilcoxon',
       n_genes=100
   )

Complete Function List
----------------------

All functions from ``scanpy.tl`` are available in ``lotus.tl``. Functions marked with ⚙️ are also available as high-level wrappers in ``lotus.workflows``, which we recommend for most users. For advanced features not covered by workflows, use the building blocks API directly.

**Clustering Functions:**

The ``leiden()`` function provides Leiden clustering algorithm (also available as ``lotus.workflows.clustering(method='leiden')``). The ``louvain()`` function provides Louvain clustering algorithm (also available as ``lotus.workflows.clustering(method='louvain')``).

**Dimensionality Reduction Functions:**

The ``umap()`` function computes UMAP embedding (also available as ``lotus.workflows.umap()``). The ``tsne()`` function provides t-SNE embedding for alternative dimensionality reduction. The ``diffmap()`` function computes diffusion map embedding. The ``draw_graph()`` function generates force-directed graph layouts such as Fruchterman-Reingold.

**Trajectory Inference Functions:**

The ``paga()`` function performs partition-based graph abstraction for trajectory inference. The ``dpt()`` function computes diffusion pseudotime analysis.

**Differential Expression Functions:**

The ``rank_genes_groups()`` function performs differential expression analysis (similar functionality available as ``lotus.workflows.marker_genes()`` but with different interface and more options).

**Gene Scoring Functions:**

The ``score_genes()`` function scores cells based on gene expression. The ``score_genes_cell_cycle()`` function scores cells based on cell cycle phase genes.

**Other Functions:**

The ``dendrogram()`` function computes hierarchical clustering dendrograms. The ``ingest()`` function maps labels from reference dataset to query dataset. The ``marker_gene_overlap()`` function analyzes marker gene overlap. The ``embedding_density()`` function computes embedding density. The ``filter_rank_genes_groups()`` function filters ranked genes. The ``sim()`` function simulates single-cell data.

For detailed documentation of each function, please refer to the `scanpy official documentation <https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.html>`__.

Note
----

All functions in ``lotus.tl`` are direct wrappers around ``scanpy.tl`` functions. They have identical signatures, parameters, and behavior. You can use them as drop-in replacements for scanpy.tl functions.

