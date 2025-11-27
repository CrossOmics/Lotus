Tools Module (lotus.tl)
========================

The ``lotus.tl`` module provides scanpy-compatible tools and analysis functions. This module is a complete wrapper around ``scanpy.tl``, providing all scanpy tools functions with identical interfaces.

Overview
--------

``lotus.tl`` is a direct wrapper around ``scanpy.tl``. All functions work exactly the same way as in scanpy, making Lotus fully compatible with scanpy workflows.

**Key Features:**

The module provides complete scanpy.tl compatibility with the same function signatures and behavior. It can be used as a drop-in replacement for scanpy.tl, and all scanpy tools functions are available.

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

All functions from ``scanpy.tl`` are available in ``lotus.tl``:

- ``leiden()`` - Leiden clustering
- ``louvain()`` - Louvain clustering
- ``umap()`` - UMAP embedding
- ``tsne()`` - t-SNE embedding
- ``diffmap()`` - Diffusion map
- ``draw_graph()`` - Force-directed graph layout
- ``paga()`` - PAGA trajectory inference
- ``dpt()`` - Diffusion pseudotime
- ``rank_genes_groups()`` - Differential expression
- ``score_genes()`` - Gene scoring
- ``score_genes_cell_cycle()`` - Cell cycle scoring
- ``dendrogram()`` - Hierarchical clustering
- ``ingest()`` - Label transfer
- ``marker_gene_overlap()`` - Marker gene overlap analysis
- ``embedding_density()`` - Embedding density
- ``filter_rank_genes_groups()`` - Filter ranked genes
- ``sim()`` - Simulate data

For detailed documentation of each function, please refer to the `scanpy official documentation <https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.html>`__.

Note
----

All functions in ``lotus.tl`` are direct wrappers around ``scanpy.tl`` functions. They have identical signatures, parameters, and behavior. You can use them as drop-in replacements for scanpy.tl functions.

