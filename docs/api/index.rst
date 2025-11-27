API Reference
=============

Complete API reference documentation for Lotus. All functions and classes are organized by module.

.. raw:: html

   <div class="lotus-cards">
       <a href="workflows.html" class="lotus-card">
           <div class="lotus-card-title">
               <span class="lotus-card-icon">⚙️</span>
               <span>Workflows Module</span>
           </div>
           <p class="lotus-card-description">Overview of all workflows functions and high-level API.</p>
       </a>
       <a href="preprocess.html" class="lotus-card">
           <div class="lotus-card-title">
               <span class="lotus-card-icon">🔧</span>
               <span>Preprocessing Module</span>
           </div>
           <p class="lotus-card-description">Preprocessing functions: QC, filtering, normalization, HVG, scaling, PCA, and neighbors.</p>
       </a>
       <a href="core_selection.html" class="lotus-card">
           <div class="lotus-card-title">
               <span class="lotus-card-icon">🎯</span>
               <span>Core Analysis Module</span>
           </div>
           <p class="lotus-card-description">Core analysis functions for identifying core cells and computing core map embeddings.</p>
       </a>
       <a href="clustering.html" class="lotus-card">
           <div class="lotus-card-title">
               <span class="lotus-card-icon">🔀</span>
               <span>Clustering Module</span>
           </div>
           <p class="lotus-card-description">Clustering analysis functions using cplearn or scanpy-compatible methods.</p>
       </a>
       <a href="visualization.html" class="lotus-card">
           <div class="lotus-card-title">
               <span class="lotus-card-icon">📊</span>
               <span>Visualization Module</span>
           </div>
           <p class="lotus-card-description">Visualization functions for UMAP plots and marker gene visualizations.</p>
       </a>
       <a href="deg.html" class="lotus-card">
           <div class="lotus-card-title">
               <span class="lotus-card-icon">🧬</span>
               <span>Differential Expression Module</span>
           </div>
           <p class="lotus-card-description">Differential expression analysis functions for finding marker genes.</p>
       </a>
       <a href="api_module.html" class="lotus-card">
           <div class="lotus-card-title">
               <span class="lotus-card-icon">🌐</span>
               <span>Web API Module</span>
           </div>
           <p class="lotus-card-description">Flask REST API endpoints for web applications and Lotus Web Demo.</p>
       </a>
       <a href="tools.html" class="lotus-card">
           <div class="lotus-card-title">
               <span class="lotus-card-icon">🔧</span>
               <span>Tools Module (lotus.tl)</span>
           </div>
           <p class="lotus-card-description">Scanpy-compatible tools: clustering, dimensionality reduction, trajectory inference, and more.</p>
       </a>
       <a href="preprocessing_compat.html" class="lotus-card">
           <div class="lotus-card-title">
               <span class="lotus-card-icon">⚙️</span>
               <span>Preprocessing Module (lotus.pp)</span>
           </div>
           <p class="lotus-card-description">Scanpy-compatible preprocessing: QC, filtering, normalization, batch correction, and more.</p>
       </a>
   </div>

Overview
--------

Lotus API is organized into three main parts:

1. **Workflows Module** (Main API): High-level, workflow-oriented functions
   - Provides complete analysis pipelines
   - Clean interface with smart defaults
   - Recommended for most users
   - Includes: preprocessing, clustering, visualization, DEG analysis, core analysis

2. **Compatibility API**: Low-level functions compatible with scanpy
   - ``lotus.tl``: Tools module (clustering, dimensionality reduction, trajectory inference)
   - ``lotus.pp``: Preprocessing module (QC, filtering, normalization, batch correction)
   - Complete wrappers around scanpy with identical interfaces
   - Can be used as drop-in replacements for scanpy

3. **Web API Module**: Flask REST API for web applications
   - ``lotus.api``: HTTP endpoints for Lotus Web Demo
   - Session-based data management
   - Designed for web applications, not direct Python calls

For most users, we recommend using functions from the **Workflows Module**, which provide a cleaner interface and smart defaults.

Compatibility API Details
-------------------------

Lotus provides scanpy-compatible APIs through ``lotus.tl`` and ``lotus.pp``. These modules are complete wrappers around scanpy:

- **lotus.tl**: All functions from ``scanpy.tl`` (clustering, UMAP, PAGA, etc.)
- **lotus.pp**: All functions from ``scanpy.pp`` (preprocessing, batch correction, etc.)

These functions work exactly the same way as scanpy. For detailed API documentation of individual functions, please refer to the `scanpy official documentation <https://scanpy.readthedocs.io/>`__.

**Example:**

.. code-block:: python

   import lotus as lt
   
   # Use lotus.tl exactly like scanpy.tl
   lt.tl.leiden(adata, resolution=0.5)
   lt.tl.umap(adata)
   lt.tl.paga(adata, groups='leiden')
   
   # Use lotus.pp exactly like scanpy.pp
   lt.pp.normalize_total(adata, target_sum=1e4)
   lt.pp.log1p(adata)
   lt.pp.highly_variable_genes(adata, n_top_genes=2000)
   lt.pp.combat(adata, key='batch')  # Batch correction
