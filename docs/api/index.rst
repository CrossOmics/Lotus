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
           <p class="lotus-card-description">Preprocessing functions.</p>
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
           <p class="lotus-card-description">Clustering analysis functions.</p>
       </a>
       <a href="visualization.html" class="lotus-card">
           <div class="lotus-card-title">
               <span class="lotus-card-icon">📊</span>
               <span>Visualization Module</span>
           </div>
           <p class="lotus-card-description">Visualization functions.</p>
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
           <p class="lotus-card-description">Flask REST API documentation for web applications and Lotus Web Demo.</p>
       </a>
       <a href="tools.html" class="lotus-card">
           <div class="lotus-card-title">
               <span class="lotus-card-icon">🔧</span>
               <span>Tools Module (lotus.tl)</span>
           </div>
           <p class="lotus-card-description">Tools functions for clustering, dimensionality reduction, trajectory inference, and more.</p>
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

Lotus API is organized into two main layers. The **Workflows Module** provides high-level, workflow-oriented functions with complete analysis pipelines covering standard single-cell analysis steps including preprocessing, clustering, visualization, differential expression analysis, and core analysis. These functions offer a clean interface with smart defaults and are recommended for most users. The workflows module internally uses ``lotus.tl`` and ``lotus.pp`` as building blocks to implement its functionality.

The **Building Blocks API** consists of ``lotus.tl`` and ``lotus.pp``, which provide all scanpy tools and preprocessing functions with identical interfaces. These modules serve as low-level building blocks that can be used directly when you need fine-grained control or want to access advanced features not covered by the workflows module, such as batch correction, trajectory inference (PAGA, DPT), t-SNE, and other specialized scanpy functions.

**Usage Recommendations:**

For standard analysis workflows including preprocessing, clustering, UMAP visualization, and differential expression analysis, we recommend using the **Workflows Module** functions such as ``lotus.workflows.preprocess()``, ``lotus.workflows.clustering()``, ``lotus.workflows.umap()``, and ``lotus.workflows.marker_genes()``. These provide optimized defaults and a streamlined interface.

For advanced features not covered by the workflows module, such as batch correction (ComBat), trajectory inference (PAGA, DPT), t-SNE dimensionality reduction, gene scoring, and other specialized scanpy functions, use the **Building Blocks API** directly: ``lotus.tl`` for analysis tools and ``lotus.pp`` for preprocessing functions. These functions work exactly the same way as scanpy, making Lotus fully compatible with scanpy workflows. For detailed API documentation of individual functions, please refer to the `scanpy official documentation <https://scanpy.readthedocs.io/>`__.

The **Web API Module** (``lotus.api``) is a Flask REST API that provides HTTP endpoints for web-based single-cell analysis. This module powers the Lotus Web Demo application and can be deployed locally for interactive data analysis through a web interface. We recommend local deployment for better performance and to support larger datasets. The Web API provides session-based data management and analysis workflows accessible via HTTP requests. For Python programming and direct function calls, use the Workflows Module or Building Blocks API instead. For detailed setup instructions and deployment guide, please refer to the `Lotus Web Demo documentation <https://github.com/CrossOmics/Lotus/tree/main/Lotus-Web-Demo>`__.

