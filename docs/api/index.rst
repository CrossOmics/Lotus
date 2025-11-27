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
               <span>Interactive API Module</span>
           </div>
           <p class="lotus-card-description">Flask REST API documentation used for building interactive analysis and visualization tool, Lotus Embedding Projector.</p>
       </a>
       <a href="tools.html" class="lotus-card">
           <div class="lotus-card-title">
               <span class="lotus-card-icon">🔧</span>
               <span>Tools Module</span>
           </div>
           <p class="lotus-card-description">Tools functions for clustering, dimensionality reduction, trajectory inference, and more.</p>
       </a>
       <a href="preprocessing_compat.html" class="lotus-card">
           <div class="lotus-card-title">
               <span class="lotus-card-icon">⚙️</span>
               <span>Preprocessing Module</span>
           </div>
           <p class="lotus-card-description">Scanpy-compatible preprocessing: QC, filtering, normalization, batch correction, and more.</p>
       </a>
   </div>

Overview
--------

Lotus API is organized into two main layers. The **Workflows Module** provides high-level, workflow-oriented functions with complete analysis pipelines covering standard single-cell analysis steps including preprocessing, clustering, visualization, differential expression analysis, and core analysis. These functions offer a clean interface with smart defaults and are recommended for most users. The workflows module internally uses Tools Module(``lotus.tl``) and Preprocessing Module(``lotus.pp``) as building blocks to implement its functionality.

The **Building Blocks API** is consists of Tools Module(``lotus.tl``) and Preprocessing Module(``lotus.pp``). These modules serve as low-level building blocks that can be used directly when you need fine-grained control or want to access advanced features not covered by the workflows module, such as batch correction and trajectory inference (PAGA, DPT).  For standard analysis workflows, we recommend using the Workflows Module functions such as ``lotus.workflows.preprocess()``, ``lotus.workflows.clustering()``, ``lotus.workflows.umap()``, and ``lotus.workflows.marker_genes()``. For advanced features, use the Building Blocks API directly. 

The **Web API Module** (``lotus.api``) is a Flask REST API that powers the Lotus Embedding Projector, an interactive tool for single-cell data analysis and visualization.  We recommend local deployment for better performance and to support larger datasets. For detailed setup instructions and deployment guide, please refer to the `Lotus Web Demo documentation <https://github.com/CrossOmics/Lotus/tree/main/Lotus-Web-Demo>`__.
