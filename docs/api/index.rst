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
   </div>

Overview
--------

Lotus API provides high-level, workflow-oriented functions with complete analysis pipelines covering standard single-cell analysis steps including preprocessing, clustering, visualization, differential expression analysis, and core analysis. These functions offer a clean interface with smart defaults and are recommended for most users.

.. toctree::
   :maxdepth: 1

   workflows
   preprocess
   core_selection
   clustering
   visualization
   deg 

