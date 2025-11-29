Lotus API Module (Web API)
============================

The ``lotus.api`` module provides Flask REST API endpoints for the Lotus Web Demo application. This module is designed for web-based interaction with Lotus workflows through HTTP requests.

Overview
--------

The Lotus API module is a Flask-based REST API that wraps Lotus workflows for web applications. It provides HTTP endpoints for data upload and management, preprocessing workflows, clustering analysis, visualization generation, differential expression analysis, core analysis, and session management.

**Important:** This module is designed for web applications (Flask REST API), not for direct Python function calls. For Python programming, use ``lotus.workflows`` instead.

Main Components
---------------

Flask Application
~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.api.app.create_app

   Creates and configures the Flask application instance for the Lotus Web Demo.

   **Parameters:**
   
   The ``static_folder`` parameter specifies the path to static files folder (for serving HTML/JS files).
   
   **Returns:**
   
   Flask application instance with all API endpoints registered.

   **Usage Example:**
   
   .. code-block:: python
   
      from lotus.api import create_app
      
      # Create Flask app
      app = create_app(static_folder='./Interactive-Lotus')
      
      # Run the server
      if __name__ == '__main__':
          app.run(debug=True, host='0.0.0.0', port=5000)

API Endpoints
~~~~~~~~~~~~~

The Lotus API provides the following REST endpoints (all prefixed with ``/api``):

**Data Management:** The ``POST /api/upload`` endpoint uploads single-cell data in h5ad format. The ``POST /api/session/create`` endpoint creates a new analysis session. The ``POST /api/session/heartbeat`` endpoint keeps a session alive. The ``POST /api/session/cleanup`` endpoint cleans up expired sessions.

**Analysis Workflows:** The ``POST /api/preprocess`` endpoint runs the preprocessing pipeline. The ``POST /api/cluster`` endpoint runs clustering analysis. The ``POST /api/visualize`` endpoint generates UMAP visualization. The ``POST /api/marker-genes`` endpoint runs differential expression analysis. The ``POST /api/core-selection`` endpoint runs core analysis.

**Information:** The ``GET /api/info/health`` endpoint checks API health status. The ``GET /api/info/metadata`` endpoint gets dataset metadata. The ``POST /api/available-cluster-keys`` endpoint gets available cluster keys.

**Ground Truth:** The ``POST /api/ground-truth/upload`` endpoint uploads ground truth labels. The ``GET /api/ground-truth/list`` endpoint lists available ground truth labels. The ``GET /api/ground-truth/get`` endpoint gets ground truth labels.

Configuration
--------------

The API module uses configuration from ``lotus.api.config``. The upload folder is a temporary directory for storing uploaded data, with a maximum file size limit of 500MB by default. Default parameters are optimized for 4GB server memory.

Session Management
------------------

The API uses session-based data management. Each analysis session has a unique ``session_id``, and data is stored in session-specific directories. Sessions expire after a configurable timeout, with automatic cleanup of expired sessions.

Compatibility
-------------

The API module works with both Lotus workflows and scanpy (fallback), supports standard AnnData format (h5ad files), is compatible with scanpy-compatible data structures, and handles memory optimization for web server environments.

Note
----

For Python programming and direct function calls, use ``lotus.workflows`` instead of ``lotus.api``. The API module is specifically designed for web applications and HTTP-based interactions.

