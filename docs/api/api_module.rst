Lotus API Module (Web API)
============================

The ``lotus.api`` module provides Flask REST API endpoints for the Lotus Web Demo application. This module is designed for web-based interaction with Lotus workflows through HTTP requests.

Overview
--------

The Lotus API module is a Flask-based REST API that wraps Lotus workflows for web applications. It provides HTTP endpoints for:

- Data upload and management
- Preprocessing workflows
- Clustering analysis
- Visualization generation
- Differential expression analysis
- Core analysis
- Session management

**Important:** This module is designed for web applications (Flask REST API), not for direct Python function calls. For Python programming, use ``lotus.workflows`` instead.

Main Components
---------------

Flask Application
~~~~~~~~~~~~~~~~~~

.. autofunction:: lotus.api.app.create_app

   Creates and configures the Flask application instance for the Lotus Web Demo.

   **Parameters:**
   
   - ``static_folder``: Path to static files folder (for serving HTML/JS files)
   
   **Returns:**
   
   Flask application instance with all API endpoints registered.

   **Usage Example:**
   
   .. code-block:: python
   
      from lotus.api import create_app
      
      # Create Flask app
      app = create_app(static_folder='./Lotus-Web-Demo')
      
      # Run the server
      if __name__ == '__main__':
          app.run(debug=True, host='0.0.0.0', port=5000)

API Endpoints
~~~~~~~~~~~~~

The Lotus API provides the following REST endpoints (all prefixed with ``/api``):

**Data Management:**
- ``POST /api/upload`` - Upload single-cell data (h5ad format)
- ``POST /api/session/create`` - Create a new analysis session
- ``POST /api/session/heartbeat`` - Keep session alive
- ``POST /api/session/cleanup`` - Clean up expired sessions

**Analysis Workflows:**
- ``POST /api/preprocess`` - Run preprocessing pipeline
- ``POST /api/cluster`` - Run clustering analysis
- ``POST /api/visualize`` - Generate UMAP visualization
- ``POST /api/marker-genes`` - Run differential expression analysis
- ``POST /api/core-selection`` - Run core analysis

**Information:**
- ``GET /api/info/health`` - Check API health status
- ``GET /api/info/metadata`` - Get dataset metadata
- ``POST /api/available-cluster-keys`` - Get available cluster keys

**Ground Truth:**
- ``POST /api/ground-truth/upload`` - Upload ground truth labels
- ``GET /api/ground-truth/list`` - List available ground truth labels
- ``GET /api/ground-truth/get`` - Get ground truth labels

Configuration
--------------

The API module uses configuration from ``lotus.api.config``:

- **Upload Folder**: Temporary directory for storing uploaded data
- **Max File Size**: 500MB default limit
- **Memory Optimization**: Default parameters optimized for 4GB server memory

Session Management
------------------

The API uses session-based data management:

- Each analysis session has a unique ``session_id``
- Data is stored in session-specific directories
- Sessions expire after a configurable timeout
- Automatic cleanup of expired sessions

Compatibility
-------------

The API module:
- Works with both Lotus workflows and scanpy (fallback)
- Supports standard AnnData format (h5ad files)
- Compatible with scanpy-compatible data structures
- Handles memory optimization for web server environments

Note
----

For Python programming and direct function calls, use ``lotus.workflows`` instead of ``lotus.api``. The API module is specifically designed for web applications and HTTP-based interactions.

