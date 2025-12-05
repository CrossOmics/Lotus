"""
Main Flask application for Lotus API
"""

# Set matplotlib backend BEFORE importing any matplotlib code
# This prevents GUI-related errors in Flask/threaded environments (especially on macOS)
# 'Agg' is a non-interactive backend that works on all platforms (Linux, Windows, macOS)
import matplotlib
matplotlib.use('Agg') 

from flask import Flask, send_from_directory
from flask_cors import CORS
from pathlib import Path

from .config import UPLOAD_FOLDER, MAX_CONTENT_LENGTH, LOTUS_AVAILABLE, SCANPY_AVAILABLE
from .get_info import bp as get_info_bp
from .upload import bp as upload_bp
from .preprocess import bp as preprocess_bp
from .cluster import bp as cluster_bp
from .visualize import bp as visualize_bp
from .core_analyze import bp as core_analyze_bp
from .deg_analyze import bp as deg_analyze_bp
from .session import bp as session_bp


def create_app(static_folder: str = None):
    """
    Create and configure Flask application
    
    Parameters:
        static_folder: Path to static files folder (for serving HTML/JS)
    
    Returns:
        Flask application instance
    """
    app = Flask(__name__, static_folder=static_folder, static_url_path='')
    CORS(app)
    
    # Configuration
    app.config['UPLOAD_FOLDER'] = str(UPLOAD_FOLDER)
    app.config['MAX_CONTENT_LENGTH'] = MAX_CONTENT_LENGTH
    
    # Register blueprints
    app.register_blueprint(get_info_bp)
    app.register_blueprint(upload_bp)
    app.register_blueprint(preprocess_bp)
    app.register_blueprint(cluster_bp)
    app.register_blueprint(visualize_bp)
    app.register_blueprint(core_analyze_bp)
    app.register_blueprint(deg_analyze_bp)
    app.register_blueprint(session_bp)
    
    # Serve static files if static_folder is provided
    if static_folder:
        @app.route('/')
        def index():
            """Serve the main HTML page"""
            response = send_from_directory(static_folder, 'index.html')
            response.headers['Cache-Control'] = 'no-cache, no-store, must-revalidate'
            response.headers['Pragma'] = 'no-cache'
            response.headers['Expires'] = '0'
            return response
        
        @app.route('/app.js')
        def app_js():
            """Serve the JavaScript file"""
            response = send_from_directory(static_folder, 'app.js')
            response.headers['Cache-Control'] = 'no-cache, no-store, must-revalidate'
            response.headers['Pragma'] = 'no-cache'
            response.headers['Expires'] = '0'
            return response
    
    return app


def main():
    """Main entry point for running the Flask server"""
    import logging
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger('werkzeug')
    logger.setLevel(logging.INFO)
    
    print("=" * 60)
    print("Lotus Web Application")
    print("=" * 60)
    print(f"Lotus available: {LOTUS_AVAILABLE}")
    print(f"Scanpy available: {SCANPY_AVAILABLE}")
    print(f"Upload folder: {UPLOAD_FOLDER}")
    
    # Try to find static folder (Interactive-Lotus directory)
    # When running from Interactive-Lotus/app.py, we need to go up one level
    # When running as module, we need to go up three levels
    current_file = Path(__file__)
    # Try different paths
    possible_paths = [
        current_file.parent.parent.parent / 'Interactive-Lotus',  # From lotus/api/app.py
        current_file.parent.parent / 'Interactive-Lotus',  # Alternative
        Path.cwd() / 'Interactive-Lotus',  # Current working directory
    ]
    
    static_folder = None
    for path in possible_paths:
        if path.exists() and path.is_dir():
            static_folder = path
            break
    
    # If still not found, try to use current directory if index.html exists
    if static_folder is None:
        cwd = Path.cwd()
        if (cwd / 'index.html').exists():
            static_folder = cwd
    
    app = create_app(static_folder=str(static_folder) if static_folder else None)
    
    # Get port from environment variable, default to 5259 (5000 is often used by AirPlay on macOS)
    port = int(os.environ.get('PORT', 5259))
    
    print(f"\nStarting server on http://localhost:{port}")
    print(f"Open http://localhost:{port} in your browser")
    print("=" * 60)
    print("\n[INFO] Server logs will appear below:\n")
    
    app.run(debug=True, host='0.0.0.0', port=port, use_reloader=False)


# Create app instance for Gunicorn production deployment
# Gunicorn will use this app object directly: gunicorn lotus.api.app:app
import os

# Check if running in production (Render sets PORT, Gunicorn sets GUNICORN)
# In Docker, serve static files; on Render/cloud, don't serve (frontend is on GitHub Pages)
if os.environ.get('PORT') or os.environ.get('GUNICORN'):
    # Check if STATIC_FOLDER is explicitly set
    static_folder_env = os.environ.get('STATIC_FOLDER')
    if static_folder_env:
        static_folder = static_folder_env
    elif os.environ.get('DOCKER') or os.path.exists('/app/Interactive-Lotus'):
        # In Docker, try to find static folder
        possible_paths = [
            Path('/app/Interactive-Lotus'),
            Path.cwd() / 'Interactive-Lotus',
            Path(__file__).parent.parent.parent / 'Interactive-Lotus',
        ]
        static_folder = None
        for path in possible_paths:
            if path.exists() and path.is_dir() and (path / 'index.html').exists():
                static_folder = str(path)
                break
    else:
        # On cloud platforms like Render, don't serve static files
        static_folder = None
    app = create_app(static_folder=static_folder)
else:
    # For development or when imported as module
    app = None

if __name__ == '__main__':
    main()

