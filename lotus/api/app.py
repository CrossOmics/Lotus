"""
Main Flask application for Lotus API
"""

from flask import Flask, send_from_directory
from flask_cors import CORS
from pathlib import Path

from .config import UPLOAD_FOLDER, MAX_CONTENT_LENGTH, LOTUS_AVAILABLE, SCANPY_AVAILABLE
from .info import bp as info_bp
from .upload import bp as upload_bp
from .preprocess import bp as preprocess_bp
from .cluster import bp as cluster_bp
from .visualize import bp as visualize_bp
from .analysis import bp as analysis_bp
from .core_selection import bp as core_selection_bp


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
    app.register_blueprint(info_bp)
    app.register_blueprint(upload_bp)
    app.register_blueprint(preprocess_bp)
    app.register_blueprint(cluster_bp)
    app.register_blueprint(visualize_bp)
    app.register_blueprint(analysis_bp)
    app.register_blueprint(core_selection_bp)
    
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
    print("Lotus Embedding Projector Web Application")
    print("=" * 60)
    print(f"Lotus available: {LOTUS_AVAILABLE}")
    print(f"Scanpy available: {SCANPY_AVAILABLE}")
    print(f"Upload folder: {UPLOAD_FOLDER}")
    print("\nStarting server on http://localhost:5000")
    print("Open http://localhost:5000 in your browser")
    print("=" * 60)
    print("\n[INFO] Server logs will appear below:\n")
    
    # Try to find static folder (Lotus-Web-Demo directory)
    # When running from Lotus-Web-Demo/app.py, we need to go up one level
    # When running as module, we need to go up three levels
    current_file = Path(__file__)
    # Try different paths
    possible_paths = [
        current_file.parent.parent.parent / 'Lotus-Web-Demo',  # From lotus/api/app.py
        current_file.parent.parent / 'Lotus-Web-Demo',  # Alternative
        Path.cwd() / 'Lotus-Web-Demo',  # Current working directory
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
    app.run(debug=True, host='0.0.0.0', port=5000, use_reloader=False)


if __name__ == '__main__':
    main()

