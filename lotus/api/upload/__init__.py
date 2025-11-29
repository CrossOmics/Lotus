"""
Upload API endpoints
"""

from flask import Blueprint, request, jsonify
import pandas as pd
from anndata import AnnData
from pathlib import Path
from ..config import UPLOAD_FOLDER, LOTUS_AVAILABLE, read, read_10x_h5, sc
from ..utils import get_session_dir, save_adata

bp = Blueprint('upload', __name__, url_prefix='/api')

# Default datasets directory - check multiple possible locations
def get_default_datasets_dir():
    """Get the directory containing default datasets"""
    # Try multiple possible locations, prioritize 'data' directory
    possible_paths = [
        Path(__file__).parent.parent.parent.parent / 'data',  # From lotus/api/upload/__init__.py -> data
        Path(__file__).parent.parent.parent / 'data',  # Alternative -> data
        Path.cwd() / 'data',  # Current working directory -> data
        Path(__file__).parent.parent.parent.parent / 'datasets',  # From lotus/api/upload/__init__.py -> datasets
        Path(__file__).parent.parent.parent / 'datasets',  # Alternative -> datasets
        Path.cwd() / 'datasets',  # Current working directory -> datasets
        UPLOAD_FOLDER / 'datasets',  # In upload folder
        UPLOAD_FOLDER,  # Directly in upload folder
    ]
    
    for path in possible_paths:
        if path.exists() and path.is_dir():
            return path
    
    # If no datasets directory found, use upload folder as fallback
    return UPLOAD_FOLDER


@bp.route('/upload', methods=['POST'])
def upload_data():
    """Upload and process data file"""
    try:
        print(f"[UPLOAD] Request received. Files: {list(request.files.keys())}, Form: {list(request.form.keys())}")
        
        if 'file' not in request.files:
            print("[UPLOAD] Error: No file in request")
            return jsonify({'error': 'No file provided'}), 400
        
        file = request.files['file']
        file_type = request.form.get('type', 'h5ad')
        
        print(f"[UPLOAD] File: {file.filename}, Type: {file_type}")
        
        if file.filename == '':
            print("[UPLOAD] Error: Empty filename")
            return jsonify({'error': 'No file selected'}), 400
        
        # Validate file type (only h5ad, csv, tsv allowed)
        allowed_types = ['h5ad', 'csv', 'tsv']
        if file_type not in allowed_types:
            print(f"[UPLOAD] Error: Unsupported file type: {file_type}")
            return jsonify({'error': f'Unsupported file format. Only .h5ad, .csv, and .tsv files are allowed.'}), 400
        
        # Additional validation: check file extension
        file_ext = Path(file.filename).suffix.lower()
        allowed_extensions = ['.h5ad', '.csv', '.tsv']
        if file_ext not in allowed_extensions:
            print(f"[UPLOAD] Error: Invalid file extension: {file_ext}")
            return jsonify({'error': f'Invalid file extension. Only .h5ad, .csv, and .tsv files are allowed.'}), 400
        
        # Save uploaded file
        filepath = UPLOAD_FOLDER / file.filename
        print(f"[UPLOAD] Saving to: {filepath}")
        file.save(filepath)
        print(f"[UPLOAD] File saved. Size: {filepath.stat().st_size} bytes")
        
        # Load data based on file type
        print(f"[UPLOAD] Loading file type: {file_type}")
        try:
            if file_type == 'h5ad':
                print(f"[UPLOAD] Loading h5ad file. Lotus available: {LOTUS_AVAILABLE}")
                if LOTUS_AVAILABLE:
                    adata = read(filepath)
                else:
                    if sc is None:
                        return jsonify({'error': 'scanpy not available. Please install scanpy.'}), 500
                    adata = sc.read(filepath)
                print(f"[UPLOAD] Loaded: {adata.shape[0]} cells, {adata.shape[1]} genes")
            elif file_type in ['csv', 'tsv']:
                # CSV/TSV - create AnnData from matrix
                print(f"[UPLOAD] Loading {file_type.upper()} file")
                df = pd.read_csv(filepath, sep='\t' if file_type == 'tsv' else ',', index_col=0)
                adata = AnnData(df.T)  # Transpose: genes as vars, cells as obs
                print(f"[UPLOAD] Loaded: {adata.shape[0]} cells, {adata.shape[1]} genes")
            else:
                return jsonify({'error': f'Unsupported file type: {file_type}. Only .h5ad, .csv, and .tsv are supported.'}), 400
        except Exception as e:
            import traceback
            error_msg = f'Failed to load file: {str(e)}'
            print(f"[UPLOAD] Load error: {error_msg}")
            print(traceback.format_exc())
            return jsonify({'error': error_msg}), 400
        
        # Store session data
        session_id = request.form.get('session_id', 'default')
        session_dir = get_session_dir(session_id)
        print(f"[UPLOAD] Session dir: {session_dir}")
        
        # Clear all previous data and results when loading new data
        print(f"[UPLOAD] Clearing previous data and results...")
        if session_dir.exists():
            import shutil
            try:
                # Remove entire session directory to clear everything
                shutil.rmtree(session_dir)
                print(f"[UPLOAD] Cleared session directory: {session_dir}")
            except Exception as e:
                print(f"[UPLOAD] Warning: Could not fully clear session directory: {e}")
                # Fallback: just remove the adata file
                adata_path = session_dir / 'data.h5ad'
                if adata_path.exists():
                    adata_path.unlink()
                    print(f"[UPLOAD] Removed old adata file: {adata_path}")
        
        # Save new AnnData (this will create the session directory)
        save_adata(adata, session_id)
        print(f"[UPLOAD] Saved successfully")
        
        result = {
            'success': True,
            'session_id': session_id,
            'shape': list(adata.shape),
            'obs_columns': list(adata.obs.columns),
            'obsm_keys': list(adata.obsm.keys()),
            'message': f'Loaded {adata.shape[0]} cells and {adata.shape[1]} genes'
        }
        print(f"[UPLOAD] Success: {result['message']}")
        return jsonify(result)
    
    except Exception as e:
        import traceback
        error_msg = f'Upload failed: {str(e)}'
        print(f"[UPLOAD] Exception: {error_msg}")
        print(traceback.format_exc())
        return jsonify({'error': error_msg}), 500


@bp.route('/load-default-dataset', methods=['POST'])
def load_default_dataset():
    """Load a default dataset by filename"""
    try:
        data = request.get_json()
        if not data:
            print("[LOAD] Error: No JSON data in request")
            return jsonify({'error': 'No data provided'}), 400
        
        filename = data.get('filename')
        if not filename:
            print("[LOAD] Error: No filename provided")
            return jsonify({'error': 'No filename provided'}), 400
        
        session_id = data.get('session_id', 'default')
        
        print(f"[LOAD] Loading default dataset: {filename}, session: {session_id}")
        
        # Validate filename (security: only allow specific filenames)
        allowed_datasets = ['demo_data.h5ad', 'pbmc3k_raw.h5ad']
        if filename not in allowed_datasets:
            print(f"[LOAD] Error: Dataset not allowed: {filename}")
            return jsonify({'error': f'Dataset not allowed. Allowed datasets: {", ".join(allowed_datasets)}'}), 400
        
        # Find the dataset file - prioritize 'data' directory
        datasets_dir = get_default_datasets_dir()
        filepath = datasets_dir / filename
        
        # Also check in data directory directly
        if not filepath.exists():
            possible_data_paths = [
                Path(__file__).parent.parent.parent.parent / 'data' / filename,
                Path(__file__).parent.parent.parent / 'data' / filename,
                Path.cwd() / 'data' / filename,
            ]
            for loc in possible_data_paths:
                if loc.exists() and loc.is_file():
                    filepath = loc
                    break
        
        # Also check in upload folder directly
        if not filepath.exists():
            filepath = UPLOAD_FOLDER / filename
        
        # Also check in parent directories
        if not filepath.exists():
            possible_locations = [
                Path(__file__).parent.parent.parent.parent / filename,
                Path(__file__).parent.parent.parent / filename,
                Path.cwd() / filename,
            ]
            for loc in possible_locations:
                if loc.exists() and loc.is_file():
                    filepath = loc
                    break
        
        if not filepath.exists():
            # Collect all searched paths for error message
            searched_paths = [
                str(datasets_dir),
                str(UPLOAD_FOLDER),
            ]
            # Add data directory paths
            data_paths = [
                Path(__file__).parent.parent.parent.parent / 'data',
                Path(__file__).parent.parent.parent / 'data',
                Path.cwd() / 'data',
            ]
            for p in data_paths:
                if str(p) not in searched_paths:
                    searched_paths.append(str(p))
            
            error_msg = f'Default dataset file not found: {filename}. Please ensure the file exists in the data directory or upload folder.'
            print(f"[LOAD] Error: {error_msg}")
            print(f"[LOAD] Searched in: {', '.join(searched_paths)}")
            return jsonify({'error': error_msg}), 404
        
        print(f"[LOAD] Found dataset file: {filepath}")
        
        # Load the dataset
        try:
            if LOTUS_AVAILABLE:
                adata = read(filepath)
            else:
                if sc is None:
                    return jsonify({'error': 'scanpy not available. Please install scanpy.'}), 500
                adata = sc.read(filepath)
            print(f"[LOAD] Loaded: {adata.shape[0]} cells, {adata.shape[1]} genes")
        except Exception as e:
            import traceback
            error_msg = f'Failed to load dataset: {str(e)}'
            print(f"[LOAD] Load error: {error_msg}")
            print(traceback.format_exc())
            return jsonify({'error': error_msg}), 400
        
        # Store session data
        session_dir = get_session_dir(session_id)
        print(f"[LOAD] Session dir: {session_dir}")
        
        # Clear all previous data and results when loading new data
        print(f"[LOAD] Clearing previous data and results...")
        if session_dir.exists():
            import shutil
            try:
                # Remove entire session directory to clear everything
                shutil.rmtree(session_dir)
                print(f"[LOAD] Cleared session directory: {session_dir}")
            except Exception as e:
                print(f"[LOAD] Warning: Could not fully clear session directory: {e}")
                # Fallback: just remove the adata file
                adata_path = session_dir / 'data.h5ad'
                if adata_path.exists():
                    adata_path.unlink()
                    print(f"[LOAD] Removed old adata file: {adata_path}")
        
        # Save new AnnData (this will create the session directory)
        save_adata(adata, session_id)
        print(f"[LOAD] Saved successfully")
        
        result = {
            'success': True,
            'session_id': session_id,
            'shape': list(adata.shape),
            'obs_columns': list(adata.obs.columns),
            'obsm_keys': list(adata.obsm.keys()),
            'message': f'Loaded {adata.shape[0]} cells and {adata.shape[1]} genes from {filename}'
        }
        print(f"[LOAD] Success: {result['message']}")
        return jsonify(result)
    
    except Exception as e:
        import traceback
        error_msg = f'Load default dataset failed: {str(e)}'
        print(f"[LOAD] Exception: {error_msg}")
        print(traceback.format_exc())
        return jsonify({'error': error_msg}), 500

