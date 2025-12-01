"""
Upload API endpoints
"""

from flask import Blueprint, request, jsonify
import pandas as pd
from anndata import AnnData
from pathlib import Path
from ..config import UPLOAD_FOLDER, LOTUS_AVAILABLE, read, read_10x_h5, read_10x_mtx, sc
from ..utils import get_session_dir, save_adata

bp = Blueprint('upload', __name__, url_prefix='/api')

# Default datasets directory - check multiple possible locations
def get_default_datasets_dir():
    """Get the directory containing default datasets"""
    # Try multiple possible locations
    possible_paths = [
        Path(__file__).parent.parent.parent.parent / 'datasets',  # From lotus/api/upload/__init__.py
        Path(__file__).parent.parent.parent / 'datasets',  # Alternative
        Path.cwd() / 'datasets',  # Current working directory
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
        
        # Validate file type (h5ad, csv, tsv, mtx allowed)
        allowed_types = ['h5ad', 'csv', 'tsv', 'mtx']
        if file_type not in allowed_types:
            print(f"[UPLOAD] Error: Unsupported file type: {file_type}")
            return jsonify({'error': f'Unsupported file format. Only .h5ad, .csv, .tsv, and .mtx files are allowed.'}), 400
        
        # Additional validation: check file extension
        file_ext = Path(file.filename).suffix.lower()
        allowed_extensions = ['.h5ad', '.csv', '.tsv', '.mtx']
        if file_ext not in allowed_extensions:
            print(f"[UPLOAD] Error: Invalid file extension: {file_ext}")
            return jsonify({'error': f'Invalid file extension. Only .h5ad, .csv, .tsv, and .mtx files are allowed.'}), 400
        
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
            elif file_type == 'mtx':
                # MTX format - 10x Genomics format
                # read_10x_mtx expects a folder containing matrix.mtx, genes.tsv/features.tsv, and barcodes.tsv
                print(f"[UPLOAD] Loading MTX file: {filepath.name}")
                if not LOTUS_AVAILABLE or read_10x_mtx is None:
                    if sc is None:
                        return jsonify({'error': 'Lotus or scanpy not available. Please install lotus or scanpy to read mtx files.'}), 500
                    # Use scanpy's read_10x_mtx
                    read_10x_mtx_func = sc.read_10x_mtx
                else:
                    read_10x_mtx_func = read_10x_mtx
                
                mtx_folder = filepath.parent
                
                # Check if uploaded file is not named matrix.mtx - create a link or copy
                expected_mtx_name = mtx_folder / 'matrix.mtx'
                expected_mtx_gz_name = mtx_folder / 'matrix.mtx.gz'
                
                # If the uploaded file is not matrix.mtx, create a link or copy
                if filepath.name not in ['matrix.mtx', 'matrix.mtx.gz']:
                    if not expected_mtx_name.exists() and not expected_mtx_gz_name.exists():
                        # Create a symlink or copy the file as matrix.mtx
                        import os
                        try:
                            if filepath.suffix == '.gz' or filepath.name.endswith('.mtx.gz'):
                                if not expected_mtx_gz_name.exists():
                                    os.link(str(filepath), str(expected_mtx_gz_name))
                                    print(f"[UPLOAD] Created link: {filepath.name} -> matrix.mtx.gz")
                            else:
                                if not expected_mtx_name.exists():
                                    os.link(str(filepath), str(expected_mtx_name))
                                    print(f"[UPLOAD] Created link: {filepath.name} -> matrix.mtx")
                        except OSError:
                            # If link fails (e.g., on Windows or cross-filesystem), copy the file
                            import shutil
                            if filepath.suffix == '.gz' or filepath.name.endswith('.mtx.gz'):
                                if not expected_mtx_gz_name.exists():
                                    shutil.copy2(str(filepath), str(expected_mtx_gz_name))
                                    print(f"[UPLOAD] Copied file: {filepath.name} -> matrix.mtx.gz")
                            else:
                                if not expected_mtx_name.exists():
                                    shutil.copy2(str(filepath), str(expected_mtx_name))
                                    print(f"[UPLOAD] Copied file: {filepath.name} -> matrix.mtx")
                
                # Check for required files
                genes_file = mtx_folder / 'genes.tsv'
                features_file = mtx_folder / 'features.tsv'
                barcodes_file = mtx_folder / 'barcodes.tsv'
                
                if not (expected_mtx_name.exists() or expected_mtx_gz_name.exists()):
                    return jsonify({'error': 'matrix.mtx file not found. Please ensure the uploaded .mtx file is accessible, or upload it in a folder with the correct structure.'}), 400
                
                if not (genes_file.exists() or features_file.exists()):
                    return jsonify({'error': 'genes.tsv or features.tsv file not found. 10x mtx format requires: matrix.mtx, genes.tsv/features.tsv, and barcodes.tsv in the same folder. Please upload all required files or use a zip file containing the complete folder structure.'}), 400
                
                if not barcodes_file.exists():
                    return jsonify({'error': 'barcodes.tsv file not found. 10x mtx format requires: matrix.mtx, genes.tsv/features.tsv, and barcodes.tsv in the same folder. Please upload all required files or use a zip file containing the complete folder structure.'}), 400
                
                print(f"[UPLOAD] Reading mtx from folder: {mtx_folder}")
                adata = read_10x_mtx_func(str(mtx_folder))
                print(f"[UPLOAD] Loaded: {adata.shape[0]} cells, {adata.shape[1]} genes")
            else:
                return jsonify({'error': f'Unsupported file type: {file_type}. Only .h5ad, .csv, .tsv, and .mtx are supported.'}), 400
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
        
        # Find the dataset file
        datasets_dir = get_default_datasets_dir()
        filepath = datasets_dir / filename
        
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
            error_msg = f'Default dataset file not found: {filename}. Please ensure the file exists in the datasets directory or upload folder.'
            print(f"[LOAD] Error: {error_msg}")
            print(f"[LOAD] Searched in: {datasets_dir}, {UPLOAD_FOLDER}")
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

