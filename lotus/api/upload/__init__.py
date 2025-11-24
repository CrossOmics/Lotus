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


@bp.route('/upload', methods=['POST'])
def upload_data():
    """Upload and process data file"""
    try:
        print(f"[UPLOAD] Request received. Files: {list(request.files.keys())}, Form: {list(request.form.keys())}")
        
        if 'file' not in request.files:
            print("[UPLOAD] Error: No file in request")
            return jsonify({'error': 'No file provided'}), 400
        
        file = request.files['file']
        file_type = request.form.get('type', 'h5ad')  # h5ad, h5, mtx, csv, tsv
        
        print(f"[UPLOAD] File: {file.filename}, Type: {file_type}")
        
        if file.filename == '':
            print("[UPLOAD] Error: Empty filename")
            return jsonify({'error': 'No file selected'}), 400
        
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
            elif file_type == 'h5':
                print(f"[UPLOAD] Loading h5 file. Lotus available: {LOTUS_AVAILABLE}")
                if LOTUS_AVAILABLE:
                    adata = read_10x_h5(filepath)
                else:
                    if sc is None:
                        return jsonify({'error': 'scanpy not available. Please install scanpy.'}), 500
                    adata = sc.read_10x_h5(filepath)
            elif file_type == 'mtx':
                return jsonify({'error': 'MTX format requires additional files'}), 400
            else:
                # CSV/TSV - create AnnData from matrix
                print(f"[UPLOAD] Loading CSV/TSV file")
                df = pd.read_csv(filepath, sep='\t' if file_type == 'tsv' else ',', index_col=0)
                adata = AnnData(df.T)  # Transpose: genes as vars, cells as obs
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

