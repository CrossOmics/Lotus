"""
Upload API endpoints
"""

from flask import Blueprint, request, jsonify
import pandas as pd
from anndata import AnnData
from pathlib import Path
import zipfile
import shutil
import tempfile
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
        
        # Validate file type (h5ad, csv, tsv, mtx, zip allowed)
        allowed_types = ['h5ad', 'csv', 'tsv', 'mtx', 'zip']
        if file_type not in allowed_types:
            print(f"[UPLOAD] Error: Unsupported file type: {file_type}")
            return jsonify({'error': f'Unsupported file format. Only .h5ad, .csv, .tsv, .mtx, and .zip (for mtx) files are allowed.'}), 400
        
        # Additional validation: check file extension
        file_ext = Path(file.filename).suffix.lower()
        allowed_extensions = ['.h5ad', '.csv', '.tsv', '.mtx', '.zip']
        if file_ext not in allowed_extensions:
            print(f"[UPLOAD] Error: Invalid file extension: {file_ext}")
            return jsonify({'error': f'Invalid file extension. Only .h5ad, .csv, .tsv, .mtx, and .zip (for mtx) files are allowed.'}), 400
        
        # For mtx format, file can be .mtx or .zip
        if file_type == 'mtx' and file_ext != '.zip' and file_ext != '.mtx':
            return jsonify({'error': 'MTX format requires a .zip file containing the mtx folder (with matrix.mtx, genes.tsv/features.tsv, and barcodes.tsv) or a single .mtx file with other files in the same folder.'}), 400
        
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
            elif file_type == 'mtx' or (file_type == 'zip' and file_ext == '.zip'):
                # MTX format - 10x Genomics format
                # read_10x_mtx expects a folder containing matrix.mtx, genes.tsv/features.tsv, and barcodes.tsv
                print(f"[UPLOAD] Loading MTX format (file: {filepath.name})")
                if not LOTUS_AVAILABLE or read_10x_mtx is None:
                    if sc is None:
                        return jsonify({'error': 'Lotus or scanpy not available. Please install lotus or scanpy to read mtx files.'}), 500
                    read_10x_mtx_func = sc.read_10x_mtx
                else:
                    read_10x_mtx_func = read_10x_mtx
                
                # Handle zip file or single mtx file
                extract_dir = None
                if file_ext == '.zip':
                    # Extract zip file to temporary directory
                    extract_dir = UPLOAD_FOLDER / f"mtx_extract_{filepath.stem}"
                    extract_dir.mkdir(exist_ok=True)
                    try:
                        with zipfile.ZipFile(filepath, 'r') as zip_ref:
                            zip_ref.extractall(extract_dir)
                        print(f"[UPLOAD] Extracted zip to: {extract_dir}")
                        mtx_folder = extract_dir
                    except Exception as e:
                        return jsonify({'error': f'Failed to extract zip file: {str(e)}'}), 400
                else:
                    # Single mtx file - use parent directory
                    mtx_folder = filepath.parent
                
                # Find mtx file and required files (support different naming patterns)
                # Use rglob for recursive search
                mtx_files = list(mtx_folder.rglob('*.mtx'))
                mtx_files.extend(list(mtx_folder.rglob('*.mtx.gz')))
                # Remove duplicates
                mtx_files = list(set(mtx_files))
                
                if not mtx_files:
                    return jsonify({'error': 'No .mtx file found in the uploaded folder. Please ensure the zip contains a matrix.mtx file (or data_matrix.mtx, counts_matrix.mtx, etc.)'}), 400
                
                # Find the mtx folder (could be root or subdirectory)
                mtx_file = mtx_files[0]
                actual_mtx_folder = mtx_file.parent
                print(f"[UPLOAD] Found mtx file: {mtx_file.name} in folder: {actual_mtx_folder}")
                
                # Check for required files with different naming patterns
                # Search in the same folder as mtx file
                genes_files = list(actual_mtx_folder.glob('genes.tsv')) + list(actual_mtx_folder.glob('features.tsv'))
                # Also check for data_* and counts_* patterns
                genes_files.extend(list(actual_mtx_folder.glob('*genes.tsv')) + list(actual_mtx_folder.glob('*features.tsv')))
                genes_files.extend(list(actual_mtx_folder.glob('data_*genes.tsv')) + list(actual_mtx_folder.glob('data_*features.tsv')))
                genes_files.extend(list(actual_mtx_folder.glob('counts_*genes.tsv')) + list(actual_mtx_folder.glob('counts_*features.tsv')))
                # Remove duplicates
                genes_files = list(set(genes_files))
                
                barcodes_files = list(actual_mtx_folder.glob('barcodes.tsv'))
                barcodes_files.extend(list(actual_mtx_folder.glob('*barcodes.tsv')))
                barcodes_files.extend(list(actual_mtx_folder.glob('data_*barcodes.tsv')))
                barcodes_files.extend(list(actual_mtx_folder.glob('counts_*barcodes.tsv')))
                # Remove duplicates
                barcodes_files = list(set(barcodes_files))
                
                print(f"[UPLOAD] Found {len(genes_files)} genes/features files, {len(barcodes_files)} barcodes files in {actual_mtx_folder}")
                if genes_files:
                    print(f"[UPLOAD] Genes/features files: {[f.name for f in genes_files]}")
                if barcodes_files:
                    print(f"[UPLOAD] Barcodes files: {[f.name for f in barcodes_files]}")
                
                # Define expected file names
                import os
                expected_mtx_name = actual_mtx_folder / 'matrix.mtx'
                expected_mtx_gz_name = actual_mtx_folder / 'matrix.mtx.gz'
                expected_genes_name = actual_mtx_folder / 'genes.tsv'
                expected_features_name = actual_mtx_folder / 'features.tsv'
                expected_barcodes_name = actual_mtx_folder / 'barcodes.tsv'
                
                # If mtx file is not named matrix.mtx, create a link/copy
                # Always create matrix.mtx (not .gz) for read_10x_mtx compatibility
                if mtx_file.name not in ['matrix.mtx', 'matrix.mtx.gz']:
                    # For read_10x_mtx, we need matrix.mtx (uncompressed) even if original is .gz
                    # But if original is .gz, we should keep it as .gz
                    if mtx_file.name.endswith('.gz'):
                        target_name = expected_mtx_gz_name
                    else:
                        target_name = expected_mtx_name
                    
                    if not target_name.exists():
                        try:
                            os.link(str(mtx_file), str(target_name))
                            print(f"[UPLOAD] Created link: {mtx_file.name} -> {target_name.name}")
                        except OSError as e:
                            print(f"[UPLOAD] Link failed ({e}), trying copy...")
                            shutil.copy2(str(mtx_file), str(target_name))
                            print(f"[UPLOAD] Copied file: {mtx_file.name} -> {target_name.name}")
                    else:
                        print(f"[UPLOAD] Target file already exists: {target_name.name}")
                elif mtx_file.name == 'matrix.mtx.gz':
                    # If it's already matrix.mtx.gz, that's fine
                    print(f"[UPLOAD] Mtx file already has correct name: {mtx_file.name}")
                else:
                    # If it's matrix.mtx, that's perfect
                    print(f"[UPLOAD] Mtx file already has correct name: {mtx_file.name}")
                
                # If genes/features files have different names, create links
                if genes_files:
                    genes_file = genes_files[0]
                    if 'features' in genes_file.name.lower():
                        target = expected_features_name
                    else:
                        target = expected_genes_name
                    if not target.exists():
                        try:
                            os.link(str(genes_file), str(target))
                            print(f"[UPLOAD] Created link: {genes_file.name} -> {target.name}")
                        except OSError:
                            shutil.copy2(str(genes_file), str(target))
                            print(f"[UPLOAD] Copied file: {genes_file.name} -> {target.name}")
                
                # If barcodes file has different name, create link
                if barcodes_files:
                    barcodes_file = barcodes_files[0]
                    if not expected_barcodes_name.exists():
                        try:
                            os.link(str(barcodes_file), str(expected_barcodes_name))
                            print(f"[UPLOAD] Created link: {barcodes_file.name} -> barcodes.tsv")
                        except OSError:
                            shutil.copy2(str(barcodes_file), str(expected_barcodes_name))
                            print(f"[UPLOAD] Copied file: {barcodes_file.name} -> barcodes.tsv")
                
                # Final check for required files
                mtx_exists = expected_mtx_name.exists() or expected_mtx_gz_name.exists()
                genes_exists = expected_genes_name.exists() or expected_features_name.exists()
                barcodes_exists = expected_barcodes_name.exists()
                
                print(f"[UPLOAD] File check - mtx: {mtx_exists}, genes/features: {genes_exists}, barcodes: {barcodes_exists}")
                print(f"[UPLOAD] Expected files in {actual_mtx_folder}:")
                print(f"[UPLOAD]   - matrix.mtx: {expected_mtx_name.exists()}")
                print(f"[UPLOAD]   - matrix.mtx.gz: {expected_mtx_gz_name.exists()}")
                print(f"[UPLOAD]   - genes.tsv: {expected_genes_name.exists()}")
                print(f"[UPLOAD]   - features.tsv: {expected_features_name.exists()}")
                print(f"[UPLOAD]   - barcodes.tsv: {expected_barcodes_name.exists()}")
                
                if not mtx_exists:
                    # List all files in the folder for debugging
                    all_files = list(actual_mtx_folder.glob('*'))
                    file_list = ', '.join([f.name for f in all_files])
                    return jsonify({'error': f'matrix.mtx file not found after processing. Files in folder: {file_list}. Please ensure the zip contains a .mtx file.'}), 400
                
                if not genes_exists:
                    return jsonify({'error': 'genes.tsv or features.tsv file not found. 10x mtx format requires: matrix.mtx, genes.tsv/features.tsv, and barcodes.tsv in the same folder.'}), 400
                
                if not barcodes_exists:
                    return jsonify({'error': 'barcodes.tsv file not found. 10x mtx format requires: matrix.mtx, genes.tsv/features.tsv, and barcodes.tsv in the same folder.'}), 400
                
                print(f"[UPLOAD] Reading mtx from folder: {actual_mtx_folder}")
                try:
                    adata = read_10x_mtx_func(str(actual_mtx_folder))
                    print(f"[UPLOAD] Loaded: {adata.shape[0]} cells, {adata.shape[1]} genes")
                finally:
                    # Clean up extracted files if from zip
                    if extract_dir is not None and extract_dir.exists():
                        try:
                            shutil.rmtree(extract_dir)
                            print(f"[UPLOAD] Cleaned up extracted files: {extract_dir}")
                        except Exception as e:
                            print(f"[UPLOAD] Warning: Could not clean up extracted files: {e}")
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

