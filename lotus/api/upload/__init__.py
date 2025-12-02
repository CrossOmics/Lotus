"""
Upload API endpoints
"""

from flask import Blueprint, request, jsonify
import pandas as pd
from anndata import AnnData
from pathlib import Path
import shutil
import tempfile
from ..config import UPLOAD_FOLDER, LOTUS_AVAILABLE, read, read_10x_h5, read_10x_mtx, sc
from ..utils import get_session_dir, save_adata

bp = Blueprint('upload', __name__, url_prefix='/api')

def handle_mtx_multi_upload():
    """Handle upload of multiple MTX files (matrix, features, barcodes)"""
    try:
        print(f"[UPLOAD] Handling MTX multi-file upload")
        
        # Check for required files
        if 'matrix_file' not in request.files:
            return jsonify({'error': 'No matrix file provided'}), 400
        if 'features_file' not in request.files:
            return jsonify({'error': 'No features file provided'}), 400
        if 'barcodes_file' not in request.files:
            return jsonify({'error': 'No barcodes file provided'}), 400
        
        matrix_file = request.files['matrix_file']
        features_file = request.files['features_file']
        barcodes_file = request.files['barcodes_file']
        session_id = request.form.get('session_id')
        
        if matrix_file.filename == '' or features_file.filename == '' or barcodes_file.filename == '':
            return jsonify({'error': 'One or more files are empty'}), 400
        
        print(f"[UPLOAD] MTX files: matrix={matrix_file.filename}, features={features_file.filename}, barcodes={barcodes_file.filename}")
        
        # Check if scanpy is available
        if sc is None:
            return jsonify({'error': 'scanpy not available. Please install scanpy to read mtx files.'}), 500
        
        # Create temporary directory for MTX files
        import tempfile
        import numpy as np
        temp_dir = Path(tempfile.mkdtemp(prefix='lotus_mtx_', dir=str(UPLOAD_FOLDER)))
        print(f"[UPLOAD] Created temp directory: {temp_dir}")
        
        try:
            # Save files
            matrix_path = temp_dir / matrix_file.filename
            features_path = temp_dir / 'features.tsv'
            barcodes_path = temp_dir / 'barcodes.tsv'
            
            matrix_file.save(matrix_path)
            features_file.save(features_path)
            barcodes_file.save(barcodes_path)
            
            # 1. 用 scanpy.read_mtx() 读取矩阵
            print(f"[UPLOAD] Reading matrix with scanpy.read_mtx: {matrix_path}")
            adata = sc.read_mtx(str(matrix_path))
            print(f"[UPLOAD] Matrix loaded: {adata.shape}")
            
            # 读取 MTX header 来确认格式
            with open(matrix_path, 'r') as f:
                header_line = f.readline().strip()
                dims_line = f.readline().strip()
                mtx_rows, mtx_cols = map(int, dims_line.split()[:2])
            print(f"[UPLOAD] MTX header: {mtx_rows} rows, {mtx_cols} columns")
            
            # 2. 判断 scanpy.read_mtx() 的格式
            # scanpy.read_mtx() 可能返回 (genes, cells) 或 (cells, genes)
            # 需要根据 MTX header 和实际 shape 来判断
            format_warning = None
            if adata.shape[0] == mtx_rows and adata.shape[1] == mtx_cols:
                # scanpy 没有转置：shape = (MTX行数, MTX列数)
                # 在 10x 格式中：MTX行数 = 基因，MTX列数 = 细胞
                # 所以：adata.shape = (genes, cells)，需要转置
                # 这可能是 Seurat 格式（行=基因，列=细胞）
                print(f"[UPLOAD] scanpy.read_mtx() returned (genes, cells) format, transposing...")
                print(f"[UPLOAD] WARNING: Detected Seurat format (genes x cells). Lotus is optimized for scanpy format (cells x genes).")
                format_warning = "Detected Seurat format (genes x cells). Data has been automatically transposed to scanpy format (cells x genes)."
                adata = adata.T
                print(f"[UPLOAD] After transpose: {adata.shape} (cells x genes)")
                n_genes_needed = adata.shape[1]
                n_cells_needed = adata.shape[0]
            elif adata.shape[0] == mtx_cols and adata.shape[1] == mtx_rows:
                # scanpy 已经转置：shape = (MTX列数, MTX行数) = (cells, genes)
                print(f"[UPLOAD] scanpy.read_mtx() returned (cells, genes) format, no transpose needed")
                n_genes_needed = adata.shape[1]
                n_cells_needed = adata.shape[0]
            else:
                # 无法确定，尝试根据 features/barcodes 行数判断
                print(f"[UPLOAD] Cannot determine format from MTX header, will infer from features/barcodes")
                n_genes_needed = None
                n_cells_needed = None
            
            # 3. 读取 features.tsv 第一列（只有一列）
            features_df = pd.read_csv(features_path, sep='\t', header=None)
            n_features = len(features_df)
            
            # 4. 读取 barcodes.tsv 第一列（只有一列）
            barcodes_df = pd.read_csv(barcodes_path, sep='\t', header=None)
            n_barcodes = len(barcodes_df)
            
            print(f"[UPLOAD] Features file: {n_features} rows")
            print(f"[UPLOAD] Barcodes file: {n_barcodes} rows")
            
            # 5. 根据 features 和 barcodes 的行数判断格式
            if n_genes_needed is None:
                # 需要推断格式
                if n_features == adata.shape[0] and n_barcodes == adata.shape[1]:
                    # shape = (genes, cells)，需要转置
                    # 这可能是 Seurat 格式
                    print(f"[UPLOAD] Inferred format: (genes, cells), transposing...")
                    print(f"[UPLOAD] WARNING: Detected Seurat format (genes x cells). Lotus is optimized for scanpy format (cells x genes).")
                    format_warning = "Detected Seurat format (genes x cells). Data has been automatically transposed to scanpy format (cells x genes)."
                    adata = adata.T
                    n_genes_needed = adata.shape[1]
                    n_cells_needed = adata.shape[0]
                elif n_features == adata.shape[1] and n_barcodes == adata.shape[0]:
                    # shape = (cells, genes)，不需要转置
                    print(f"[UPLOAD] Inferred format: (cells, genes), no transpose needed")
                    n_genes_needed = adata.shape[1]
                    n_cells_needed = adata.shape[0]
                else:
                    return jsonify({
                        'error': f'Cannot determine matrix format. Features has {n_features} rows, barcodes has {n_barcodes} rows, but matrix shape is {adata.shape}. Please check your data files.'
                    }), 400
            
            print(f"[UPLOAD] Final format: {adata.shape} (cells x genes)")
            print(f"[UPLOAD] Need {n_genes_needed} genes, {n_cells_needed} cells")
            
            # 6. 检查 features.tsv 行数
            if n_features < n_genes_needed:
                return jsonify({
                    'error': f'Data mismatch: features.tsv has {n_features} rows, but matrix has {n_genes_needed} genes. Please check your data files.'
                }), 400
            elif n_features > n_genes_needed:
                # 矩阵可能是过滤后的，只取前 n_genes_needed 行
                print(f"[UPLOAD] Features has {n_features} rows, matrix has {n_genes_needed} genes. Using first {n_genes_needed} rows.")
                features_df = features_df.iloc[:n_genes_needed].reset_index(drop=True)
                n_features = len(features_df)
            
            # 读取第一列
            gene_names = features_df.iloc[:, 0].values.astype(str)
            
            # 检查内容是否为空/无效，如果是则按顺序补齐（0, 1, 2, ...）
            # 检查：全部为空、全部是NaN、全部是空字符串、或全部相同（可能是占位符）
            is_empty = len(gene_names) == 0
            is_all_na = all(pd.isna(gene_names)) if len(gene_names) > 0 else True
            is_all_empty_str = all(gene_names == '') if len(gene_names) > 0 else True
            is_all_same = len(set(gene_names)) == 1 if len(gene_names) > 0 else True
            
            if is_empty or is_all_na or is_all_empty_str or (is_all_same and len(gene_names) > 0 and (pd.isna(gene_names[0]) or gene_names[0] == '')):
                print(f"[UPLOAD] Features file appears empty or invalid, using sequential labels (0-{n_genes_needed-1})")
                gene_names = np.array([str(i) for i in range(n_genes_needed)])
            else:
                # 检查是否有空值，如果有则用索引补齐
                mask = pd.isna(gene_names) | (gene_names == '')
                if mask.any():
                    print(f"[UPLOAD] Found {mask.sum()} empty values in features.tsv, filling with sequential indices")
                    gene_names = np.where(mask, [str(i) for i in range(n_genes_needed)], gene_names)
            
            adata.var_names = gene_names
            
            # 7. 检查 barcodes.tsv 行数
            if n_barcodes < n_cells_needed:
                return jsonify({
                    'error': f'Data mismatch: barcodes.tsv has {n_barcodes} rows, but matrix has {n_cells_needed} cells. Please check your data files.'
                }), 400
            elif n_barcodes > n_cells_needed:
                # 矩阵可能是过滤后的，只取前 n_cells_needed 行
                print(f"[UPLOAD] Barcodes has {n_barcodes} rows, matrix has {n_cells_needed} cells. Using first {n_cells_needed} rows.")
                barcodes_df = barcodes_df.iloc[:n_cells_needed].reset_index(drop=True)
                n_barcodes = len(barcodes_df)
            
            # 读取第一列
            cell_names = barcodes_df.iloc[:, 0].values.astype(str)
            
            # 检查内容是否为空/无效，如果是则按顺序补齐（0, 1, 2, ...）
            # 检查：全部为空、全部是NaN、全部是空字符串、或全部相同（可能是占位符）
            is_empty = len(cell_names) == 0
            is_all_na = all(pd.isna(cell_names)) if len(cell_names) > 0 else True
            is_all_empty_str = all(cell_names == '') if len(cell_names) > 0 else True
            is_all_same = len(set(cell_names)) == 1 if len(cell_names) > 0 else True
            
            if is_empty or is_all_na or is_all_empty_str or (is_all_same and len(cell_names) > 0 and (pd.isna(cell_names[0]) or cell_names[0] == '')):
                print(f"[UPLOAD] Barcodes file appears empty or invalid, using sequential labels (0-{n_cells_needed-1})")
                cell_names = np.array([str(i) for i in range(n_cells_needed)])
            else:
                # 检查是否有空值，如果有则用索引补齐
                mask = pd.isna(cell_names) | (cell_names == '')
                if mask.any():
                    print(f"[UPLOAD] Found {mask.sum()} empty values in barcodes.tsv, filling with sequential indices")
                    cell_names = np.where(mask, [str(i) for i in range(n_cells_needed)], cell_names)
            
            adata.obs_names = cell_names
            
            # 4. 确保基因名唯一（scanpy要求）
            if len(adata.var_names) != len(set(adata.var_names)):
                print(f"[UPLOAD] Making gene names unique...")
                seen = {}
                unique_names = []
                for name in adata.var_names:
                    if name in seen:
                        seen[name] += 1
                        unique_names.append(f"{name}-{seen[name]}")
                    else:
                        seen[name] = 0
                        unique_names.append(name)
                adata.var_names = unique_names
            
            print(f"[UPLOAD] Loaded: {adata.shape[0]} cells, {adata.shape[1]} genes")
            
            # Save to session
            save_adata(adata, session_id)
            
            result = {
                'success': True,
                'session_id': session_id,
                'shape': list(adata.shape),
                'obs_columns': list(adata.obs.columns),
                'obsm_keys': list(adata.obsm.keys()),
                'message': f'Loaded {adata.shape[0]} cells and {adata.shape[1]} genes from MTX files'
            }
            
            # Add warning if Seurat format was detected
            if format_warning:
                result['warning'] = format_warning
            
            return jsonify(result)
            
        finally:
            # Clean up temporary directory AFTER reading is complete
            if temp_dir.exists():
                try:
                    shutil.rmtree(temp_dir)
                    print(f"[UPLOAD] Cleaned up temp directory: {temp_dir}")
                except Exception as e:
                    print(f"[UPLOAD] Warning: Could not clean up temp directory: {e}")
    
    except Exception as e:
        print(f"[UPLOAD] Error in handle_mtx_multi_upload: {e}")
        import traceback
        traceback.print_exc()
        return jsonify({'error': f'Failed to process MTX files: {str(e)}'}), 500

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
        
        file_type = request.form.get('type', 'h5ad')
        
        # Handle multiple MTX files upload
        if file_type == 'mtx_multi':
            return handle_mtx_multi_upload()
        
        # Handle single file upload
        if 'file' not in request.files:
            print("[UPLOAD] Error: No file in request")
            return jsonify({'error': 'No file provided'}), 400
        
        file = request.files['file']
        
        print(f"[UPLOAD] File: {file.filename}, Type: {file_type}")
        
        if file.filename == '':
            print("[UPLOAD] Error: Empty filename")
            return jsonify({'error': 'No file selected'}), 400
        
        # Validate file type (h5ad, csv, tsv allowed; mtx should use multi-file upload)
        allowed_types = ['h5ad', 'csv', 'tsv']
        if file_type not in allowed_types:
            print(f"[UPLOAD] Error: Unsupported file type: {file_type}")
            return jsonify({'error': f'Unsupported file format. Only .h5ad, .csv, and .tsv files are allowed. For MTX format, please use the separate file upload option.'}), 400
        
        # Additional validation: check file extension
        file_ext = Path(file.filename).suffix.lower()
        allowed_extensions = ['.h5ad', '.csv', '.tsv']
        if file_ext not in allowed_extensions:
            print(f"[UPLOAD] Error: Invalid file extension: {file_ext}")
            return jsonify({'error': f'Invalid file extension. Only .h5ad, .csv, and .tsv files are allowed. For MTX format, please use the separate file upload option.'}), 400
        
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
                return jsonify({'error': f'Unsupported file type: {file_type}. For MTX format, please use the separate file upload option.'}), 400
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

