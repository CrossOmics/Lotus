"""
Info API endpoints (health, metadata, cluster-info, file serving)
"""

from flask import Blueprint, request, jsonify, send_from_directory
from pathlib import Path
from ..utils import load_adata, get_session_dir
from ..config import LOTUS_AVAILABLE

bp = Blueprint('get_info', __name__, url_prefix='/api')


@bp.route('/health')
def health():
    """Health check endpoint"""
    return jsonify({
        'status': 'ok',
        'lotus_available': LOTUS_AVAILABLE
    })


@bp.route('/metadata', methods=['GET'])
def get_metadata():
    """Get available metadata columns"""
    try:
        session_id = request.args.get('session_id', 'default')
        adata = load_adata(session_id)
        
        if adata is None:
            return jsonify({'error': 'No data found.'}), 400
        
        return jsonify({
            'obs_columns': list(adata.obs.columns),
            'obsm_keys': list(adata.obsm.keys()),
            'var_columns': list(adata.var.columns) if hasattr(adata, 'var') else []
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@bp.route('/cluster-info', methods=['POST'])
def get_cluster_info():
    """Get information about a specific cluster"""
    try:
        import numpy as np
        
        data = request.json
        session_id = data.get('session_id', 'default')
        cluster_id = data.get('cluster_id')
        cluster_key = data.get('cluster_key', 'cplearn')
        
        adata = load_adata(session_id)
        
        if adata is None:
            return jsonify({'error': 'No data found.'}), 400
        
        if cluster_key not in adata.obs:
            return jsonify({'error': f'Cluster key {cluster_key} not found'}), 400
        
        # Filter cells in this cluster
        cluster_mask = adata.obs[cluster_key].astype(str) == str(cluster_id)
        cluster_cells = adata[cluster_mask]
        
        # Get cell type information (if available)
        cell_type_cols = [col for col in cluster_cells.obs.columns 
                         if 'type' in col.lower() or 'cell_type' in col.lower()]
        
        cell_types = {}
        if cell_type_cols:
            for col in cell_type_cols:
                cell_types[col] = cluster_cells.obs[col].value_counts().to_dict()
        
        # Get statistics
        stats = {
            'n_cells': int(cluster_mask.sum()),
            'percentage': float(cluster_mask.sum() / len(adata) * 100)
        }
        
        # Get top expressed genes (if expression data available)
        top_genes = None
        if hasattr(cluster_cells, 'X') and cluster_cells.X is not None:
            try:
                # Average expression per gene in this cluster
                avg_expr = np.array(cluster_cells.X.mean(axis=0)).flatten()
                top_gene_indices = np.argsort(avg_expr)[-10:][::-1]
                top_genes = [cluster_cells.var_names[i] for i in top_gene_indices]
            except:
                pass
        
        return jsonify({
            'success': True,
            'cluster_id': str(cluster_id),
            'stats': stats,
            'cell_types': cell_types,
            'top_genes': top_genes
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@bp.route('/files/<session_id>/<filename>')
def serve_file(session_id: str, filename: str):
    """Serve files from session directory (e.g., visualization HTML files)"""
    try:
        session_dir = get_session_dir(session_id)
        file_path = session_dir / filename
        
        if not file_path.exists():
            return jsonify({'error': 'File not found'}), 404
        
        # Security check: ensure file is within session directory
        try:
            file_path.resolve().relative_to(session_dir.resolve())
        except ValueError:
            return jsonify({'error': 'Invalid file path'}), 403
        
        return send_from_directory(str(session_dir), filename)
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

