"""
Preprocess API endpoints
"""

from flask import Blueprint, request, jsonify
from ..config import LOTUS_AVAILABLE, SCANPY_AVAILABLE, preprocess, sc
from ..utils import load_adata, save_adata

bp = Blueprint('preprocess', __name__, url_prefix='/api')


@bp.route('/preprocess', methods=['POST'])
def run_preprocess():
    """Run preprocessing pipeline"""
    try:
        if not LOTUS_AVAILABLE and not SCANPY_AVAILABLE:
            return jsonify({'error': 'Neither Lotus nor scanpy available'}), 500
        
        data = request.json
        session_id = data.get('session_id', 'default')
        n_pcs = data.get('n_pcs', 20)
        n_neighbors = data.get('n_neighbors', 15)
        
        adata = load_adata(session_id)
        
        if adata is None:
            return jsonify({'error': 'No data found. Please upload data first.'}), 400
        
        # Run preprocessing
        if LOTUS_AVAILABLE:
            preprocess(
                adata,
                n_pcs=n_pcs,
                n_neighbors=n_neighbors,
                save_raw=True
            )
        else:
            # Use scanpy for preprocessing
            sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata, n_top_genes=2000)
            sc.pp.scale(adata, max_value=10)
            sc.tl.pca(adata, n_comps=n_pcs)
            sc.pp.neighbors(adata, n_neighbors=n_neighbors)
        
        # Save updated AnnData
        save_adata(adata, session_id)
        
        return jsonify({
            'success': True,
            'message': 'Preprocessing complete',
            'obsm_keys': list(adata.obsm.keys()),
            'has_pca': 'X_pca' in adata.obsm,
            'has_neighbors': 'neighbors' in adata.uns
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

