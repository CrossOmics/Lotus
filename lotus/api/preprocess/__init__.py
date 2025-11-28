"""
Preprocess API endpoints
"""

import gc
import sys
from flask import Blueprint, request, jsonify
from ..config import LOTUS_AVAILABLE, SCANPY_AVAILABLE, preprocess, sc, DEFAULT_N_PCS, DEFAULT_N_NEIGHBORS, DEFAULT_N_TOP_GENES
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
        n_pcs = data.get('n_pcs', DEFAULT_N_PCS)  # Optimized default for 4GB memory
        n_neighbors = data.get('n_neighbors', DEFAULT_N_NEIGHBORS)  # Optimized default
        target_sum = data.get('target_sum', 1e4)
        # Use optimized default if not provided, but allow None for auto-detection
        n_top_genes = data.get('n_top_genes') if 'n_top_genes' in data else DEFAULT_N_TOP_GENES
        use_rep = data.get('use_rep', 'X_pca')
        save_raw = data.get('save_raw', True)
        raw_layer = data.get('raw_layer', 'raw_counts')
        min_genes = data.get('min_genes', None)
        min_cells = data.get('min_cells', None)
        # New filtering parameters
        min_counts = data.get('min_counts', None)
        max_counts = data.get('max_counts', None)
        max_genes = data.get('max_genes', None)
        pct_mt_max = data.get('pct_mt_max', None)
        # HVG flavor
        hvg_flavor = data.get('hvg_flavor', 'seurat_v3')
        # Batch correction parameters
        batch_key = data.get('batch_key', None)
        regress_out_keys = data.get('regress_out_keys', None)
        use_combat = data.get('use_combat', False)
        
        adata = load_adata(session_id)
        
        if adata is None:
            return jsonify({'error': 'No data found. Please upload data first.'}), 400
        
        # Run preprocessing
        if LOTUS_AVAILABLE:
            preprocess(
                adata,
                n_pcs=n_pcs,
                n_neighbors=n_neighbors,
                target_sum=target_sum,
                n_top_genes=n_top_genes,
                use_rep=use_rep,
                save_raw=save_raw,
                raw_layer=raw_layer,
                min_genes=min_genes,
                min_cells=min_cells,
                min_counts=min_counts,
                max_counts=max_counts,
                max_genes=max_genes,
                pct_mt_max=pct_mt_max,
                hvg_flavor=hvg_flavor,
                batch_key=batch_key,
                regress_out_keys=regress_out_keys,
                use_combat=use_combat,
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
        
        # Get results before freeing memory
        obsm_keys = list(adata.obsm.keys())
        has_pca = 'X_pca' in adata.obsm
        has_neighbors = 'neighbors' in adata.uns
        
        # Extract PCA variance ratio data for visualization
        pca_variance_ratio = None
        if 'pca' in adata.uns and 'variance_ratio' in adata.uns['pca']:
            variance_ratio = adata.uns['pca']['variance_ratio']
            # Convert numpy array to list for JSON serialization
            pca_variance_ratio = variance_ratio.tolist() if hasattr(variance_ratio, 'tolist') else list(variance_ratio)
        
        # Free memory after preprocessing
        del adata
        gc.collect()
        sys.stdout.flush()
        
        return jsonify({
            'success': True,
            'message': 'Preprocessing complete',
            'obsm_keys': obsm_keys,
            'has_pca': has_pca,
            'has_neighbors': has_neighbors,
            'pca_variance_ratio': pca_variance_ratio
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

