"""
Visualize API endpoints
"""

import numpy as np
from flask import Blueprint, request, jsonify
from ..config import LOTUS_AVAILABLE, SCANPY_AVAILABLE, lt, sc
from ..utils import load_adata, save_adata

bp = Blueprint('visualize', __name__, url_prefix='/api')


@bp.route('/visualize', methods=['POST'])
def run_visualization():
    """Run visualization (UMAP)"""
    try:
        if not LOTUS_AVAILABLE and not SCANPY_AVAILABLE:
            return jsonify({'error': 'Neither Lotus nor scanpy available'}), 500
        
        data = request.json
        session_id = data.get('session_id', 'default')
        cluster_key = data.get('cluster_key', 'cplearn')
        n_components = data.get('n_components', 2)
        min_dist = data.get('min_dist', 0.5)
        spread = data.get('spread', 1.0)
        use_rep = data.get('use_rep', None)  # Optional: use specific embedding from obsm
        
        adata = load_adata(session_id)
        
        if adata is None:
            return jsonify({'error': 'No data found.'}), 400
        
        # If use_rep is specified and exists, use it directly
        if use_rep and use_rep in adata.obsm:
            print(f"[VISUALIZE] Using existing embedding: {use_rep}")
            embedding = adata.obsm[use_rep]
            # Handle NaN values - filter them out or use valid points
            valid_mask = ~np.isnan(embedding).any(axis=1)
            if embedding.shape[1] >= n_components:
                # Use first n_components dimensions
                coords = embedding[:, :n_components].tolist()
            else:
                # Pad or use available dimensions
                coords = embedding.tolist()
        else:
            # Compute UMAP if not exists - Always use Lotus if available
            if 'X_umap' not in adata.obsm or adata.obsm['X_umap'].shape[1] != n_components:
                if LOTUS_AVAILABLE:
                    print(f"[VISUALIZE] Using Lotus UMAP computation")
                    lt.tl.umap(
                        adata,
                        n_components=n_components,
                        min_dist=min_dist,
                        spread=spread
                    )
                else:
                    if sc is None:
                        return jsonify({'error': 'scanpy not available for UMAP computation'}), 500
                    print(f"[VISUALIZE] Using scanpy UMAP (Lotus not available)")
                    sc.tl.umap(
                        adata,
                        n_components=n_components,
                        min_dist=min_dist,
                        spread=spread
                    )
            
            # Save updated AnnData
            save_adata(adata, session_id)
            
            # Extract coordinates
            coords = adata.obsm['X_umap'].tolist()
        
        # Extract metadata
        metadata = {}
        if cluster_key in adata.obs:
            metadata[cluster_key] = adata.obs[cluster_key].astype(str).tolist()
        
        # Add other obs columns
        for col in adata.obs.columns:
            if col != cluster_key:
                metadata[col] = adata.obs[col].astype(str).tolist()
        
        return jsonify({
            'success': True,
            'coordinates': coords,
            'metadata': metadata,
            'n_points': len(coords),
            'n_components': n_components
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

