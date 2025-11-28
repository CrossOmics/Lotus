"""
Visualize API endpoints
"""

import numpy as np
import hashlib
import json as json_lib
from typing import Optional
from flask import Blueprint, request, jsonify
from pathlib import Path
from ..config import LOTUS_AVAILABLE, SCANPY_AVAILABLE, lt, sc
from ..utils import load_adata, save_adata, get_session_dir

bp = Blueprint('visualize', __name__, url_prefix='/api')

# Cache for coremap visualizations: {cache_key: html_content}
_coremap_cache = {}


@bp.route('/visualize', methods=['POST'])
def run_visualization():
    """Run visualization (UMAP or Coremap)"""
    try:
        if not LOTUS_AVAILABLE and not SCANPY_AVAILABLE:
            return jsonify({'error': 'Neither Lotus nor scanpy available'}), 500
        
        # Handle both JSON and FormData requests
        if request.content_type and 'multipart/form-data' in request.content_type:
            # FormData request (for file uploads)
            data = request.form
            session_id = data.get('session_id', 'default')
            cluster_key = data.get('cluster_key', 'cplearn')
            n_components = int(data.get('n_components', 2))
            min_dist = float(data.get('min_dist', 0.5))
            spread = float(data.get('spread', 1.0))
            use_rep = data.get('use_rep', None)
            method = data.get('method', 'umap')
            
            # Handle ground truth file if provided
            truth_json_content = None
            if 'truth_json_file' in request.files:
                json_file = request.files['truth_json_file']
                if json_file.filename:
                    try:
                        truth_json_content = json_file.read().decode('utf-8')
                    except Exception as e:
                        return jsonify({'error': f'Failed to read JSON file: {str(e)}'}), 400
            
            truth_key = data.get('truth_key', None)
        else:
            # JSON request
            if not request.is_json:
                return jsonify({'error': 'Request must be JSON or multipart/form-data'}), 400
            
            data = request.json
            session_id = data.get('session_id', 'default')
            cluster_key = data.get('cluster_key', 'cplearn')
            n_components = data.get('n_components', 2)
            min_dist = data.get('min_dist', 0.5)
            spread = data.get('spread', 1.0)
            use_rep = data.get('use_rep', None)
            method = data.get('method', 'umap')
            truth_key = data.get('truth_key', None)
            truth_json_content = None
        
        adata = load_adata(session_id)
        
        if adata is None:
            return jsonify({'error': 'No data found.'}), 400
        
        # Handle coremap visualization
        if method == 'coremap' or (use_rep and 'coremap' in use_rep.lower()):
            # Get fast_view parameter
            if request.content_type and 'multipart/form-data' in request.content_type:
                fast_view_str = data.get('fast_view', 'true')
                fast_view = fast_view_str.lower() == 'true' if isinstance(fast_view_str, str) else bool(fast_view_str)
            else:
                fast_view = data.get('fast_view', True)
            
            return _generate_coremap_visualization_api(adata, session_id, cluster_key, n_components, use_rep, truth_key, truth_json_content, fast_view)
        
        # Handle UMAP visualization (existing logic)
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
        import traceback
        print(f"[VISUALIZE] Error: {str(e)}")
        print(traceback.format_exc())
        return jsonify({'error': str(e)}), 500


def _get_coremap_cache_key(session_id, cluster_key, n_components, use_rep, truth_key=None, truth_json_content=None, fast_view=True):
    """Generate a cache key for coremap visualization based on parameters"""
    # Create a unique key based on all parameters
    key_parts = [
        session_id,
        cluster_key or 'cplearn',
        str(n_components),
        use_rep or 'X_cplearn_coremap',
        truth_key or 'none',
        'fastview' if fast_view else 'full'
    ]
    
    # Add hash of truth_json_content if provided
    if truth_json_content:
        content_hash = hashlib.md5(truth_json_content.encode('utf-8')).hexdigest()[:8]
        key_parts.append(f'truth_json_{content_hash}')
    
    cache_key = '_'.join(key_parts)
    return cache_key


def _generate_coremap_visualization_api(adata, session_id, cluster_key, n_components, use_rep, truth_key=None, truth_json_content=None, fast_view=True):
    """Generate coremap visualization and return HTML content"""
    try:
        if not LOTUS_AVAILABLE:
            return jsonify({'error': 'Lotus not available for coremap visualization'}), 500
        
        # Find coremap key
        coremap_key = use_rep if use_rep and use_rep in adata.obsm else 'X_cplearn_coremap'
        if coremap_key not in adata.obsm:
            return jsonify({'error': f'Coremap embedding "{coremap_key}" not found. Please run Core Analysis first.'}), 400
        
        # Handle ground truth labels if provided (needed for cache key generation)
        truth_labels = None
        final_truth_key = truth_key
        if truth_key:
            # Reference existing ground truth in adata.obs
            if truth_key not in adata.obs:
                return jsonify({
                    'error': f'Ground truth key "{truth_key}" not found in adata.obs',
                    'available_keys': list(adata.obs.columns)
                }), 404
            final_truth_key = truth_key
            print(f"[VISUALIZE] Using existing ground truth from adata.obs['{truth_key}']")
        elif truth_json_content:
            # Parse ground truth from JSON
            try:
                truth_data = json_lib.loads(truth_json_content)
                if isinstance(truth_data, list):
                    truth_labels = truth_data
                elif isinstance(truth_data, dict):
                    if 'key' in truth_data:
                        final_truth_key = truth_data['key']
                        print(f"[VISUALIZE] Using truth_key from JSON: {final_truth_key}")
                    if 'labels' in truth_data:
                        truth_labels = truth_data['labels']
                    elif 'ground_truth' in truth_data:
                        truth_labels = truth_data['ground_truth']
                    else:
                        truth_labels = list(truth_data.values())
                print(f"[VISUALIZE] Loaded {len(truth_labels)} ground truth labels from JSON")
            except json_lib.JSONDecodeError as e:
                return jsonify({'error': f'Invalid JSON format: {str(e)}'}), 400
        
        # Generate cache key after processing truth_key (may have changed from JSON)
        cache_key = _get_coremap_cache_key(session_id, cluster_key, n_components, coremap_key, final_truth_key, truth_json_content, fast_view)
        
        # Check cache first (before expensive operations like model recreation)
        if cache_key in _coremap_cache:
            print(f"[VISUALIZE] Using cached coremap visualization (key: {cache_key[:50]}...)")
            cached_html = _coremap_cache[cache_key]
            
            # Also check if file exists on disk (for persistence)
            session_dir = get_session_dir(session_id)
            viz_path = session_dir / f'coremap_{coremap_key}_{cache_key}.html'
            if viz_path.exists():
                return jsonify({
                    'success': True,
                    'html_content': cached_html,
                    'file_url': f'/api/files/{session_id}/{viz_path.name}',
                    'method': 'coremap',
                    'cached': True
                })
            else:
                # Cache exists but file doesn't, regenerate
                print(f"[VISUALIZE] Cache exists but file missing, regenerating...")
                del _coremap_cache[cache_key]
        
        # Get model from adata.uns if available - need to recreate it
        model_info_key = f'{cluster_key}_model_info'
        model = None
        if model_info_key in adata.uns:
            # Recreate model by re-running clustering
            from ..config import cluster
            model_info = adata.uns[model_info_key]
            rep_to_use = model_info.get('use_rep') or 'X_pca'
            resolution = model_info.get('resolution', 1.2)
            
            print(f"[VISUALIZE] Recreating model for coremap visualization (resolution={resolution})...")
            model = cluster(
                adata,
                method='cplearn',
                use_rep=rep_to_use,
                key_added=cluster_key,
                cluster_resolution=resolution,
                print_summary=False
            )
        
        # Handle ground truth labels if provided
        if truth_labels is not None:
            if len(truth_labels) != adata.n_obs:
                return jsonify({
                    'error': f'Truth labels length ({len(truth_labels)}) does not match number of cells ({adata.n_obs})',
                    'n_cells': adata.n_obs,
                    'n_labels': len(truth_labels)
                }), 400
            
            # Save truth labels to adata.obs
            if final_truth_key is None:
                final_truth_key = 'ground_truth'
            
            import pandas as pd
            adata.obs[final_truth_key] = pd.Categorical(truth_labels)
            print(f"[VISUALIZE] Saved ground truth labels to adata.obs['{final_truth_key}']")
            
            # Save updated AnnData
            save_adata(adata, session_id)
        
        # Generate visualization using the local function
        viz_result = _generate_coremap_visualization(
            adata=adata,
            model=model,
            coremap_key=coremap_key,
            cluster_key=cluster_key,
            truth_key=final_truth_key,
            fast_view=fast_view,
            use_webgl=True,
            save_viz=True,
            session_id=session_id,
            use_rep=None,
            cache_key=cache_key  # Pass cache key for file naming
        )
        
        if viz_result and 'file_url' in viz_result:
            # Read the HTML file and return it
            session_dir = get_session_dir(session_id)
            # Use cache key in filename for better organization
            viz_path = session_dir / f'coremap_{coremap_key}_{cache_key}.html'
            
            # Fallback to old filename if new one doesn't exist
            if not viz_path.exists():
                old_viz_path = session_dir / f'coremap_{coremap_key}.html'
                if old_viz_path.exists():
                    viz_path = old_viz_path
            
            if viz_path.exists():
                with open(viz_path, 'r', encoding='utf-8') as f:
                    html_content = f.read()
                
                # Store in cache
                _coremap_cache[cache_key] = html_content
                print(f"[VISUALIZE] Cached coremap visualization (key: {cache_key[:50]}...)")
                
                return jsonify({
                    'success': True,
                    'html_content': html_content,
                    'file_url': viz_result['file_url'],
                    'method': 'coremap',
                    'cached': False
                })
            else:
                return jsonify({'error': 'Visualization file not found after generation'}), 500
        else:
            return jsonify({'error': 'Failed to generate coremap visualization'}), 500
            
    except Exception as e:
        import traceback
        print(f"[VISUALIZE] Coremap visualization error: {str(e)}")
        print(traceback.format_exc())
        return jsonify({'error': f'Coremap visualization failed: {str(e)}'}), 500


def _generate_coremap_visualization(
    adata,
    model,
    coremap_key: str,
    cluster_key: str,
    truth_key: Optional[str],
    fast_view: bool,
    use_webgl: bool,
    save_viz: bool,
    session_id: str,
    use_rep: Optional[str],
    cache_key: Optional[str] = None
):
    """Generate coremap visualization using cplearn's visualize_coremap
    
    This function is now part of the Visualization module, completely decoupled from Core Analysis.
    """
    try:
        import umap
        from cplearn.coremap import Coremap
        from cplearn.coremap.vizualizer import visualize_coremap
        from pathlib import Path
        from ..utils import get_session_dir
        import numpy as np
    except ImportError as e:
        return {'error': f'Required packages not available: {str(e)}'}
    
    # Step 1: Generate UMAP skeleton if not exists
    if 'X_umap' not in adata.obsm or adata.obsm['X_umap'].shape[1] != 2:
        print("[VISUALIZE] Computing UMAP skeleton...")
        # Use the representation for UMAP
        if use_rep and use_rep in adata.obsm:
            X_for_umap = adata.obsm[use_rep]
        elif 'X_pca' in adata.obsm:
            X_for_umap = adata.obsm['X_pca']
        elif 'X_latent' in adata.obsm:
            X_for_umap = adata.obsm['X_latent']
        else:
            X_for_umap = adata.X
            if hasattr(X_for_umap, 'toarray'):
                X_for_umap = X_for_umap.toarray()
        
        reducer = umap.UMAP(n_components=2)
        X_umap = reducer.fit_transform(X_for_umap)
        adata.obsm['X_umap_skeleton'] = X_umap
        print(f"[VISUALIZE] UMAP skeleton computed: shape {X_umap.shape}")
    else:
        X_umap = adata.obsm['X_umap']
        print("[VISUALIZE] Using existing UMAP embedding")
    
    # Step 2: Initialize Coremap module
    print("[VISUALIZE] Initializing Coremap...")
    cmap = Coremap(model, global_umap=X_umap, fast_view=fast_view)
    
    # Step 3: Prepare labels for visualization
    # Use truth labels if provided, otherwise use cluster labels
    labels_to_use = None
    label_source = None
    
    if truth_key and truth_key in adata.obs:
        labels_to_use = np.asarray(adata.obs[truth_key].values)
        label_source = truth_key
        print(f"[VISUALIZE] Using ground truth labels from '{truth_key}'")
    elif cluster_key in adata.obs:
        labels_to_use = np.asarray(adata.obs[cluster_key].values, dtype=int)
        label_source = cluster_key
        print(f"[VISUALIZE] Using cluster labels from '{cluster_key}'")
    elif hasattr(model, 'labels_') and model is not None:
        labels_to_use = np.asarray(model.labels_, dtype=int)
        label_source = 'model.labels_'
        print("[VISUALIZE] Using labels from model")
    else:
        print("[VISUALIZE] Warning: No labels found, visualization may not show clusters")
    
    # Step 4: Generate visualization
    print("[VISUALIZE] Generating layer-wise visualization...")
    fig = visualize_coremap(cmap, labels_to_use, use_webgl=use_webgl)
    
    # Step 5: Save visualization if requested
    viz_path = None
    if save_viz:
        session_dir = get_session_dir(session_id)
        session_dir.mkdir(exist_ok=True, parents=True)
        
        # Use cache key in filename if provided, otherwise use default
        if cache_key:
            viz_path = session_dir / f'coremap_{coremap_key}_{cache_key}.html'
        else:
            viz_path = session_dir / f'coremap_{coremap_key}.html'
        
        fig.write_html(str(viz_path))
        print(f"[VISUALIZE] Visualization saved to: {viz_path}")
    
    # Return visualization info
    result = {
        'success': True,
        'label_source': label_source,
        'fast_view': fast_view,
        'use_webgl': use_webgl
    }
    
    if viz_path:
        result['file_path'] = str(viz_path)
        result['file_url'] = f'/api/files/{session_id}/{viz_path.name}'
    
    return result

