"""
Core Selection API endpoints
"""

from typing import Optional
import numpy as np
from flask import Blueprint, request, jsonify
from ..config import LOTUS_AVAILABLE, clustering
from ..utils import load_adata, save_adata

bp = Blueprint('core_selection', __name__, url_prefix='/api')


@bp.route('/core-selection', methods=['POST'])
def run_core_selection():
    """Run core selection (core map embedding) with visualization
    
    This endpoint can run independently of clustering:
    - If clustering already exists, it will use the existing clustering
    - If no clustering exists, it will automatically run clustering first
    - Users can optionally provide clustering parameters
    - Always generates visualization
    - Supports ground truth labels via JSON file upload
    
    Parameters:
    - cluster_resolution: Clustering resolution (default: 1.2)
    - cluster_key: Key name for cluster labels (default: 'cplearn')
    - truth_key: Key name for ground truth labels in adata.obs (optional)
    - truth_json: JSON file with ground truth labels (optional)
    """
    try:
        if not LOTUS_AVAILABLE:
            return jsonify({'error': 'Lotus not available'}), 500
        
        from lotus.workflows import core_analysis
        from lotus.workflows.preprocess import neighbors as compute_neighbors
        from pathlib import Path
        import json as json_lib
        from ..utils import get_session_dir
        
        # Handle both JSON and form data (for file upload)
        truth_key = None  # Will be extracted from JSON if present
        if request.is_json:
            data = request.json
            session_id = data.get('session_id', 'default')
            use_rep = data.get('use_rep', None)
            key_added = data.get('key_added', 'X_cplearn_coremap')
            cluster_resolution = data.get('cluster_resolution', 1.2)
            cluster_key = data.get('cluster_key', 'cplearn')
            fast_view = data.get('fast_view', True)
            truth_json_content = data.get('truth_json', None)  # JSON content as string
        else:
            # Form data (for file upload)
            data = request.form
            session_id = data.get('session_id', 'default')
            use_rep = data.get('use_rep', None)
            if use_rep == '':  # Empty string means None
                use_rep = None
            key_added = data.get('key_added', 'X_cplearn_coremap')
            cluster_resolution = float(data.get('cluster_resolution', 1.2))
            cluster_key = data.get('cluster_key', 'cplearn')
            fast_view_str = data.get('fast_view', 'true')
            fast_view = fast_view_str.lower() == 'true' if isinstance(fast_view_str, str) else bool(fast_view_str)
            truth_json_content = None
            
            # Handle JSON file upload
            if 'truth_json_file' in request.files:
                json_file = request.files['truth_json_file']
                if json_file.filename:
                    try:
                        truth_json_content = json_file.read().decode('utf-8')
                    except Exception as e:
                        return jsonify({'error': f'Failed to read JSON file: {str(e)}'}), 400
        
        # Default visualization parameters (always enabled)
        visualize = True  # Always generate visualization
        use_webgl = True
        save_viz = True
        
        # Parse ground truth labels from JSON
        truth_labels = None
        if truth_json_content:
            try:
                truth_data = json_lib.loads(truth_json_content)
                # Handle different JSON formats
                if isinstance(truth_data, list):
                    # Direct list of labels
                    truth_labels = truth_data
                elif isinstance(truth_data, dict):
                    # Extract key from JSON if present
                    if 'key' in truth_data:
                        truth_key = truth_data['key']
                        print(f"[CORE] Using truth_key from JSON: {truth_key}")
                    
                    # Extract labels
                    if 'labels' in truth_data:
                        truth_labels = truth_data['labels']
                    elif 'ground_truth' in truth_data:
                        truth_labels = truth_data['ground_truth']
                    else:
                        # Assume values are labels, keys are cell IDs (we'll use values)
                        truth_labels = list(truth_data.values())
                else:
                    return jsonify({'error': 'Invalid JSON format for ground truth labels'}), 400
                
                print(f"[CORE] Loaded {len(truth_labels)} ground truth labels from JSON")
                if truth_key:
                    print(f"[CORE] Using key '{truth_key}' from JSON for adata.obs")
            except json_lib.JSONDecodeError as e:
                return jsonify({'error': f'Invalid JSON format: {str(e)}'}), 400
        
        adata = load_adata(session_id)
        
        if adata is None:
            return jsonify({'error': 'No data found.'}), 400
        
        # Auto-detect representation if not specified
        if use_rep is None:
            if 'X_latent' in adata.obsm:
                use_rep = 'X_latent'
            elif 'X_pca' in adata.obsm:
                use_rep = 'X_pca'
            else:
                # Check if we have any representation
                available_reps = list(adata.obsm.keys())
                if available_reps:
                    use_rep = available_reps[0]
                    print(f"[CORE] Warning: Using first available representation: {use_rep}")
                else:
                    return jsonify({
                        'error': 'No representation found. Please run preprocessing first to generate PCA or latent representation.',
                        'suggestion': 'Run /api/preprocess endpoint first.'
                    }), 400
        
        # Check if we need to run clustering
        needs_clustering = False
        model = None
        
        # Check if clustering already exists
        if cluster_key in adata.obs:
            print(f"[CORE] Found existing clustering in adata.obs['{cluster_key}']")
            # Check if we have model info to recreate the model
            model_info_key = f'{cluster_key}_model_info'
            if model_info_key in adata.uns:
                print("[CORE] Found model info, recreating model...")
                model_info = adata.uns[model_info_key]
                rep_to_use = use_rep or model_info.get('use_rep') or use_rep
                resolution = model_info.get('resolution', cluster_resolution)
                
                # Ensure neighbors graph exists - get n_neighbors from existing neighbors if available
                if 'neighbors' not in adata.uns:
                    print("[CORE] Neighbors graph not found. Please run preprocessing first.")
                    return jsonify({
                        'error': 'Neighbors graph not found. Please run preprocessing first.',
                        'suggestion': 'Run /api/preprocess endpoint to compute neighbors graph.'
                    }), 400
                else:
                    # Get n_neighbors from existing neighbors graph
                    neighbors_params = adata.uns['neighbors']['params']
                    n_neighbors = neighbors_params.get('n_neighbors', 15)
                    print(f"[CORE] Using existing neighbors graph (n_neighbors={n_neighbors})")
                
                # Re-run clustering to get model
                print(f"[CORE] Re-running clustering to obtain model (resolution={resolution})...")
                # Use same parameters as lotus_workflow.py for consistency
                model = clustering(
                    adata,
                    method='cplearn',
                    use_rep=rep_to_use,
                    key_added=cluster_key,
                    cluster_resolution=resolution,
                    stable_core_frac=0.25,  # Match lotus_workflow.py
                    stable_ng_num=8,  # Match lotus_workflow.py
                    fine_grained=False,  # Match lotus_workflow.py
                    propagate=True,  # Match lotus_workflow.py
                    print_summary=False  # Don't print in API
                )
            else:
                # No model info, need to run clustering
                needs_clustering = True
        else:
            # No clustering exists, need to run it
            needs_clustering = True
        
        # Run clustering if needed
        if needs_clustering:
            print(f"[CORE] No existing clustering found. Running clustering first...")
            
            # Ensure neighbors graph exists - get from preprocessing results
            if 'neighbors' not in adata.uns:
                print("[CORE] Neighbors graph not found. Please run preprocessing first.")
                return jsonify({
                    'error': 'Neighbors graph not found. Please run preprocessing first.',
                    'suggestion': 'Run /api/preprocess endpoint to compute neighbors graph.'
                }), 400
            else:
                # Get n_neighbors from existing neighbors graph
                neighbors_params = adata.uns['neighbors']['params']
                n_neighbors = neighbors_params.get('n_neighbors', 15)
                print(f"[CORE] Using existing neighbors graph from preprocessing (n_neighbors={n_neighbors})")
                print(f"[CORE] Parameters: resolution={cluster_resolution}, use_rep={use_rep}")
            
            # Run clustering to get model (use cplearn method)
            print(f"[CORE] Running clustering with resolution={cluster_resolution}...")
            # Use same parameters as lotus_workflow.py for consistency
            model = clustering(
                adata,
                method='cplearn',
                use_rep=use_rep,
                key_added=cluster_key,
                cluster_resolution=cluster_resolution,
                stable_core_frac=0.25,  # Match lotus_workflow.py
                stable_ng_num=8,  # Match lotus_workflow.py
                fine_grained=False,  # Match lotus_workflow.py
                propagate=True,  # Match lotus_workflow.py
                print_summary=False  # Don't print in API
            )
            
            # Save model info for future use
            adata.uns[f'{cluster_key}_model_info'] = {
                'method': 'cplearn',
                'resolution': cluster_resolution,
                'use_rep': use_rep
            }
            save_adata(adata, session_id)
        
        # Now run core analysis
        print(f"[CORE] Running core analysis with use_rep={use_rep}, key_added={key_added}")
        core_analysis(
            adata,
            model=model,
            use_rep=use_rep,
            key_added=key_added,
            cluster_key=cluster_key
        )
        
        # Handle ground truth labels if provided
        if truth_labels is not None:
            if len(truth_labels) != adata.n_obs:
                return jsonify({
                    'error': f'Truth labels length ({len(truth_labels)}) does not match number of cells ({adata.n_obs})'
                }), 400
            
            # Save truth labels to adata.obs
            if truth_key is None:
                truth_key = 'ground_truth'
            
            import pandas as pd
            adata.obs[truth_key] = pd.Categorical(truth_labels)
            print(f"[CORE] Saved ground truth labels to adata.obs['{truth_key}']")
        
        # Save updated AnnData
        save_adata(adata, session_id)
        
        # Get assignment statistics
        embedding = adata.obsm[key_added]
        assigned = np.sum(~np.isnan(embedding).any(axis=1))
        n_core = 0
        if f"{key_added}_is_core" in adata.obs:
            n_core = int(np.sum(adata.obs[f"{key_added}_is_core"] == True))
        
        # Generate visualization if requested
        viz_result = None
        if visualize:
            try:
                print("[CORE] Generating visualization...")
                viz_result = _generate_coremap_visualization(
                    adata=adata,
                    model=model,
                    coremap_key=key_added,
                    cluster_key=cluster_key,
                    truth_key=truth_key,
                    fast_view=fast_view,
                    use_webgl=use_webgl,
                    save_viz=save_viz,
                    session_id=session_id,
                    use_rep=use_rep
                )
                print("[CORE] Visualization generated successfully")
            except Exception as viz_error:
                import traceback
                print(f"[CORE] Visualization failed: {str(viz_error)}")
                print(traceback.format_exc())
                # Don't fail the whole request if visualization fails
                viz_result = {'error': str(viz_error), 'warning': 'Visualization failed but core selection completed'}
        
        response = {
            'success': True,
            'key_added': key_added,
            'obsm_keys': list(adata.obsm.keys()),
            'assigned_points': int(assigned),
            'total_points': int(adata.n_obs),
            'core_cells': n_core,
            'clustering_auto_run': needs_clustering,
            'cluster_key': cluster_key,
            'message': f'Core selection complete. Embedding stored in obsm["{key_added}"] ({assigned}/{adata.n_obs} points assigned, {n_core} core cells)'
        }
        
        if viz_result:
            response['visualization'] = viz_result
        
        if truth_key and truth_key in adata.obs:
            response['truth_key'] = truth_key
        
        return jsonify(response)
    
    except Exception as e:
        import traceback
        error_msg = f'Core selection failed: {str(e)}'
        print(f"[CORE] Error: {error_msg}")
        print(traceback.format_exc())
        return jsonify({'error': error_msg}), 500


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
    use_rep: Optional[str]
):
    """Generate coremap visualization using cplearn's visualize_coremap"""
    try:
        import umap
        from cplearn.coremap import Coremap
        from cplearn.coremap.vizualizer import visualize_coremap
        from pathlib import Path
        from ..utils import get_session_dir
    except ImportError as e:
        return {'error': f'Required packages not available: {str(e)}'}
    
    # Step 1: Generate UMAP skeleton if not exists
    if 'X_umap' not in adata.obsm or adata.obsm['X_umap'].shape[1] != 2:
        print("[CORE-VIZ] Computing UMAP skeleton...")
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
        
        reducer = umap.UMAP(n_components=2, random_state=42)
        X_umap = reducer.fit_transform(X_for_umap)
        adata.obsm['X_umap_skeleton'] = X_umap
        print(f"[CORE-VIZ] UMAP skeleton computed: shape {X_umap.shape}")
    else:
        X_umap = adata.obsm['X_umap']
        print("[CORE-VIZ] Using existing UMAP embedding")
    
    # Step 2: Initialize Coremap module
    print("[CORE-VIZ] Initializing Coremap...")
    cmap = Coremap(model, global_umap=X_umap, fast_view=fast_view)
    
    # Step 3: Prepare labels for visualization
    # Use truth labels if provided, otherwise use cluster labels
    labels_to_use = None
    label_source = None
    
    if truth_key and truth_key in adata.obs:
        labels_to_use = np.asarray(adata.obs[truth_key].values)
        label_source = truth_key
        print(f"[CORE-VIZ] Using ground truth labels from '{truth_key}'")
    elif cluster_key in adata.obs:
        labels_to_use = np.asarray(adata.obs[cluster_key].values, dtype=int)
        label_source = cluster_key
        print(f"[CORE-VIZ] Using cluster labels from '{cluster_key}'")
    elif hasattr(model, 'labels_'):
        labels_to_use = np.asarray(model.labels_, dtype=int)
        label_source = 'model.labels_'
        print("[CORE-VIZ] Using labels from model")
    else:
        print("[CORE-VIZ] Warning: No labels found, visualization may not show clusters")
    
    # Step 4: Generate visualization
    print("[CORE-VIZ] Generating layer-wise visualization...")
    fig = visualize_coremap(cmap, labels_to_use, use_webgl=use_webgl)
    
    # Step 5: Save visualization if requested
    viz_path = None
    if save_viz:
        session_dir = get_session_dir(session_id)
        viz_path = session_dir / f'coremap_{coremap_key}.html'
        session_dir.mkdir(exist_ok=True, parents=True)
        
        fig.write_html(str(viz_path))
        print(f"[CORE-VIZ] Visualization saved to: {viz_path}")
    
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

