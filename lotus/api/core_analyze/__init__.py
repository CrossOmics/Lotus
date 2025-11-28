"""
Core Analysis API endpoints
"""

import gc
import sys
import traceback
from flask import Blueprint, request, jsonify
import numpy as np
from ..config import LOTUS_AVAILABLE
from ..utils import load_adata, save_adata

bp = Blueprint('core_analyze', __name__, url_prefix='/api')


@bp.route('/core-analyze', methods=['POST'])
def run_core_analysis():
    """Run core analysis (core map embedding) with visualization
    
    This endpoint can run independently of clustering:
    - If clustering already exists, it will use the existing clustering
    - If no clustering exists, it will automatically run clustering first
    - Users can optionally provide clustering parameters
    - Note: Visualization is now handled by Visualization module
    - Note: Ground Truth is now handled by Visualization module when using coremap
    
    Parameters:
    - cluster_resolution: Clustering resolution (default: 1.2)
    - cluster_key: Key name for cluster labels (default: 'cplearn')
    """
    import sys
    import traceback
    import gc
    
    try:
        if not LOTUS_AVAILABLE:
            return jsonify({'error': 'Lotus not available'}), 500
        
        from lotus.workflows import core_analyze
        from lotus.workflows.preprocessing import neighbors as compute_neighbors
        from lotus.workflows.clustering import cluster
        from pathlib import Path
        
        print(f"[CORE] Starting core analysis request...")
        sys.stdout.flush()  # Ensure logs are flushed immediately
        
        # Handle JSON request (Ground Truth is now handled by Visualization module)
        if not request.is_json:
            return jsonify({'error': 'Request must be JSON. Ground Truth should be handled via Visualization API.'}), 400
        
        data = request.json
        session_id = data.get('session_id', 'default')
        use_rep = data.get('use_rep', None)
        key_added = data.get('key_added', 'X_cplearn_coremap')
        cluster_resolution = data.get('cluster_resolution', 1.2)
        cluster_key = data.get('cluster_key', 'cplearn')
        # Note: fast_view is now handled by Visualization module
        # Note: Ground Truth is now handled by Visualization module
        # cplearn clustering parameters (used when auto-running clustering)
        stable_core_frac = data.get('stable_core_frac', 0.25)
        stable_ng_num = data.get('stable_ng_num', 8)
        fine_grained = data.get('fine_grained', False)
        propagate = data.get('propagate', True)
        
        # Load adata first to validate session
        adata = load_adata(session_id)
        
        if adata is None:
            return jsonify({'error': 'No data found.'}), 400
        
        # Note: Ground Truth is now handled by Visualization module when using coremap
        
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
                # Use parameters from request
                model = cluster(
                    adata,
                    method='cplearn',
                    use_rep=rep_to_use,
                    key_added=cluster_key,
                    cluster_resolution=resolution,
                    stable_core_frac=stable_core_frac,
                    stable_ng_num=stable_ng_num,
                    fine_grained=fine_grained,
                    propagate=propagate,
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
            # Use parameters from request
            model = cluster(
                adata,
                method='cplearn',
                use_rep=use_rep,
                key_added=cluster_key,
                cluster_resolution=cluster_resolution,
                stable_core_frac=stable_core_frac,
                stable_ng_num=stable_ng_num,
                fine_grained=fine_grained,
                propagate=propagate,
                print_summary=False  # Don't print in API
            )
            
            # Save model info for future use
            adata.uns[f'{cluster_key}_model_info'] = {
                'method': 'cplearn',
                'resolution': cluster_resolution,
                'use_rep': use_rep
            }
            # Save intermediate result (clustering) before core analysis
            print(f"[CORE] Saving intermediate result (clustering)...")
            save_adata(adata, session_id)
            sys.stdout.flush()
            gc.collect()  # Free memory after saving
        
        # Now run core analysis
        print(f"[CORE] Running core analysis with use_rep={use_rep}, key_added={key_added}")
        sys.stdout.flush()
        
        try:
            core_analyze(
                adata,
                model=model,
                use_rep=use_rep,
                key_added=key_added,
                cluster_key=cluster_key
            )
            print(f"[CORE] Core analysis completed successfully")
            sys.stdout.flush()
            
            # Save intermediate result (core analysis) before visualization
            print(f"[CORE] Saving intermediate result (core analysis)...")
            save_adata(adata, session_id)
            sys.stdout.flush()
            gc.collect()  # Free memory after saving
        except Exception as core_error:
            # Save what we have so far
            print(f"[CORE] Core analysis failed, but saving intermediate results...")
            try:
                save_adata(adata, session_id)
            except:
                pass
            raise core_error
        
        # Note: Ground Truth is now handled by Visualization module when using coremap
        
        # Get assignment statistics before freeing memory
        embedding = adata.obsm[key_added]
        assigned = int(np.sum(~np.isnan(embedding).any(axis=1)))
        n_core = 0
        if f"{key_added}_is_core" in adata.obs:
            n_core = int(np.sum(adata.obs[f"{key_added}_is_core"] == True))
        total_points = int(adata.n_obs)
        obsm_keys = list(adata.obsm.keys())
        
        # Build response before freeing memory
        # Note: Visualization is now handled by the Visualization module
        # Note: Ground Truth is now handled by the Visualization module
        response = {
            'success': True,
            'key_added': key_added,
            'obsm_keys': obsm_keys,
            'assigned_points': assigned,
            'total_points': total_points,
            'core_cells': n_core,
            'clustering_auto_run': needs_clustering,
            'cluster_key': cluster_key,
            'message': f'Core analysis complete. Embedding stored in obsm["{key_added}"] ({assigned}/{total_points} points assigned, {n_core} core cells)'
        }
        
        # Free memory before returning
        del adata
        if model is not None:
            del model
        if 'embedding' in locals():
            del embedding
        gc.collect()
        sys.stdout.flush()
        
        return jsonify(response)
    
    except MemoryError as e:
        error_msg = 'Out of memory: The dataset is too large for the current server configuration. Please try with a smaller dataset or contact support to upgrade server resources.'
        print(f"[CORE] Memory Error: {error_msg}")
        print(traceback.format_exc())
        sys.stdout.flush()
        return jsonify({
            'error': error_msg,
            'error_type': 'MemoryError',
            'suggestion': 'Try with a smaller dataset or use data subsampling'
        }), 500
    except TimeoutError as e:
        error_msg = 'Operation timeout: The computation took too long. Please try with a smaller dataset or lower resolution parameters.'
        print(f"[CORE] Timeout Error: {error_msg}")
        print(traceback.format_exc())
        sys.stdout.flush()
        return jsonify({
            'error': error_msg,
            'error_type': 'TimeoutError',
            'suggestion': 'Try with a smaller dataset or reduce cluster_resolution'
        }), 500
    except Exception as e:
        error_msg = f'Core analysis failed: {str(e)}'
        print(f"[CORE] Error: {error_msg}")
        print(traceback.format_exc())
        sys.stdout.flush()
        
        # Provide more helpful error messages
        error_type = type(e).__name__
        suggestion = None
        
        if 'neighbors' in str(e).lower():
            suggestion = 'Please run preprocessing first to compute neighbors graph.'
        elif 'representation' in str(e).lower() or 'X_pca' in str(e) or 'X_latent' in str(e):
            suggestion = 'Please run preprocessing first to generate PCA or latent representation.'
        elif 'clustering' in str(e).lower():
            suggestion = 'Please run clustering first, or ensure clustering parameters are correct.'
        
        return jsonify({
            'error': error_msg,
            'error_type': error_type,
            'suggestion': suggestion
        }), 500


def _validate_session(session_id: str) -> tuple[bool, str]:
    """Validate that a session exists and has data
    
    Returns:
        (is_valid, error_message)
    """
    try:
        adata = load_adata(session_id)
        if adata is None:
            return False, f'Session {session_id} not found or has no data'
        return True, None
    except Exception as e:
        return False, f'Error validating session: {str(e)}'


@bp.route('/cluster-keys', methods=['POST'])
def list_cluster_keys():
    """List all cluster keys for a session
    
    Request:
        POST: {"session_id": "xxx"}
    
    Returns:
        {
            "success": true,
            "session_id": "string",
            "cluster_keys": ["key1", "key2", ...],
            "n_keys": int
        }
    """
    try:
        data = request.json
        session_id = data.get('session_id', 'default')
        
        from ..utils import load_adata
        adata = load_adata(session_id)
        
        if adata is None:
            return jsonify({'error': 'No data found.'}), 400
        
        # Find all cluster-related keys in adata.obs
        # Common cluster keys: cplearn, leiden, louvain, kmeans, etc.
        cluster_keys_list = []
        cluster_keys_info = []
        
        # Exclude patterns that are NOT cluster keys
        exclude_patterns = [
            '_is_core',  # Boolean markers for core cells
            '_model_info',  # Model metadata
            'n_genes',  # QC metrics
            'n_counts',  # QC metrics
            'total_counts',  # QC metrics
            'pct_counts',  # QC metrics
            'log1p_',  # Normalized QC metrics
        ]
        
        # Known non-cluster keys (metadata, QC, etc.)
        exclude_keys = [
            'truth',  # Ground truth labels (not clustering results)
            'ground_truth',  # Ground truth labels
            'cell_type',  # Cell type annotations (not clustering)
            'celltype',  # Cell type annotations
        ]
        
        for key in adata.obs.columns:
            # Skip excluded keys
            if key in exclude_keys:
                continue
            
            # Skip keys matching exclude patterns
            if any(pattern in key for pattern in exclude_patterns):
                continue
            
            # Check if it looks like a cluster key
            # (categorical or has integer/string values that could be cluster labels)
            is_cluster_key = False
            method = 'unknown'
            
            # Known cluster key names
            if key in ['cplearn', 'leiden', 'louvain', 'kmeans']:
                is_cluster_key = True
                method = key
            elif adata.obs[key].dtype.name == 'category':
                # Categorical columns might be cluster keys
                unique_values = adata.obs[key].unique()
                n_unique = len(unique_values)
                
                # Additional checks to avoid false positives:
                # 1. Must have reasonable number of unique values (3-100 for clusters, not just 2)
                # 2. Should not be boolean-like (only True/False)
                # 3. Should have integer-like or string labels
                if 3 <= n_unique <= 100:
                    # Check if values look like cluster labels (integers or strings, not just booleans)
                    sample_values = list(unique_values[:10])  # Sample first 10 values
                    has_non_bool = False
                    for val in sample_values:
                        # Check if value is not just boolean
                        if val not in [True, False, 'True', 'False', 'true', 'false']:
                            has_non_bool = True
                            break
                    
                    if has_non_bool:
                        is_cluster_key = True
                        # Try to detect method from key name
                        key_lower = key.lower()
                        if 'cplearn' in key_lower:
                            method = 'cplearn'
                        elif 'leiden' in key_lower:
                            method = 'leiden'
                        elif 'louvain' in key_lower:
                            method = 'louvain'
                        elif 'kmeans' in key_lower:
                            method = 'kmeans'
                        else:
                            method = 'unknown'
            
            if is_cluster_key:
                cluster_keys_list.append(key)
                n_clusters = len(adata.obs[key].unique())
                cluster_keys_info.append({
                    'key': key,
                    'method': method,
                    'n_clusters': n_clusters
                })
        
        # Remove duplicates and sort
        cluster_keys_list = sorted(list(set(cluster_keys_list)))
        cluster_keys_info = sorted(cluster_keys_info, key=lambda x: x['key'])
        
        return jsonify({
            'success': True,
            'session_id': session_id,
            'cluster_keys': cluster_keys_info,  # Return detailed info for frontend
            'n_keys': len(cluster_keys_info)
        })
    
    except Exception as e:
        import traceback
        error_msg = f'Failed to list cluster keys: {str(e)}'
        print(f"[CORE] Error: {error_msg}")
        print(traceback.format_exc())
        return jsonify({
            'error': error_msg,
            'details': str(e)
        }), 500
