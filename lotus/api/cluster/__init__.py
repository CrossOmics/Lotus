"""
Cluster API endpoints
"""

import gc
import sys
from flask import Blueprint, request, jsonify
import numpy as np
import pandas as pd
from ..config import LOTUS_AVAILABLE, clustering, sc
from ..utils import load_adata, save_adata

bp = Blueprint('cluster', __name__, url_prefix='/api')


@bp.route('/cluster', methods=['POST'])
def run_clustering():
    """Run clustering analysis"""
    if not LOTUS_AVAILABLE:
        return jsonify({'error': 'Lotus not available'}), 500
    
    try:
        data = request.json
        session_id = data.get('session_id', 'default')
        method = data.get('method', 'cplearn')  # 'cplearn', 'leiden', 'louvain', 'kmeans'
        resolution = data.get('resolution', 1.2)
        use_rep = data.get('use_rep', None)
        key_added = data.get('key_added', None)  # If None, uses method-specific default
        # cplearn-specific parameters
        stable_core_frac = data.get('stable_core_frac', 0.25)
        stable_ng_num = data.get('stable_ng_num', 4)
        fine_grained = data.get('fine_grained', False)
        propagate = data.get('propagate', True)
        # scanpy-specific parameters
        random_state = data.get('random_state', 0)
        
        adata = load_adata(session_id)
        
        if adata is None:
            return jsonify({'error': 'No data found. Please upload and preprocess data first.'}), 400
        
        # Run clustering based on method
        model = None
        if method == 'cplearn':
            cluster_key_to_use = key_added or 'cplearn'
            
            # Check if cplearn clustering already exists (from core selection)
            # Check both the cluster key and model_info
            has_existing_clustering = False
            if cluster_key_to_use in adata.obs:
                model_info_key = f'{cluster_key_to_use}_model_info'
                if model_info_key in adata.uns:
                    model_info = adata.uns[model_info_key]
                    # Check if parameters match (allow some flexibility)
                    existing_resolution = model_info.get('resolution', None)
                    existing_use_rep = model_info.get('use_rep', None)
                    
                    # Auto-detect use_rep for comparison
                    rep_to_check = use_rep or ('X_latent' if 'X_latent' in adata.obsm else 'X_pca' if 'X_pca' in adata.obsm else None)
                    
                    # If resolution matches (or close) and use_rep matches, reuse existing clustering
                    if existing_resolution is not None:
                        resolution_match = abs(existing_resolution - resolution) < 0.1  # Allow 0.1 tolerance
                        rep_match = (existing_use_rep == rep_to_check) or (existing_use_rep is None and rep_to_check is None)
                        
                        if resolution_match and rep_match:
                            print(f"[CLUSTER] Found existing cplearn clustering from core selection (resolution={existing_resolution}, use_rep={existing_use_rep})")
                            print(f"[CLUSTER] Reusing existing clustering instead of re-running")
                            has_existing_clustering = True
                            cluster_key = cluster_key_to_use
                        else:
                            print(f"[CLUSTER] Existing clustering found but parameters differ:")
                            print(f"[CLUSTER]   Existing: resolution={existing_resolution}, use_rep={existing_use_rep}")
                            print(f"[CLUSTER]   Requested: resolution={resolution}, use_rep={rep_to_check}")
                            print(f"[CLUSTER]   Will re-run clustering with new parameters")
                    else:
                        # If model_info exists but no resolution info, assume compatible
                        print(f"[CLUSTER] Found existing cplearn clustering from core selection")
                        print(f"[CLUSTER] Reusing existing clustering")
                        has_existing_clustering = True
                        cluster_key = cluster_key_to_use
            
            # If existing clustering found and compatible, skip re-running
            if has_existing_clustering:
                # Get cluster information from existing clustering
                clusters = adata.obs[cluster_key].unique()
                cluster_counts = adata.obs[cluster_key].value_counts().to_dict()
                n_clusters = len(clusters)
                clusters_list = [str(c) for c in clusters]
                cluster_counts_dict = {str(k): int(v) for k, v in cluster_counts.items()}
                
                # Check if model exists (for core selection compatibility)
                model_info_key = f'{cluster_key}_model_info'
                has_model = model_info_key in adata.uns
                
                # Free memory before returning
                del adata
                gc.collect()
                sys.stdout.flush()
                
                return jsonify({
                    'success': True,
                    'cluster_key': cluster_key,
                    'n_clusters': n_clusters,
                    'clusters': clusters_list,
                    'cluster_counts': cluster_counts_dict,
                    'has_model': has_model,
                    'reused': True,  # Indicate that existing clustering was reused
                    'message': f'Using existing cplearn clustering: {n_clusters} clusters found (from core selection)'
                })
            
            # No existing clustering found, proceed with normal clustering
            # Ensure neighbors graph exists before cplearn clustering
            if 'neighbors' not in adata.uns:
                print("[CLUSTER] Neighbors graph not found, computing it first...")
                if LOTUS_AVAILABLE:
                    from lotus.workflows.preprocess import neighbors as compute_neighbors
                    rep_to_use = use_rep or ('X_pca' if 'X_pca' in adata.obsm else 'X_latent' if 'X_latent' in adata.obsm else None)
                    if rep_to_use:
                        compute_neighbors(adata, use_rep=rep_to_use, n_neighbors=15)
                    else:
                        sc.pp.neighbors(adata, n_neighbors=15)
                else:
                    sc.pp.neighbors(adata, n_neighbors=15)
            
            print(f"[CLUSTER] Running cplearn clustering with resolution={resolution}, use_rep={use_rep}")
            try:
                # Ensure data is in correct format for cplearn
                rep_to_use = use_rep or ('X_latent' if 'X_latent' in adata.obsm else 'X_pca' if 'X_pca' in adata.obsm else None)
                if rep_to_use and rep_to_use in adata.obsm:
                    # Check for NaN values
                    data_matrix = adata.obsm[rep_to_use]
                    if np.isnan(data_matrix).any():
                        print(f"[CLUSTER] Warning: NaN values found in {rep_to_use}, filling with 0")
                        data_matrix = np.nan_to_num(data_matrix, nan=0.0)
                        adata.obsm[rep_to_use] = data_matrix
                
                # Use parameters from request (cluster_key_to_use already defined above)
                model = clustering(
                    adata,
                    method='cplearn',
                    use_rep=rep_to_use,
                    key_added=cluster_key_to_use,
                    cluster_resolution=resolution,
                    stable_core_frac=stable_core_frac,
                    stable_ng_num=stable_ng_num,
                    fine_grained=fine_grained,
                    propagate=propagate,
                    print_summary=False  # Don't print in API
                )
                cluster_key = cluster_key_to_use
            except Exception as e:
                import traceback
                error_msg = f'cplearn clustering failed: {str(e)}'
                print(f"[CLUSTER] Error: {error_msg}")
                print(traceback.format_exc())
                # Fallback to scanpy if cplearn fails
                print("[CLUSTER] Falling back to scanpy Leiden...")
                if sc is None:
                    return jsonify({'error': error_msg, 'suggestion': 'Try using scanpy method instead'}), 500
                sc.tl.leiden(adata, resolution=resolution, key_added='leiden')
                cluster_key = 'leiden'
        
        elif method == 'leiden':
            if sc is None:
                return jsonify({'error': 'scanpy not available'}), 500
            try:
                # Ensure neighbors graph exists
                if 'neighbors' not in adata.uns:
                    print("[CLUSTER] Computing neighbors graph for Leiden...")
                    sc.pp.neighbors(adata, use_rep=use_rep or 'X_pca', n_neighbors=15)
                cluster_key_to_use = key_added or 'leiden'
                print(f"[CLUSTER] Running scanpy Leiden with resolution={resolution}, key_added={cluster_key_to_use}")
                sc.tl.leiden(adata, resolution=resolution, key_added=cluster_key_to_use, random_state=random_state)
                cluster_key = cluster_key_to_use
            except Exception as e:
                error_msg = f'Leiden clustering failed: {str(e)}'
                print(f"[CLUSTER] ERROR: {error_msg}")
                import traceback
                traceback.print_exc()
                return jsonify({
                    'error': error_msg,
                    'details': str(e)
                }), 500
        
        elif method == 'louvain':
            if sc is None:
                return jsonify({'error': 'scanpy not available'}), 500
            try:
                # Check if louvain module is available (pre-check)
                try:
                    import louvain
                    print("[CLUSTER] louvain module is available")
                except (ImportError, ModuleNotFoundError):
                    print("[CLUSTER] WARNING: louvain module not found, but will try anyway (scanpy may handle it)")
                
                # Ensure neighbors graph exists
                if 'neighbors' not in adata.uns:
                    print("[CLUSTER] Computing neighbors graph for Louvain...")
                    sc.pp.neighbors(adata, use_rep=use_rep or 'X_pca', n_neighbors=15)
                cluster_key_to_use = key_added or 'louvain'
                print(f"[CLUSTER] Running scanpy Louvain with resolution={resolution}, key_added={cluster_key_to_use}")
                sc.tl.louvain(adata, resolution=resolution, key_added=cluster_key_to_use, random_state=random_state)
                cluster_key = cluster_key_to_use
            except (ImportError, ModuleNotFoundError) as e:
                error_msg = f'Louvain clustering failed: louvain module not installed'
                print(f"[CLUSTER] ERROR: {error_msg}")
                print(f"[CLUSTER] Error details: {str(e)}")
                return jsonify({
                    'error': error_msg,
                    'details': str(e),
                    'suggestion': 'Please install louvain: pip install louvain (Note: this is different from python-louvain)'
                }), 500
            except Exception as e:
                error_msg = f'Louvain clustering failed: {str(e)}'
                print(f"[CLUSTER] ERROR: {error_msg}")
                import traceback
                traceback.print_exc()
                # Check if it's actually a module import error
                if 'louvain' in str(e).lower() or 'No module named' in str(e):
                    return jsonify({
                        'error': 'louvain module not installed',
                        'details': str(e),
                        'suggestion': 'Please install louvain: pip install louvain (Note: this is different from python-louvain)'
                    }), 500
                return jsonify({
                    'error': error_msg,
                    'details': str(e)
                }), 500
        
        elif method == 'kmeans':
            if sc is None:
                return jsonify({'error': 'scanpy not available'}), 500
            try:
                from sklearn.cluster import KMeans
                print(f"[CLUSTER] Running KMeans clustering")
                # Determine number of clusters from resolution (rough estimate)
                n_clusters = max(2, int(resolution * 5))
                # Use representation
                rep_to_use = use_rep or ('X_pca' if 'X_pca' in adata.obsm else 'X_latent' if 'X_latent' in adata.obsm else None)
                if rep_to_use and rep_to_use in adata.obsm:
                    X = adata.obsm[rep_to_use]
                else:
                    X = adata.X
                cluster_key_to_use = key_added or 'kmeans'
                # Run KMeans
                kmeans = KMeans(n_clusters=n_clusters, random_state=random_state, n_init=10)
                labels = kmeans.fit_predict(X)
                adata.obs[cluster_key_to_use] = pd.Categorical(labels)
                cluster_key = cluster_key_to_use
            except ImportError as e:
                error_msg = f'KMeans clustering failed: {str(e)}'
                print(f"[CLUSTER] ERROR: {error_msg}")
                return jsonify({
                    'error': error_msg,
                    'suggestion': 'Please install scikit-learn: pip install scikit-learn'
                }), 500
            except Exception as e:
                error_msg = f'KMeans clustering failed: {str(e)}'
                print(f"[CLUSTER] ERROR: {error_msg}")
                import traceback
                traceback.print_exc()
                return jsonify({
                    'error': error_msg,
                    'details': str(e)
                }), 500
        
        else:
            return jsonify({'error': f'Unknown clustering method: {method}'}), 400
        
        # Save updated AnnData
        save_adata(adata, session_id)
        gc.collect()  # Free memory after first save
        
        # Store model for core selection (if cplearn was used)
        if model is not None:
            # Store model metadata in uns for later use
            adata.uns[f'{cluster_key}_model_info'] = {
                'method': 'cplearn',
                'resolution': resolution,
                'use_rep': use_rep
            }
            save_adata(adata, session_id)
            gc.collect()  # Free memory after second save
        
        # Get cluster information before freeing memory
        clusters = adata.obs[cluster_key].unique()
        cluster_counts = adata.obs[cluster_key].value_counts().to_dict()
        n_clusters = len(clusters)
        clusters_list = [str(c) for c in clusters]
        cluster_counts_dict = {str(k): int(v) for k, v in cluster_counts.items()}
        has_model = model is not None
        
        # Free memory
        del adata
        if model is not None:
            del model
        gc.collect()
        sys.stdout.flush()
        
        return jsonify({
            'success': True,
            'cluster_key': cluster_key,
            'n_clusters': n_clusters,
            'clusters': clusters_list,
            'cluster_counts': cluster_counts_dict,
            'has_model': has_model,  # Indicates if core selection is available
            'message': f'Clustering complete: {n_clusters} clusters found using {method}'
        })
    
    except Exception as e:
        import traceback
        error_msg = f'Clustering failed: {str(e)}'
        print(f"[CLUSTER] ERROR: {error_msg}")
        traceback.print_exc()
        return jsonify({
            'error': error_msg,
            'details': str(e)
        }), 500

