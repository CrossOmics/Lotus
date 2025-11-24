"""
Analysis API endpoints (marker-genes)
"""

import numpy as np
from flask import Blueprint, request, jsonify
from ..config import LOTUS_AVAILABLE, clustering
from ..utils import load_adata, save_adata

bp = Blueprint('analysis', __name__, url_prefix='/api')


@bp.route('/available-cluster-keys', methods=['POST'])
def get_available_cluster_keys():
    """Get all available cluster keys from adata.obs"""
    try:
        data = request.json
        session_id = data.get('session_id', 'default')
        
        from ..utils import load_adata
        adata = load_adata(session_id)
        
        if adata is None:
            return jsonify({'error': 'No data found.'}), 400
        
        # Find all cluster-related keys in adata.obs
        # Common cluster keys: cplearn, leiden, louvain, kmeans, etc.
        cluster_keys = []
        for key in adata.obs.columns:
            # Check if it looks like a cluster key
            # (categorical or has integer/string values that could be cluster labels)
            if key in ['cplearn', 'leiden', 'louvain', 'kmeans']:
                unique_vals = adata.obs[key].unique()
                if len(unique_vals) > 1:  # At least 2 clusters
                    cluster_keys.append({
                        'key': key,
                        'n_clusters': len(unique_vals),
                        'clusters': [str(v) for v in sorted(unique_vals)],
                        'method': _infer_clustering_method(key)
                    })
            # Also check for other potential cluster keys
            elif key.endswith('_labels') or key.endswith('_clusters'):
                unique_vals = adata.obs[key].unique()
                if len(unique_vals) > 1:
                    cluster_keys.append({
                        'key': key,
                        'n_clusters': len(unique_vals),
                        'clusters': [str(v) for v in sorted(unique_vals)],
                        'method': _infer_clustering_method(key)
                    })
        
        return jsonify({
            'success': True,
            'cluster_keys': cluster_keys,
            'message': f'Found {len(cluster_keys)} cluster key(s)'
        })
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


def _infer_clustering_method(key: str) -> str:
    """Infer clustering method from key name"""
    if 'cplearn' in key.lower():
        return 'cplearn'
    elif 'leiden' in key.lower():
        return 'leiden'
    elif 'louvain' in key.lower():
        return 'louvain'
    elif 'kmeans' in key.lower():
        return 'kmeans'
    else:
        return 'unknown'


@bp.route('/marker-genes', methods=['POST'])
def run_marker_genes():
    """Run marker gene identification (DEG analysis)"""
    try:
        if not LOTUS_AVAILABLE:
            return jsonify({'error': 'Lotus not available'}), 500
        
        from lotus.workflows import marker_genes
        
        data = request.json
        session_id = data.get('session_id', 'default')
        cluster_key = data.get('cluster_key', 'cplearn')
        
        adata = load_adata(session_id)
        
        if adata is None:
            return jsonify({'error': 'No data found.'}), 400
        
        # Check if cluster_key exists in adata.obs
        if cluster_key not in adata.obs:
            available_keys = [k for k in adata.obs.columns if k in ['cplearn', 'leiden', 'louvain', 'kmeans'] or k.endswith('_labels') or k.endswith('_clusters')]
            error_msg = f'Cluster key "{cluster_key}" not found in adata.obs. Please run clustering first.'
            if available_keys:
                error_msg += f' Available cluster keys: {", ".join(available_keys)}'
            print(f"[DEG] ERROR: {error_msg}")
            return jsonify({
                'error': error_msg,
                'available_keys': available_keys
            }), 400
        
        # Debug: Check cluster labels and data
        labels = adata.obs[cluster_key]
        unique_labels = sorted([l for l in labels.unique() if l != -1])
        label_counts = labels.value_counts().to_dict()
        print(f"[DEG] Cluster key '{cluster_key}' found with {len(unique_labels)} unique labels: {unique_labels}")
        print(f"[DEG] Label counts: {label_counts}")
        print(f"[DEG] Data shape: {adata.shape}, n_obs={adata.n_obs}, n_vars={adata.n_vars}")
        print(f"[DEG] Available layers: {list(adata.layers.keys())}")
        
        if len(unique_labels) < 2:
            error_msg = f'Need at least 2 clusters for DEG analysis. Found {len(unique_labels)} cluster(s): {unique_labels}'
            print(f"[DEG] ERROR: {error_msg}")
            return jsonify({
                'error': error_msg,
                'unique_labels': unique_labels,
                'label_counts': label_counts
            }), 400
        
        # Check if cplearn layers info exists in adata.uns
        # If not, create a default layers structure (all cells in one layer)
        # Note: load_adata() should have already restored layers_json if present
        layers_key = f"{cluster_key}_cplearn"
        if layers_key not in adata.uns:
            print(f"[DEG] Creating default layers structure for {layers_key}")
            # Create a default layers structure where all cells are in one layer
            adata.uns[layers_key] = {
                'layers': [[i for i in range(adata.n_obs)]]
            }
            # Save the updated adata
            save_adata(adata, session_id)
        else:
            # Check if layers exist, if not create default
            layers_meta = adata.uns[layers_key]
            if not isinstance(layers_meta, dict):
                print(f"[DEG] Converting {layers_key} to dict format")
                adata.uns[layers_key] = {'layers': [[i for i in range(adata.n_obs)]]}
                save_adata(adata, session_id)
            elif 'layers' not in layers_meta or layers_meta.get('layers') is None:
                print(f"[DEG] Adding default layers to existing {layers_key}")
                layers_meta['layers'] = [[i for i in range(adata.n_obs)]]
                save_adata(adata, session_id)
            elif isinstance(layers_meta.get('layers'), str):
                # If layers is still a JSON string (shouldn't happen after load_adata fix, but just in case)
                print(f"[DEG] Restoring layers from JSON string in {layers_key}")
                import json as json_lib
                try:
                    layers_meta['layers'] = json_lib.loads(layers_meta['layers'])
                    save_adata(adata, session_id)
                except Exception as e:
                    print(f"[DEG] Warning: Could not restore layers from JSON: {e}")
                    layers_meta['layers'] = [[i for i in range(adata.n_obs)]]
                    save_adata(adata, session_id)
        
        # Run marker genes analysis
        print(f"[DEG] Running marker genes analysis with cluster_key={cluster_key}")
        # Auto-detect layer (prefer raw_counts for DEG analysis, matching lotus_workflow.py)
        layer = None
        for layer_name in ["raw_counts", "raw", "counts"]:
            if layer_name in adata.layers:
                layer = layer_name
                print(f"[DEG] Using layer: {layer_name}")
                break
        
        # Use exact same parameters as lotus_workflow.py
        print(f"[DEG] Calling marker_genes with: cluster_key={cluster_key}, layer={layer}, min_detect_pct=0.0, min_cells_per_group=5")
        
        # Additional validation before calling marker_genes
        labels_check = adata.obs[cluster_key]
        unique_check = sorted([l for l in labels_check.unique() if l != -1])
        print(f"[DEG] Pre-call validation: {len(unique_check)} unique labels, labels={unique_check}")
        
        # Check layer data
        if layer:
            layer_data = adata.layers[layer]
            print(f"[DEG] Layer '{layer}' shape: {layer_data.shape}, dtype: {layer_data.dtype}")
            # Handle sparse matrices properly
            try:
                from scipy.sparse import issparse
                if issparse(layer_data):
                    non_zero_count = layer_data.nnz
                    total_size = layer_data.shape[0] * layer_data.shape[1]
                else:
                    non_zero_count = np.count_nonzero(layer_data)
                    total_size = layer_data.size
                print(f"[DEG] Layer '{layer}' non-zero: {non_zero_count}/{total_size}")
            except Exception as e:
                print(f"[DEG] Could not count non-zero values: {e}")
        
        try:
            result = marker_genes(
                adata,
                cluster_key=cluster_key,
                layer=layer,  # Explicitly pass layer (should be "raw_counts" if available)
                min_detect_pct=0.0,  # Match lotus_workflow.py exactly
                min_cells_per_group=5,  # Match lotus_workflow.py exactly
                auto_pick_groups=True  # Match lotus_workflow.py exactly
            )
            print(f"[DEG] marker_genes returned: type={type(result)}, length={len(result) if hasattr(result, '__len__') else 'N/A'}")
            if hasattr(result, '__len__') and len(result) == 0:
                print(f"[DEG] WARNING: Empty result! Checking why...")
                # Re-check labels after marker_genes call
                labels_after = adata.obs[cluster_key]
                unique_after = sorted([l for l in labels_after.unique() if l != -1])
                print(f"[DEG] Labels after call: {len(unique_after)} unique labels, labels={unique_after}")
        except Exception as e:
            error_msg = f'Marker genes analysis failed: {str(e)}'
            print(f"[DEG] ERROR in marker_genes call: {error_msg}")
            import traceback
            traceback.print_exc()
            return jsonify({
                'error': error_msg,
                'details': str(e)
            }), 500
        
        # Save updated AnnData
        save_adata(adata, session_id)
        
        # Format results for JSON response (matching lotus_workflow.py output format)
        markers_dict = {}
        stats = {
            'total': 0,
            'significant_p05': 0,
            'significant_p01': 0
        }
        
        if result is not None and hasattr(result, 'to_dict'):
            print(f"[DEG] Result type: {type(result)}, length: {len(result) if hasattr(result, '__len__') else 'N/A'}")
            if len(result) > 0:
                print(f"[DEG] Result columns: {list(result.columns)}")
                
                # Calculate statistics matching lotus_workflow.py
                stats['total'] = len(result)
                if 'p_adj' in result.columns:
                    stats['significant_p05'] = int(np.sum(result['p_adj'] < 0.05))
                    stats['significant_p01'] = int(np.sum(result['p_adj'] < 0.01))
                
                print(f"[DEG] Statistics: total={stats['total']}, p_adj<0.05={stats['significant_p05']}, p_adj<0.01={stats['significant_p01']}")
                
                # Get top genes (matching lotus_workflow.py which shows top 10)
                # Note: result columns are: gene, log2fc, z_score, pvalue, p_adj, mean_a, mean_b, pct_expr_a, pct_expr_b, ...
                for idx, row in result.head(50).iterrows():
                    gene_name = row.get('gene', idx) if 'gene' in row else str(idx)
                    log2fc_val = float(row.get('log2fc', 0)) if 'log2fc' in row else 0.0
                    pval_val = float(row.get('pvalue', 1)) if 'pvalue' in row else 1.0
                    pval_adj_val = float(row.get('p_adj', 1)) if 'p_adj' in row else 1.0
                    
                    markers_dict[str(idx)] = {
                        'gene': gene_name,
                        'log2fc': log2fc_val,
                        'pval': pval_val,
                        'pval_adj': pval_adj_val,
                    }
                print(f"[DEG] Formatted {len(markers_dict)} markers for response")
            else:
                print("[DEG] WARNING: Result DataFrame is empty!")
        else:
            print(f"[DEG] WARNING: Result is None or not a DataFrame! Type: {type(result)}")
        
        response_data = {
            'success': True,
            'markers': markers_dict,
            'n_markers': len(markers_dict),
            'stats': stats,  # Add statistics matching lotus_workflow.py
            'message': f'Found {len(markers_dict)} marker genes (Total: {stats["total"]}, p_adj<0.05: {stats["significant_p05"]}, p_adj<0.01: {stats["significant_p01"]})'
        }
        print(f"[DEG] Response: n_markers={response_data['n_markers']}, stats={stats}, markers keys={list(markers_dict.keys())[:5] if markers_dict else 'empty'}")
        return jsonify(response_data)
    
    except Exception as e:
        import traceback
        error_msg = f'Marker genes analysis failed: {str(e)}'
        print(f"[DEG] ERROR: {error_msg}")
        traceback.print_exc()
        return jsonify({
            'error': error_msg,
            'details': str(e)
        }), 500

