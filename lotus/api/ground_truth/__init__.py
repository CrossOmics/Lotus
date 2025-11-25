"""
Ground Truth Management API endpoints
Handles uploading, querying, and managing ground truth labels for user sessions
"""

from typing import Optional, List, Dict, Any
import json
import pandas as pd
import numpy as np
from flask import Blueprint, request, jsonify
from ..utils import load_adata, save_adata, get_session_dir
from ..session import _session_access_times, _session_lock
import time

bp = Blueprint('ground_truth', __name__, url_prefix='/api')


def _validate_session(session_id: str):
    """Validate that session exists and has data
    
    Returns:
        (is_valid, error_message)
    """
    if not session_id or session_id == '':
        return False, 'session_id is required'
    
    # Check if session directory exists
    session_dir = get_session_dir(session_id)
    if not session_dir.exists():
        return False, f'Session {session_id} does not exist. Please upload data first.'
    
    # Check if data exists
    adata = load_adata(session_id)
    if adata is None:
        return False, f'No data found for session {session_id}. Please upload data first.'
    
    # Update session access time
    with _session_lock:
        _session_access_times[session_id] = time.time()
    
    return True, None


def _parse_truth_json(truth_json_content: str):
    """Parse ground truth JSON content
    
    Supports multiple formats:
    1. Direct list: ["label1", "label2", ...]
    2. Dict with labels: {"labels": [...]}
    3. Dict with key and labels: {"key": "cell_type", "labels": [...]}
    4. Dict with ground_truth: {"ground_truth": [...]}
    5. Dict mapping cell IDs to labels: {"cell1": "label1", "cell2": "label2", ...}
    
    Returns:
        (truth_labels, truth_key, error_message)
    """
    try:
        truth_data = json.loads(truth_json_content)
        
        truth_labels = None
        truth_key = None
        
        if isinstance(truth_data, list):
            # Format 1: Direct list of labels
            truth_labels = truth_data
        elif isinstance(truth_data, dict):
            # Extract key from JSON if present
            if 'key' in truth_data:
                truth_key = truth_data['key']
            
            # Extract labels
            if 'labels' in truth_data:
                truth_labels = truth_data['labels']
            elif 'ground_truth' in truth_data:
                truth_labels = truth_data['ground_truth']
            else:
                # Format 5: Assume values are labels, keys are cell IDs
                # We'll use values as labels (assuming they're in the same order as cells)
                truth_labels = list(truth_data.values())
        
        if truth_labels is None:
            return None, None, 'Invalid JSON format: Could not extract labels. Expected list or dict with "labels" or "ground_truth" field.'
        
        if not isinstance(truth_labels, list):
            return None, None, f'Invalid JSON format: Labels must be a list, got {type(truth_labels).__name__}'
        
        if len(truth_labels) == 0:
            return None, None, 'Invalid JSON format: Labels list is empty'
        
        return truth_labels, truth_key, None
    
    except json.JSONDecodeError as e:
        return None, None, f'Invalid JSON format: {str(e)}'
    except Exception as e:
        return None, None, f'Error parsing JSON: {str(e)}'


@bp.route('/ground-truth/upload', methods=['POST'])
def upload_ground_truth():
    """Upload ground truth labels for a session
    
    Supports multiple input formats:
    1. JSON file upload (multipart/form-data)
    2. JSON string in request body (application/json)
    3. Direct labels array in request body
    
    Request formats:
    - Form data:
        - session_id: str (required)
        - truth_json_file: file (optional, JSON file)
        - truth_key: str (optional, default: 'ground_truth')
    
    - JSON body:
        {
            "session_id": "string" (required),
            "labels": ["label1", "label2", ...] (optional),
            "truth_json": "string" (optional, JSON string),
            "truth_key": "string" (optional, default: "ground_truth")
        }
    
    Returns:
        {
            "success": true,
            "session_id": "string",
            "truth_key": "string",
            "n_labels": int,
            "message": "string"
        }
    """
    try:
        # Handle both JSON and form data
        if request.is_json:
            data = request.json
            session_id = data.get('session_id')
            truth_key = data.get('truth_key', 'ground_truth')
            truth_json_content = data.get('truth_json', None)
            labels = data.get('labels', None)
            
            # If labels are provided directly, use them
            if labels is not None:
                truth_labels = labels
            elif truth_json_content:
                truth_labels, json_truth_key, error = _parse_truth_json(truth_json_content)
                if error:
                    return jsonify({'error': error}), 400
                # Use truth_key from JSON if provided, otherwise use request parameter
                if json_truth_key:
                    truth_key = json_truth_key
            else:
                return jsonify({'error': 'Either "labels", "truth_json", or "truth_json_file" must be provided'}), 400
        else:
            # Form data (for file upload)
            data = request.form
            session_id = data.get('session_id')
            truth_key = data.get('truth_key', 'ground_truth')
            truth_labels = None
            
            # Handle JSON file upload
            if 'truth_json_file' in request.files:
                json_file = request.files['truth_json_file']
                if json_file.filename:
                    try:
                        truth_json_content = json_file.read().decode('utf-8')
                        truth_labels, json_truth_key, error = _parse_truth_json(truth_json_content)
                        if error:
                            return jsonify({'error': error}), 400
                        # Use truth_key from JSON if provided, otherwise use form parameter
                        if json_truth_key:
                            truth_key = json_truth_key
                    except Exception as e:
                        return jsonify({'error': f'Failed to read JSON file: {str(e)}'}), 400
                else:
                    return jsonify({'error': 'No file provided in truth_json_file'}), 400
            else:
                return jsonify({'error': 'Either "truth_json_file" (form data) or "labels"/"truth_json" (JSON body) must be provided'}), 400
        
        # Validate session
        is_valid, error_msg = _validate_session(session_id)
        if not is_valid:
            return jsonify({'error': error_msg}), 400
        
        # Load adata
        adata = load_adata(session_id)
        
        # Validate labels length
        if len(truth_labels) != adata.n_obs:
            return jsonify({
                'error': f'Labels length ({len(truth_labels)}) does not match number of cells ({adata.n_obs})',
                'n_cells': adata.n_obs,
                'n_labels': len(truth_labels)
            }), 400
        
        # Validate truth_key
        if not truth_key or truth_key.strip() == '':
            truth_key = 'ground_truth'
        
        # Save truth labels to adata.obs
        adata.obs[truth_key] = pd.Categorical(truth_labels)
        save_adata(adata, session_id)
        
        print(f"[GROUND_TRUTH] Saved ground truth labels to adata.obs['{truth_key}'] for session {session_id}")
        
        return jsonify({
            'success': True,
            'session_id': session_id,
            'truth_key': truth_key,
            'n_labels': len(truth_labels),
            'n_cells': adata.n_obs,
            'message': f'Ground truth labels saved successfully to adata.obs["{truth_key}"]'
        })
    
    except Exception as e:
        import traceback
        error_msg = f'Failed to upload ground truth: {str(e)}'
        print(f"[GROUND_TRUTH] Error: {error_msg}")
        print(traceback.format_exc())
        return jsonify({'error': error_msg}), 500


@bp.route('/ground-truth/list', methods=['GET', 'POST'])
def list_ground_truth():
    """List all ground truth keys for a session
    
    Request:
        GET: ?session_id=xxx
        POST: {"session_id": "xxx"}
    
    Returns:
        {
            "success": true,
            "session_id": "string",
            "truth_keys": ["key1", "key2", ...],
            "n_keys": int
        }
    """
    try:
        # Get session_id from query params or JSON body
        if request.method == 'GET':
            session_id = request.args.get('session_id')
        else:
            data = request.json or {}
            session_id = data.get('session_id')
        
        if not session_id:
            return jsonify({'error': 'session_id is required'}), 400
        
        # Validate session
        is_valid, error_msg = _validate_session(session_id)
        if not is_valid:
            return jsonify({'error': error_msg}), 400
        
        # Load adata
        adata = load_adata(session_id)
        
        # Find all ground truth keys (keys that look like ground truth)
        # Common patterns: ground_truth, truth, cell_type, cell_type_truth, etc.
        truth_keys = []
        for key in adata.obs.columns:
            # Check if key looks like a ground truth label
            key_lower = key.lower()
            if any(pattern in key_lower for pattern in ['truth', 'ground', 'label', 'type', 'annotation']):
                # Check if it's categorical (typical for labels)
                if isinstance(adata.obs[key].dtype, pd.CategoricalDtype):
                    truth_keys.append(key)
        
        # If no obvious ground truth keys found, return empty list
        # (User might have custom key names)
        
        return jsonify({
            'success': True,
            'session_id': session_id,
            'truth_keys': truth_keys,
            'n_keys': len(truth_keys),
            'all_obs_keys': list(adata.obs.columns)  # Also return all obs keys for reference
        })
    
    except Exception as e:
        import traceback
        error_msg = f'Failed to list ground truth: {str(e)}'
        print(f"[GROUND_TRUTH] Error: {error_msg}")
        print(traceback.format_exc())
        return jsonify({'error': error_msg}), 500


@bp.route('/ground-truth/get', methods=['GET', 'POST'])
def get_ground_truth():
    """Get ground truth labels for a specific key
    
    Request:
        GET: ?session_id=xxx&truth_key=xxx
        POST: {"session_id": "xxx", "truth_key": "xxx"}
    
    Returns:
        {
            "success": true,
            "session_id": "string",
            "truth_key": "string",
            "labels": ["label1", "label2", ...],
            "n_labels": int,
            "unique_labels": ["label1", "label2", ...],
            "n_unique": int
        }
    """
    try:
        # Get parameters from query params or JSON body
        if request.method == 'GET':
            session_id = request.args.get('session_id')
            truth_key = request.args.get('truth_key', 'ground_truth')
        else:
            data = request.json or {}
            session_id = data.get('session_id')
            truth_key = data.get('truth_key', 'ground_truth')
        
        if not session_id:
            return jsonify({'error': 'session_id is required'}), 400
        
        if not truth_key:
            return jsonify({'error': 'truth_key is required'}), 400
        
        # Validate session
        is_valid, error_msg = _validate_session(session_id)
        if not is_valid:
            return jsonify({'error': error_msg}), 400
        
        # Load adata
        adata = load_adata(session_id)
        
        # Check if truth_key exists
        if truth_key not in adata.obs:
            return jsonify({
                'error': f'Ground truth key "{truth_key}" not found in adata.obs',
                'available_keys': list(adata.obs.columns)
            }), 404
        
        # Get labels
        labels = adata.obs[truth_key].tolist()
        unique_labels = sorted(adata.obs[truth_key].cat.categories.tolist()) if isinstance(adata.obs[truth_key].dtype, pd.CategoricalDtype) else sorted(list(set(labels)))
        
        return jsonify({
            'success': True,
            'session_id': session_id,
            'truth_key': truth_key,
            'labels': labels,
            'n_labels': len(labels),
            'unique_labels': unique_labels,
            'n_unique': len(unique_labels)
        })
    
    except Exception as e:
        import traceback
        error_msg = f'Failed to get ground truth: {str(e)}'
        print(f"[GROUND_TRUTH] Error: {error_msg}")
        print(traceback.format_exc())
        return jsonify({'error': error_msg}), 500


@bp.route('/ground-truth/delete', methods=['POST'])
def delete_ground_truth():
    """Delete ground truth labels for a specific key
    
    Request:
        {
            "session_id": "string" (required),
            "truth_key": "string" (required)
        }
    
    Returns:
        {
            "success": true,
            "session_id": "string",
            "truth_key": "string",
            "message": "string"
        }
    """
    try:
        data = request.json or {}
        session_id = data.get('session_id')
        truth_key = data.get('truth_key')
        
        if not session_id:
            return jsonify({'error': 'session_id is required'}), 400
        
        if not truth_key:
            return jsonify({'error': 'truth_key is required'}), 400
        
        # Validate session
        is_valid, error_msg = _validate_session(session_id)
        if not is_valid:
            return jsonify({'error': error_msg}), 400
        
        # Load adata
        adata = load_adata(session_id)
        
        # Check if truth_key exists
        if truth_key not in adata.obs:
            return jsonify({
                'error': f'Ground truth key "{truth_key}" not found in adata.obs',
                'available_keys': list(adata.obs.columns)
            }), 404
        
        # Delete the key
        del adata.obs[truth_key]
        save_adata(adata, session_id)
        
        print(f"[GROUND_TRUTH] Deleted ground truth key '{truth_key}' from session {session_id}")
        
        return jsonify({
            'success': True,
            'session_id': session_id,
            'truth_key': truth_key,
            'message': f'Ground truth key "{truth_key}" deleted successfully'
        })
    
    except Exception as e:
        import traceback
        error_msg = f'Failed to delete ground truth: {str(e)}'
        print(f"[GROUND_TRUTH] Error: {error_msg}")
        print(traceback.format_exc())
        return jsonify({'error': error_msg}), 500

