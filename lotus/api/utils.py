"""
Utility functions for Lotus API
"""

from pathlib import Path
from anndata import AnnData
from typing import Optional

from .config import UPLOAD_FOLDER, LOTUS_AVAILABLE


def get_session_dir(session_id: str = 'default') -> Path:
    """Get session directory path"""
    return UPLOAD_FOLDER / session_id


def get_adata_path(session_id: str = 'default') -> Path:
    """Get AnnData file path for a session"""
    return get_session_dir(session_id) / 'data.h5ad'


def load_adata(session_id: str = 'default') -> Optional[AnnData]:
    """Load AnnData from session directory
    
    Automatically restores JSON-serialized layers and core_layers that were
    converted to JSON strings during saving (due to h5py limitations).
    """
    import json
    
    from .config import read, sc
    
    adata_path = get_adata_path(session_id)
    if not adata_path.exists():
        return None
    
    if LOTUS_AVAILABLE:
        adata = read(adata_path)
    else:
        if sc is None:
            raise ImportError('scanpy not available')
        adata = sc.read(adata_path)
    
    # Restore JSON-serialized layers and core_layers
    for key in list(adata.uns.keys()):
        if 'cplearn' in key.lower() and isinstance(adata.uns[key], dict):
            cplearn_info = adata.uns[key]
            
            # Restore 'layers' from JSON string if present
            if 'layers_json' in cplearn_info and 'layers' not in cplearn_info:
                try:
                    layers = json.loads(cplearn_info['layers_json'])
                    cplearn_info['layers'] = layers
                    del cplearn_info['layers_json']
                    print(f"[LOAD] Restored 'layers' from JSON in {key}")
                except (json.JSONDecodeError, TypeError) as e:
                    print(f"[LOAD] Warning: Could not restore layers from JSON in {key}: {e}")
            
            # Restore 'core_layers' from JSON string if present
            if 'core_layers_json' in cplearn_info and 'core_layers' not in cplearn_info:
                try:
                    core_layers = json.loads(cplearn_info['core_layers_json'])
                    cplearn_info['core_layers'] = core_layers
                    del cplearn_info['core_layers_json']
                    print(f"[LOAD] Restored 'core_layers' from JSON in {key}")
                except (json.JSONDecodeError, TypeError) as e:
                    print(f"[LOAD] Warning: Could not restore core_layers from JSON in {key}: {e}")
    
    return adata


def save_adata(adata: AnnData, session_id: str = 'default') -> None:
    """Save AnnData to session directory with memory optimization
    
    Handles serialization of complex data structures in uns that may not be
    directly compatible with h5py (e.g., lists of different lengths).
    h5py cannot save arrays with inhomogeneous shapes, so we convert such
    structures to JSON strings before saving.
    """
    import json
    import numpy as np
    
    # Fix variable-length lists in uns that cause h5py errors
    # Specifically handle 'layers' and 'core_layers' in cplearn entries
    for key in list(adata.uns.keys()):
        if 'cplearn' in key.lower() and isinstance(adata.uns[key], dict):
            cplearn_info = adata.uns[key]
            
            # Handle 'layers' field
            if 'layers' in cplearn_info and isinstance(cplearn_info['layers'], list):
                layers = cplearn_info['layers']
                if len(layers) > 0:
                    # Check if layers have variable lengths
                    try:
                        lengths = [len(layer) if hasattr(layer, '__len__') else 1 for layer in layers]
                        if len(set(lengths)) > 1:
                            # Variable lengths - convert to JSON string
                            print(f"[SAVE] Converting variable-length 'layers' in {key} to JSON string")
                            try:
                                # Ensure all elements are lists (not numpy arrays)
                                serializable = [
                                    list(layer) if isinstance(layer, (list, np.ndarray)) else layer
                                    for layer in layers
                                ]
                                cplearn_info['layers_json'] = json.dumps(serializable)
                                del cplearn_info['layers']
                            except (TypeError, ValueError) as e:
                                print(f"[SAVE] Warning: Could not serialize layers in {key}: {e}")
                                # Remove problematic entry
                                del cplearn_info['layers']
                    except Exception as e:
                        print(f"[SAVE] Warning: Error processing layers in {key}: {e}")
            
            # Handle 'core_layers' field
            if 'core_layers' in cplearn_info and isinstance(cplearn_info['core_layers'], list):
                core_layers = cplearn_info['core_layers']
                if len(core_layers) > 0:
                    try:
                        lengths = [len(layer) if hasattr(layer, '__len__') else 1 for layer in core_layers]
                        if len(set(lengths)) > 1:
                            # Variable lengths - convert to JSON string
                            print(f"[SAVE] Converting variable-length 'core_layers' in {key} to JSON string")
                            try:
                                serializable = [
                                    list(layer) if isinstance(layer, (list, np.ndarray)) else layer
                                    for layer in core_layers
                                ]
                                cplearn_info['core_layers_json'] = json.dumps(serializable)
                                del cplearn_info['core_layers']
                            except (TypeError, ValueError) as e:
                                print(f"[SAVE] Warning: Could not serialize core_layers in {key}: {e}")
                                del cplearn_info['core_layers']
                    except Exception as e:
                        print(f"[SAVE] Warning: Error processing core_layers in {key}: {e}")
            
            # Handle 'flowrank_score' field - h5py requires string keys, not integer keys
            if 'flowrank_score' in cplearn_info and isinstance(cplearn_info['flowrank_score'], dict):
                flowrank_score = cplearn_info['flowrank_score']
                if flowrank_score is not None and len(flowrank_score) > 0:
                    # Check if keys are integers (which h5py can't handle)
                    first_key = next(iter(flowrank_score.keys()))
                    if isinstance(first_key, (int, np.integer)):
                        # Convert integer keys to strings
                        print(f"[SAVE] Converting integer keys in 'flowrank_score' in {key} to strings")
                        cplearn_info['flowrank_score'] = {
                            str(k): v for k, v in flowrank_score.items()
                        }
                    # Also check if values are complex types that need conversion
                    try:
                        # Try to ensure values are serializable
                        for k, v in cplearn_info['flowrank_score'].items():
                            if isinstance(v, (list, np.ndarray)):
                                cplearn_info['flowrank_score'][k] = list(v) if isinstance(v, np.ndarray) else v
                    except Exception as e:
                        print(f"[SAVE] Warning: Error processing flowrank_score values in {key}: {e}")
                        # Convert to JSON string as fallback
                        try:
                            cplearn_info['flowrank_score_json'] = json.dumps(cplearn_info['flowrank_score'])
                            del cplearn_info['flowrank_score']
                        except:
                            # Last resort: remove it
                            print(f"[SAVE] Removing problematic flowrank_score in {key}")
                            del cplearn_info['flowrank_score']
    
    session_dir = get_session_dir(session_id)
    session_dir.mkdir(exist_ok=True, parents=True)
    adata_path = get_adata_path(session_id)
    
    try:
        adata.write(adata_path)
    except (ValueError, TypeError) as e:
        error_str = str(e).lower()
        if "inhomogeneous" in error_str or "layers" in error_str or "flowrank" in error_str or "pathlike" in error_str:
            print(f"[SAVE] Save error detected, attempting additional fix: {e}")
            # Additional fix: remove or convert any remaining problematic entries
            for key in list(adata.uns.keys()):
                if isinstance(adata.uns[key], dict):
                    # Fix layers
                    for subkey in ['layers', 'core_layers']:
                        if subkey in adata.uns[key]:
                            print(f"[SAVE] Removing problematic {key}/{subkey}")
                            del adata.uns[key][subkey]
                    
                    # Fix flowrank_score - convert integer keys to strings or remove
                    if 'flowrank_score' in adata.uns[key] and isinstance(adata.uns[key]['flowrank_score'], dict):
                        flowrank = adata.uns[key]['flowrank_score']
                        if flowrank:
                            try:
                                # Convert integer keys to strings
                                adata.uns[key]['flowrank_score'] = {
                                    str(k): v for k, v in flowrank.items()
                                }
                            except:
                                # If conversion fails, convert to JSON or remove
                                try:
                                    adata.uns[key]['flowrank_score_json'] = json.dumps(flowrank)
                                    del adata.uns[key]['flowrank_score']
                                except:
                                    print(f"[SAVE] Removing problematic {key}/flowrank_score")
                                    del adata.uns[key]['flowrank_score']
            
            # Try saving again
            try:
                adata.write(adata_path)
            except Exception as e2:
                print(f"[SAVE] Failed to save after fix attempt: {e2}")
                raise
        else:
            raise
    
    # Force garbage collection after saving to free memory
    import gc
    gc.collect()

