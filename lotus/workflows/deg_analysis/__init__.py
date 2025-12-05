from __future__ import annotations

from typing import Any, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData

from lotus.methods.cplearn import cplearn
from lotus.methods.scanpy import tl as sc_tl
from lotus.methods.scanpy import get as sc_get


def _marker_genes_scanpy(
    adata: AnnData,
    *,
    cluster_key: str,
    groups_a: set[int] | set[str] | None = None,
    groups_b: set[int] | set[str] | None = None,
    layer: str | None = None,
    auto_pick_groups: bool = True,
    scanpy_method: str | None = None,
    use_raw: bool | None = None,
) -> pd.DataFrame:
    """
    Identify marker genes using scanpy's rank_genes_groups method.
    
    Parameters:
        adata: AnnData object
        cluster_key: Key name for cluster labels in adata.obs
        groups_a: Set of cluster labels for first group
        groups_b: Set of cluster labels for second group
        layer: Layer to use for analysis
        auto_pick_groups: Auto-select first two non-negative clusters
        scanpy_method: Statistical method (wilcoxon, t-test, logreg, etc.)
        use_raw: Whether to use raw data
        
    Returns:
        pd.DataFrame: Differential expression analysis results
    """
    # Get cluster labels
    labels_raw = adata.obs[cluster_key]
    
    # Ensure labels are in Categorical format with string categories (scanpy standard)
    # This is important for cplearn clusters which might be integers
    # scanpy works best with string Categorical categories
    def get_categories_safe(series: pd.Series) -> list:
        """Safely get categories from a pandas Series."""
        if isinstance(series.dtype, pd.CategoricalDtype):
            try:
                return list(series.cat.categories)
            except (AttributeError, TypeError):
                try:
                    return list(series.dtype.categories)
                except (AttributeError, TypeError):
                    return sorted(set(series.unique()))
        else:
            return sorted(set(series.unique()))
    
    # Always ensure labels are in string Categorical format (scanpy standard)
    # This is important for cplearn clusters which might be integers
    # scanpy works best with string Categorical categories
    original_dtype = labels_raw.dtype
    print(f"[DEG] Original label dtype: {original_dtype}, Is Categorical: {isinstance(original_dtype, pd.CategoricalDtype)}")
    
    # Check if we need to convert to string Categorical
    needs_conversion = False
    if not isinstance(labels_raw.dtype, pd.CategoricalDtype):
        # Not Categorical, needs conversion
        needs_conversion = True
        print(f"[DEG] Not Categorical, will convert to string Categorical")
    else:
        # Already Categorical, check if categories are all strings
        categories = get_categories_safe(labels_raw)
        print(f"[DEG] Already Categorical, categories: {categories}, types: {[type(c).__name__ for c in categories]}")
        if not categories or not all(isinstance(c, str) for c in categories):
            # Categories are not all strings (e.g., integers), need conversion
            needs_conversion = True
            print(f"[DEG] Categories are not all strings, will convert to string Categorical")
        else:
            print(f"[DEG] Categories are already strings, no conversion needed")
    
    if needs_conversion:
        # Convert to string first, then to Categorical
        labels_str = [str(x) for x in labels_raw]
        # Create temporary Categorical column with string categories
        temp_key = f"{cluster_key}_scanpy_temp"
        adata.obs[temp_key] = pd.Categorical(labels_str)
        cluster_key_to_use = temp_key
        labels_raw = adata.obs[temp_key]
        print(f"[DEG] Created temporary Categorical column '{temp_key}' with string categories: {list(labels_raw.cat.categories)}")
    else:
        cluster_key_to_use = cluster_key
        print(f"[DEG] Using original cluster key '{cluster_key}' with categories: {list(labels_raw.cat.categories)}")
    
    # Convert labels to strings (scanpy uses string labels)
    labels_str = [str(x) for x in labels_raw]
    
    # Filter out invalid labels (e.g., -1, NaN, empty strings)
    unique_labels = sorted(set([l for l in labels_str if l and str(l).strip() and str(l) != '-1' and str(l).lower() != 'nan']))
    
    # Auto-pick groups if needed
    if auto_pick_groups and (groups_a is None or groups_b is None):
        # Pick first two groups
        if len(unique_labels) < 2:
            raise ValueError(
                f"At least two valid clusters are required for DEG analysis. "
                f"Found {len(unique_labels)} valid cluster(s): {unique_labels}. "
                f"Total unique labels (including invalid): {len(set(labels_str))}"
            )
        if len(unique_labels) == 0:
            raise ValueError("No valid clusters found in cluster_key. Please check your clustering results.")
        groups_a = {unique_labels[0]}
        groups_b = {unique_labels[1]}
    
    if groups_a is None or groups_b is None:
        raise ValueError("Must specify groups_a and groups_b, or set auto_pick_groups=True")
    
    # Convert to string sets for scanpy
    groups_a_str = {str(g) for g in groups_a}
    groups_b_str = {str(g) for g in groups_b}
    
    # Validate groups are not empty
    if not groups_a_str:
        raise ValueError("groups_a cannot be empty after conversion to strings")
    if not groups_b_str:
        raise ValueError("groups_b cannot be empty after conversion to strings")
    
    # Get the first group for comparison
    group_a_str = list(groups_a_str)[0]
    
    # Get the actual Categorical values to ensure correct matching
    # Categorical might have integer categories but string representations
    # For pandas Series, safely access categories
    def get_categories(series: pd.Series) -> list:
        """Safely get categories from a pandas Series."""
        if isinstance(series.dtype, pd.CategoricalDtype):
            try:
                # Try accessing through .cat accessor (standard way)
                return list(series.cat.categories)
            except (AttributeError, TypeError):
                try:
                    # Fallback to dtype.categories
                    return list(series.dtype.categories)
                except (AttributeError, TypeError):
                    # Final fallback: use unique values
                    return sorted(set(series.unique()))
        else:
            # Not categorical, just return unique values
            return sorted(set(series.unique()))
    
    categorical_categories = get_categories(labels_raw)
    categorical_str_values = [str(cat) for cat in categorical_categories]
    
    # Find the matching group in Categorical categories
    # After conversion, categories should all be strings
    # Use string type for group name to match scanpy's expectation
    group_a = None
    for cat_val in categorical_categories:
        if str(cat_val) == group_a_str:
            # Use string representation to match scanpy's expectation
            group_a = str(cat_val)
            break
    
    # If not found in categories, use the string value directly
    if group_a is None:
        # Check if group_a_str matches any category by string comparison
        for cat_val in categorical_categories:
            if str(cat_val) == group_a_str:
                group_a = str(cat_val)
                break
        if group_a is None:
            # Use the string value directly (should match after conversion)
            group_a = group_a_str
    
    # Validate that group_a exists in the cluster labels
    # Check both string and actual values
    if group_a not in unique_labels and str(group_a) not in unique_labels:
        # Try to find a match by comparing string representations
        found_match = False
        for ul in unique_labels:
            if str(ul) == str(group_a) or str(ul) == group_a:
                group_a = str(ul)
                found_match = True
                break
        if not found_match:
            # Try matching with Categorical categories
            for cat_val in categorical_categories:
                if str(cat_val) in unique_labels or str(cat_val) == group_a_str:
                    group_a = str(cat_val)
                    found_match = True
                    break
        if not found_match:
            raise ValueError(
                f"Group '{group_a}' (from '{group_a_str}') not found in cluster_key '{cluster_key}'. "
                f"Available groups: {unique_labels}, Categorical categories: {categorical_categories}"
            )
    
    print(f"[DEG] Selected group_a: '{group_a}' (from unique_labels: {unique_labels})")
    
    # Auto-detect layer or use_raw
    # Important: When layer is specified, use_raw should be False
    if layer is not None:
        use_raw = False  # Cannot use both layer and use_raw
        print(f"[DEG] Using layer '{layer}' for analysis (use_raw=False)")
    elif use_raw is None:
        # Check if raw data exists
        use_raw = adata.raw is not None
        if use_raw:
            print(f"[DEG] Using raw data (adata.raw)")
        else:
            print(f"[DEG] Using adata.X (no layer or raw data)")
    
    # Debug: Print group information before calling rank_genes_groups
    print(f"[DEG] Calling rank_genes_groups with group_a={group_a!r} (type: {type(group_a).__name__}), groupby='{cluster_key_to_use}'")
    print(f"[DEG] Unique labels in data: {unique_labels}")
    # Safely print categories information
    print(f"[DEG] Label dtype: {labels_raw.dtype}, categories: {categorical_categories}")
    print(f"[DEG] Category types: {[type(c).__name__ for c in categorical_categories]}")
    # Check actual values in the Categorical
    sample_values = labels_raw.head(5)
    print(f"[DEG] Sample label values: {sample_values.tolist()}, types: {[type(v).__name__ for v in sample_values]}")
    print(f"[DEG] Layer: {layer}, use_raw: {use_raw}")
    
    # Run scanpy's rank_genes_groups
    # Compare group_a vs rest (which includes group_b)
    key_added = "rank_genes_groups"
    
    # Try multiple strategies if first attempt fails
    # Note: scanpy's rank_genes_groups expects log-transformed data
    # Priority: use_raw (usually log-transformed) > adata.X (may be log-transformed) > layer (may be raw counts)
    strategies = []
    
    # Strategy 1: Use raw data if available (usually contains log-transformed data)
    # This is preferred because adata.raw typically contains normalized+log1p data
    if adata.raw is not None:
        strategies.append({
            'layer': None,
            'use_raw': True,
            'name': 'use_raw=True (log-transformed)'
        })
    
    # Strategy 2: Use adata.X (may already be log-transformed)
    strategies.append({
        'layer': None,
        'use_raw': False,
        'name': 'adata.X'
    })
    
    # Strategy 3: Use layer if provided (may be raw counts - will trigger warning)
    # Only use layer as fallback since it may contain raw counts
    if layer is not None:
        strategies.append({
            'layer': layer,
            'use_raw': False,
            'name': f'layer={layer} (may trigger log warning)'
        })
    
    # Special handling for logistic regression:
    # - Logistic regression with groups=[group_a] fails with "single cluster" error
    # - Must use groups='all' for logistic regression, but this may not provide pvals/logfoldchanges
    # - We'll try groups='all' for logistic regression, and accept missing stats if necessary
    is_logreg = (scanpy_method or "wilcoxon") == "logreg"
    
    last_error = None
    for strategy_idx, strategy in enumerate(strategies):
        try:
            print(f"[DEG] Trying strategy {strategy_idx + 1}/{len(strategies)}: {strategy['name']}")
            
            # Prepare rank_genes_groups parameters
            rank_params = {
                'groupby': cluster_key_to_use,  # Use the compatible cluster key
                'method': scanpy_method or "wilcoxon",
                'use_raw': strategy['use_raw'],
                'layer': strategy['layer'],
                'key_added': key_added,
                'pts': True,  # Include percentage of cells expressing gene
            }
            
            # For logistic regression, use groups='all' (required by scanpy)
            # For other methods, use groups=[group_a] with reference='rest'
            if is_logreg:
                rank_params['groups'] = 'all'
                print(f"[DEG] Using groups='all' for logistic regression method (required)")
            else:
                rank_params['groups'] = [group_a]  # Use the matched group name
                rank_params['reference'] = "rest"
                print(f"[DEG] Using groups=[{group_a}], reference='rest' for {scanpy_method or 'wilcoxon'} method")
            
            rank_genes_groups(adata, **rank_params)
            print(f"[DEG] Success with strategy: {strategy['name']}")
            break  # Success, exit loop
        except Exception as e:
            error_msg = str(e)
            last_error = e
            print(f"[DEG] Strategy {strategy_idx + 1} failed: {error_msg}")
            
            # Check if this is a fatal error (not just data format issue)
            if "index" in error_msg.lower() and ("out of bounds" in error_msg.lower() or "size 0" in error_msg.lower()):
                # Try next strategy if available
                if strategy_idx < len(strategies) - 1:
                    print(f"[DEG] Index error, trying next strategy...")
                    continue
                else:
                    # Last strategy failed, return empty result
                    print(f"[DEG] All strategies failed. Last error: {error_msg}")
                    if cluster_key_to_use != cluster_key and cluster_key_to_use in adata.obs:
                        del adata.obs[cluster_key_to_use]
                    return pd.DataFrame(columns=['gene', 'log2fc', 'pvalue', 'p_adj', 'z_score', 'mean_a', 'mean_b', 'pct_expr_a', 'pct_expr_b'])
            elif "single cluster" in error_msg.lower() or "cannot perform logistic regression" in error_msg.lower():
                # This should not happen with groups='all', but handle it anyway
                print(f"[DEG] Logistic regression error: {error_msg}")
                if cluster_key_to_use != cluster_key and cluster_key_to_use in adata.obs:
                    del adata.obs[cluster_key_to_use]
                raise ValueError(
                    f"Logistic regression method requires at least 2 clusters with sufficient cells. "
                    f"Found {len(unique_labels)} unique clusters: {unique_labels}. "
                    f"Please try a different statistical method (e.g., 'wilcoxon' or 't-test') or ensure you have at least 2 clusters with enough cells."
                ) from e
            else:
                # Other errors, try next strategy
                if strategy_idx < len(strategies) - 1:
                    continue
                else:
                    # All strategies exhausted
                    raise ValueError(f"Marker genes analysis failed with all strategies. Last error: {error_msg}") from e
    else:
        # All strategies failed
        if last_error is not None:
            error_msg = str(last_error)
            if cluster_key_to_use != cluster_key and cluster_key_to_use in adata.obs:
                del adata.obs[cluster_key_to_use]
            if "index" in error_msg.lower() and ("out of bounds" in error_msg.lower() or "size 0" in error_msg.lower()):
                print(f"[DEG] All strategies failed with index error: {error_msg}")
                return pd.DataFrame(columns=['gene', 'log2fc', 'pvalue', 'p_adj', 'z_score', 'mean_a', 'mean_b', 'pct_expr_a', 'pct_expr_b'])
            raise ValueError(f"Marker genes analysis failed: {error_msg}") from last_error
    
    # Check if results were actually stored
    if key_added not in adata.uns:
        # Clean up temporary cluster key if created
        if cluster_key_to_use != cluster_key and cluster_key_to_use in adata.obs:
            del adata.obs[cluster_key_to_use]
        print(f"[DEG] Warning: rank_genes_groups did not store results in adata.uns['{key_added}']")
        return pd.DataFrame(columns=['gene', 'log2fc', 'pvalue', 'p_adj', 'z_score', 'mean_a', 'mean_b', 'pct_expr_a', 'pct_expr_b'])
    
    # Check if results are empty before trying to extract
    rank_result = adata.uns[key_added]
    if 'names' not in rank_result:
        # Clean up temporary cluster key if created
        if cluster_key_to_use != cluster_key and cluster_key_to_use in adata.obs:
            del adata.obs[cluster_key_to_use]
        print(f"[DEG] Warning: rank_genes_groups results missing 'names' field")
        return pd.DataFrame(columns=['gene', 'log2fc', 'pvalue', 'p_adj', 'z_score', 'mean_a', 'mean_b', 'pct_expr_a', 'pct_expr_b'])
    
    # Check if results are actually empty
    names = rank_result['names']
    # Check if names is empty
    if hasattr(names, 'size') and names.size == 0:
        # Clean up temporary cluster key if created
        if cluster_key_to_use != cluster_key and cluster_key_to_use in adata.obs:
            del adata.obs[cluster_key_to_use]
        print(f"[DEG] Warning: rank_genes_groups results are empty (no genes found)")
        return pd.DataFrame(columns=['gene', 'log2fc', 'pvalue', 'p_adj', 'z_score', 'mean_a', 'mean_b', 'pct_expr_a', 'pct_expr_b'])
    elif hasattr(names, '__len__') and len(names) == 0:
        # Clean up temporary cluster key if created
        if cluster_key_to_use != cluster_key and cluster_key_to_use in adata.obs:
            del adata.obs[cluster_key_to_use]
        print(f"[DEG] Warning: rank_genes_groups results are empty (no genes found)")
        return pd.DataFrame(columns=['gene', 'log2fc', 'pvalue', 'p_adj', 'z_score', 'mean_a', 'mean_b', 'pct_expr_a', 'pct_expr_b'])
    
    # Convert results to DataFrame
    try:
        result = _convert_scanpy_rank_genes_groups_to_dataframe(
            adata,
            cluster_key=cluster_key,
            groups_a=groups_a_str,
            groups_b=groups_b_str,
            key=key_added,
        )
    except Exception as e:
        error_msg = str(e)
        # If conversion fails, return empty DataFrame instead of raising error
        print(f"[DEG] Warning: Failed to convert rank_genes_groups results to DataFrame: {error_msg}")
        # Clean up temporary cluster key if created
        if cluster_key_to_use != cluster_key and cluster_key_to_use in adata.obs:
            del adata.obs[cluster_key_to_use]
        return pd.DataFrame(columns=['gene', 'log2fc', 'pvalue', 'p_adj', 'z_score', 'mean_a', 'mean_b', 'pct_expr_a', 'pct_expr_b'])
    finally:
        # Clean up temporary cluster key if created
        if cluster_key_to_use != cluster_key and cluster_key_to_use in adata.obs:
            del adata.obs[cluster_key_to_use]
    
    return result


def pick_groups(labels: Sequence[int]) -> tuple[set[int], set[int]]:
    """
    Select two groups for differential expression analysis
    
    Parameters:
        labels: Sequence of cluster labels
    
    Returns:
        Two sets of group labels (groups_a, groups_b)
    """
    unique = [label for label in sorted(set(labels)) if label != -1]
    if len(unique) < 2:
        raise ValueError("At least two non-negative cluster labels are required for differential expression analysis")
    return {unique[0]}, {unique[1]}


def _convert_scanpy_rank_genes_groups_to_dataframe(
    adata: AnnData,
    cluster_key: str,
    groups_a: set[int] | set[str],
    groups_b: set[int] | set[str],
    key: str = "rank_genes_groups",
) -> pd.DataFrame:
    """
    Convert scanpy rank_genes_groups results to unified DataFrame format.
    
    This function extracts results from adata.uns[key] and converts them to
    a DataFrame format compatible with cplearn.de_from_adata output.
    
    Parameters:
        adata: AnnData object with rank_genes_groups results
        cluster_key: Key name for cluster labels
        groups_a: Set of group labels for first group
        groups_b: Set of group labels for second group
        key: Key in adata.uns containing rank_genes_groups results
        
    Returns:
        pd.DataFrame: Unified format DataFrame with columns: gene, log2fc, pvalue, p_adj, etc.
    """
    if key not in adata.uns:
        raise KeyError(f"rank_genes_groups results not found in adata.uns['{key}']")
    
    # Convert group sets to strings for scanpy (it uses string labels)
    groups_a_str = {str(g) for g in groups_a}
    
    # Get unique groups (take first from each set)
    group_a = list(groups_a_str)[0] if groups_a_str else None
    
    if group_a is None:
        raise ValueError("groups_a must contain at least one element")
    
    # First, validate the results structure and find the correct group name
    rank_result = adata.uns[key]
    
    # Check if results are empty before proceeding
    if 'names' not in rank_result:
        print(f"[DEG] Warning: rank_genes_groups results missing 'names' field")
        return pd.DataFrame(columns=['gene', 'log2fc', 'pvalue', 'p_adj', 'z_score', 'mean_a', 'mean_b', 'pct_expr_a', 'pct_expr_b'])
    
    names = rank_result['names']
    
    # Check if names array is empty
    if hasattr(names, 'size') and names.size == 0:
        print(f"[DEG] Warning: rank_genes_groups results are empty (names.size == 0)")
        return pd.DataFrame(columns=['gene', 'log2fc', 'pvalue', 'p_adj', 'z_score', 'mean_a', 'mean_b', 'pct_expr_a', 'pct_expr_b'])
    elif hasattr(names, '__len__') and len(names) == 0:
        print(f"[DEG] Warning: rank_genes_groups results are empty (len(names) == 0)")
        return pd.DataFrame(columns=['gene', 'log2fc', 'pvalue', 'p_adj', 'z_score', 'mean_a', 'mean_b', 'pct_expr_a', 'pct_expr_b'])
    
    # Check available groups in the results
    available_groups = []
    # scanpy stores results as structured arrays with group names as field names
    if hasattr(names, 'dtype') and hasattr(names.dtype, 'names') and names.dtype.names:
        available_groups = list(names.dtype.names)
    elif isinstance(names, dict):
        available_groups = list(names.keys())
    
    # If no groups found, results might be empty
    if not available_groups:
        print(f"[DEG] Warning: No groups found in rank_genes_groups results")
        return pd.DataFrame(columns=['gene', 'log2fc', 'pvalue', 'p_adj', 'z_score', 'mean_a', 'mean_b', 'pct_expr_a', 'pct_expr_b'])
    
    # Try to match group_a with available groups
    group_to_use = group_a
    if group_a not in available_groups:
        # Try string matching
        for ag in available_groups:
            if str(ag) == str(group_a):
                group_to_use = ag
                break
        # If still not found, use first available group
        if group_to_use not in available_groups:
            print(f"[DEG] Warning: Group '{group_a}' not found in available groups {available_groups}, using first available group")
            group_to_use = available_groups[0]
    
    # Extract results directly from adata.uns['rank_genes_groups']
    # rank_genes_groups_df may not return all columns for logistic regression
    print(f"[DEG] Extracting results directly from adata.uns['{key}'] for group '{group_to_use}'")
    try:
        rank_result = adata.uns[key]
        
        # Check the structure of rank_result
        print(f"[DEG] rank_result keys: {list(rank_result.keys())}")
        
        # Handle different data structures
        # scanpy stores results as structured arrays or dicts
        if isinstance(rank_result['names'], np.ndarray):
            # Structured array format - access by field name
            if hasattr(rank_result['names'].dtype, 'names') and group_to_use in rank_result['names'].dtype.names:
                names = rank_result['names'][group_to_use]
            else:
                # Try accessing as structured array
                names = rank_result['names'][group_to_use] if group_to_use in rank_result['names'].dtype.names else rank_result['names']
        elif isinstance(rank_result['names'], dict):
            # Dict format - access by group name
            names = rank_result['names'][group_to_use]
        else:
            # Fallback: try direct access
            names = rank_result['names']
        
        # Get scores
        if 'scores' in rank_result:
            if isinstance(rank_result['scores'], np.ndarray):
                if hasattr(rank_result['scores'].dtype, 'names') and group_to_use in rank_result['scores'].dtype.names:
                    scores = rank_result['scores'][group_to_use]
                else:
                    scores = rank_result['scores']
            elif isinstance(rank_result['scores'], dict):
                scores = rank_result['scores'].get(group_to_use, np.zeros(len(names)))
            else:
                scores = rank_result['scores']
        else:
            scores = np.zeros(len(names))
        
        # Check if this is logistic regression with groups='all' (which doesn't provide pvals/logfoldchanges)
        # We'll use placeholders instead of default values
        use_placeholders = 'logfoldchanges' not in rank_result or 'pvals' not in rank_result
        
        # Get logfoldchanges
        if 'logfoldchanges' in rank_result:
            if isinstance(rank_result['logfoldchanges'], np.ndarray):
                if hasattr(rank_result['logfoldchanges'].dtype, 'names') and group_to_use in rank_result['logfoldchanges'].dtype.names:
                    logfoldchanges = rank_result['logfoldchanges'][group_to_use]
                else:
                    logfoldchanges = rank_result['logfoldchanges']
            elif isinstance(rank_result['logfoldchanges'], dict):
                logfoldchanges = rank_result['logfoldchanges'].get(group_to_use, np.zeros(len(names)))
            else:
                logfoldchanges = rank_result['logfoldchanges']
        else:
            logfoldchanges = None  # Will use placeholder
            print(f"[DEG] Warning: 'logfoldchanges' not found in rank_result, will use placeholder")
        
        # Get pvals
        if 'pvals' in rank_result:
            if isinstance(rank_result['pvals'], np.ndarray):
                if hasattr(rank_result['pvals'].dtype, 'names') and group_to_use in rank_result['pvals'].dtype.names:
                    pvals = rank_result['pvals'][group_to_use]
                else:
                    pvals = rank_result['pvals']
            elif isinstance(rank_result['pvals'], dict):
                pvals = rank_result['pvals'].get(group_to_use, np.ones(len(names)))
            else:
                pvals = rank_result['pvals']
        else:
            pvals = None  # Will use placeholder
            print(f"[DEG] Warning: 'pvals' not found in rank_result, will use placeholder")
        
        # Get pvals_adj
        if 'pvals_adj' in rank_result:
            if isinstance(rank_result['pvals_adj'], np.ndarray):
                if hasattr(rank_result['pvals_adj'].dtype, 'names') and group_to_use in rank_result['pvals_adj'].dtype.names:
                    pvals_adj = rank_result['pvals_adj'][group_to_use]
                else:
                    pvals_adj = rank_result['pvals_adj']
            elif isinstance(rank_result['pvals_adj'], dict):
                pvals_adj = rank_result['pvals_adj'].get(group_to_use, np.ones(len(names)))
            else:
                pvals_adj = rank_result['pvals_adj']
        else:
            pvals_adj = None  # Will use placeholder
            print(f"[DEG] Warning: 'pvals_adj' not found in rank_result, will use placeholder")
        
        # Ensure all arrays have the same length
        min_len = len(names)
        names = names[:min_len]
        scores = scores[:min_len] if len(scores) >= min_len else np.zeros(min_len)
        
        # Use placeholders for missing values (e.g., logistic regression with groups='all')
        if logfoldchanges is None:
            logfoldchanges = np.array(['N/A'] * min_len, dtype=object)
        else:
            logfoldchanges = logfoldchanges[:min_len] if len(logfoldchanges) >= min_len else np.array(['N/A'] * min_len, dtype=object)
        
        if pvals is None:
            pvals = np.array(['N/A'] * min_len, dtype=object)
        else:
            pvals = pvals[:min_len] if len(pvals) >= min_len else np.array(['N/A'] * min_len, dtype=object)
        
        if pvals_adj is None:
            pvals_adj = np.array(['N/A'] * min_len, dtype=object)
        else:
            pvals_adj = pvals_adj[:min_len] if len(pvals_adj) >= min_len else np.array(['N/A'] * min_len, dtype=object)
        
        # Create DataFrame
        df_result = pd.DataFrame({
            'names': names,
            'scores': scores,
            'logfoldchanges': logfoldchanges,
            'pvals': pvals,
            'pvals_adj': pvals_adj,
        })
        print(f"[DEG] Direct extraction: {len(df_result)} rows, columns: {list(df_result.columns)}")
        if len(df_result) > 0:
            first_row = df_result.iloc[0]
            log2fc_val = first_row['logfoldchanges']
            pval_val = first_row['pvals']
            pval_adj_val = first_row['pvals_adj']
            log2fc_str = str(log2fc_val) if log2fc_val == 'N/A' else f"{log2fc_val:.4f}"
            pval_str = str(pval_val) if pval_val == 'N/A' else f"{pval_val:.4f}"
            pval_adj_str = str(pval_adj_val) if pval_adj_val == 'N/A' else f"{pval_adj_val:.4f}"
            print(f"[DEG] First row: gene={first_row['names']}, log2fc={log2fc_str}, pval={pval_str}, pval_adj={pval_adj_str}")
    except Exception as e:
        error_msg = str(e)
        print(f"[DEG] Error extracting directly from adata.uns: {error_msg}")
        import traceback
        traceback.print_exc()
        return pd.DataFrame(columns=['gene', 'log2fc', 'pvalue', 'p_adj', 'z_score', 'mean_a', 'mean_b', 'pct_expr_a', 'pct_expr_b'])
    
    # Check if result is empty or None
    if df_result is None or len(df_result) == 0:
        return pd.DataFrame(columns=['gene', 'log2fc', 'pvalue', 'p_adj', 'z_score', 'mean_a', 'mean_b', 'pct_expr_a', 'pct_expr_b'])
    
    # Format the DataFrame to match expected format
    return _format_scanpy_dataframe(df_result, group_a)


def _format_scanpy_dataframe(df_result: pd.DataFrame, group_a: str) -> pd.DataFrame:
    """
    Format scanpy rank_genes_groups_df result to unified format.
    
    Parameters:
        df_result: DataFrame from rank_genes_groups_df
        group_a: Group name (for reference)
        
    Returns:
        pd.DataFrame: Formatted DataFrame with expected columns
    """
    result_list = []
    
    # Check if DataFrame is empty
    if df_result is None or len(df_result) == 0:
        return pd.DataFrame(columns=['gene', 'log2fc', 'pvalue', 'p_adj', 'z_score', 'mean_a', 'mean_b', 'pct_expr_a', 'pct_expr_b'])
    
    # Extract data from DataFrame
    # Note: rank_genes_groups_df column names may vary by method
    # Common columns: 'names', 'logfoldchanges', 'pvals', 'pvals_adj', 'scores'
    for idx, row in df_result.iterrows():
        try:
            # Get gene name - try different possible column names
            gene = None
            for gene_col in ['names', 'gene', 'gene_name', 'gene_id']:
                if gene_col in row.index:
                    gene = str(row[gene_col])
                    break
            if gene is None:
                gene = str(idx)
            
            # Get log2fc - try different possible column names
            # Use placeholder if value is 'N/A' or missing
            log2fc = 'N/A'
            for fc_col in ['logfoldchanges', 'log2fc', 'log_fold_change', 'fold_change']:
                if fc_col in row.index:
                    val = row[fc_col]
                    if val != 'N/A' and val is not None:
                        try:
                            log2fc = float(val)
                        except (ValueError, TypeError):
                            log2fc = 'N/A'
                    break
            
            # Get pvalue - try different possible column names
            # Use placeholder if value is 'N/A' or missing
            pvalue = 'N/A'
            for pval_col in ['pvals', 'pvalue', 'p_value', 'pval']:
                if pval_col in row.index:
                    val = row[pval_col]
                    if val != 'N/A' and val is not None:
                        try:
                            pvalue = float(val)
                        except (ValueError, TypeError):
                            pvalue = 'N/A'
                    break
            
            # Get p_adj - try different possible column names
            # Use placeholder if value is 'N/A' or missing
            p_adj = 'N/A'
            for padj_col in ['pvals_adj', 'p_adj', 'pval_adj', 'padj', 'p_adjusted']:
                if padj_col in row.index:
                    val = row[padj_col]
                    if val != 'N/A' and val is not None:
                        try:
                            p_adj = float(val)
                        except (ValueError, TypeError):
                            p_adj = 'N/A'
                    break
            
            # Get score - try different possible column names
            score = 0.0
            for score_col in ['scores', 'score', 'z_score', 'statistic']:
                if score_col in row.index:
                    score = float(row[score_col])
                    break
            
            result_list.append({
                'gene': gene,
                'log2fc': log2fc,
                'pvalue': pvalue,
                'p_adj': p_adj,
                'z_score': score,
                'mean_a': 0.0,  # Not directly available from rank_genes_groups
                'mean_b': 0.0,  # Not directly available from rank_genes_groups
                'pct_expr_a': 0.0,  # Could be extracted from pts if needed
                'pct_expr_b': 0.0,
            })
        except (ValueError, TypeError, KeyError) as e:
            # Skip rows with invalid data
            print(f"[DEG] Warning: Skipping row {idx} due to error: {e}")
            continue
    
    if not result_list:
        return pd.DataFrame(columns=['gene', 'log2fc', 'pvalue', 'p_adj', 'z_score', 'mean_a', 'mean_b', 'pct_expr_a', 'pct_expr_b'])
    
    result_df = pd.DataFrame(result_list)
    
    # Sort by p_adj and log2fc
    # Handle placeholders ('N/A') in sorting - treat 'N/A' as high values for p_adj, low values for log2fc
    if len(result_df) > 0:
        def safe_float_sort(val, default_for_na=float('inf')):
            """Convert value to float for sorting, handling 'N/A' placeholder."""
            if val == 'N/A' or val is None:
                return default_for_na
            try:
                return float(val)
            except (ValueError, TypeError):
                return default_for_na
        
        # Create temporary sort keys
        sort_p_adj = result_df['p_adj'].apply(lambda x: safe_float_sort(x, float('inf')))
        sort_log2fc = result_df['log2fc'].apply(lambda x: safe_float_sort(x, float('-inf')))
        
        # Sort by p_adj first (ascending), then by log2fc (descending)
        # Use numpy argsort for stable sorting
        sort_indices = np.lexsort([-sort_log2fc.values, sort_p_adj.values])
        result_df = result_df.iloc[sort_indices]
    
    return result_df


def marker_genes(
    adata: AnnData,
    *,
    cluster_key: str | None = None,
    groups_a: set[int] | None = None,
    groups_b: set[int] | None = None,
    layer: str | None = None,
    min_detect_pct: float = 0.0,
    min_cells_per_group: int = 5,
    auto_pick_groups: bool = True,
    method: str = "cplearn",
    scanpy_method: str | None = None,
    use_raw: bool | None = None,
) -> pd.DataFrame:
    """
    Identify differentially expressed genes (marker genes) between groups.
    
    Parameters:
        adata (AnnData): AnnData object
        cluster_key (str | None): Key name for cluster labels in adata.obs. Default: None (auto-detect)
        groups_a (set[int] | None): Set of cluster labels for first group. Default: None
        groups_b (set[int] | None): Set of cluster labels for second group. Default: None
        layer (str | None): Layer to use for analysis. Default: None (auto-detect)
        min_detect_pct (float): Minimum detection percentage. Default: 0.0 (only used with cplearn method)
        min_cells_per_group (int): Minimum number of cells per group. Default: 5 (only used with cplearn method)
        auto_pick_groups (bool): Auto-select first two non-negative clusters. Default: True
        method (str): DEG analysis method. Options: "cplearn" (default), "scanpy"
        scanpy_method (str | None): Statistical method for scanpy. Options: "wilcoxon" (default), "t-test", "logreg", etc.
        use_raw (bool | None): Whether to use raw data (scanpy method only). Default: None (auto-detect)
    
    Returns:
        pd.DataFrame: Differential expression analysis results
    """
    # Auto-detect cluster key if not specified
    if cluster_key is None:
        for key in ["cplearn", "leiden", "louvain"]:
            if key in adata.obs:
                cluster_key = key
                break
        if cluster_key is None:
            raise ValueError(
                "No cluster key found. Please specify cluster_key or run clustering first. "
                "Compatible keys: 'cplearn', 'leiden', 'louvain'"
            )
    
    # Handle scanpy method
    if method == "scanpy":
        return _marker_genes_scanpy(
            adata,
            cluster_key=cluster_key,
            groups_a=groups_a,
            groups_b=groups_b,
            layer=layer,
            auto_pick_groups=auto_pick_groups,
            scanpy_method=scanpy_method,
            use_raw=use_raw,
        )
    
    # Auto-detect layer if not specified
    if layer is None:
        for layer_name in ["raw_counts", "raw", "counts"]:
            if layer_name in adata.layers:
                layer = layer_name
                break
        # If no layer found, use None (will use adata.X)
        if layer is None:
            layer = None
    
    # Handle categorical cluster labels (scanpy format)
    labels_raw = adata.obs[cluster_key]
    if isinstance(labels_raw.dtype, pd.CategoricalDtype):
        # Convert categorical to numeric for comparison
        labels = np.asarray([int(x) if str(x).isdigit() else -1 for x in labels_raw])
    else:
        labels = np.asarray(labels_raw.astype(int))
    
    if auto_pick_groups and (groups_a is None or groups_b is None):
        groups_a, groups_b = pick_groups(labels)
    
    if groups_a is None or groups_b is None:
        raise ValueError("Must specify groups_a and groups_b, or set auto_pick_groups=True")
    
    groups_a = {int(x) for x in groups_a}
    groups_b = {int(x) for x in groups_b}
    print(f"Comparing groups {groups_a} vs {groups_b}")
    
    # Check if cplearn layers info exists, if not create default
    layers_key = f"{cluster_key}_cplearn"
    if layers_key not in adata.uns:
        # Create default layers structure (all cells in one layer)
        # This allows DEG analysis to work even without cplearn clustering
        print(f"[DEG] Creating default layers structure for {layers_key}")
        adata.uns[layers_key] = {
            'layers': [[i for i in range(adata.n_obs)]]
        }
    
    # For scanpy methods (leiden, louvain), labels are strings ('0', '1', etc.)
    # but de_from_adata expects integer labels. Create a temporary integer label column
    temp_key = f"{cluster_key}_int"
    needs_update = False
    
    if temp_key not in adata.obs:
        needs_update = True
    elif adata.obs[temp_key].dtype != 'int64':
        needs_update = True
    else:
        # Check if the integer labels match the current string labels
        # Convert current labels to integers and compare
        current_labels_int = np.asarray([int(x) if str(x).isdigit() else -1 for x in adata.obs[cluster_key]])
        stored_labels_int = adata.obs[temp_key].values
        
        # Check if they match (same length and same values)
        if len(current_labels_int) != len(stored_labels_int) or not np.array_equal(current_labels_int, stored_labels_int):
            needs_update = True
    
    if needs_update:
        # Convert labels to integers and store temporarily
        adata.obs[temp_key] = labels
        use_key = temp_key
    else:
        use_key = temp_key  # Use existing integer labels
    
    options = cplearn.DEOptions(
        min_detect_pct=min_detect_pct,
        min_cells_per_group=min_cells_per_group,
    )
    
    # Use the integer label key for de_from_adata
    result = cplearn.de_from_adata(
        adata,
        key=use_key,  # Use integer label key
        layers_key=layers_key,  # Explicitly pass layers_key to use original cluster_key
        groups_a=groups_a,
        groups_b=groups_b,
        layer=layer,
        options=options,
    )
    
    # Clean up temporary integer label column if we created it
    if use_key == temp_key and temp_key in adata.obs:
        # Keep it for now in case it's needed, but could remove it
        pass
    
    # Ensure result has the expected column names and order
    # cplearn.de_from_adata returns columns: gene, statistic, pvalue, log2fc, effect_r, 
    # mean_a, mean_b, pct_expr_a, pct_expr_b, z_score, p_adj, ...
    # Reorder to match expected format: gene, log2fc, z_score, pvalue, p_adj, mean_a, mean_b, pct_expr_a, pct_expr_b
    expected_cols = ['gene', 'log2fc', 'z_score', 'pvalue', 'p_adj', 'mean_a', 'mean_b', 'pct_expr_a', 'pct_expr_b']
    available_cols = [col for col in expected_cols if col in result.columns]
    other_cols = [col for col in result.columns if col not in expected_cols]
    
    # Reorder columns: expected columns first, then others
    result = result[available_cols + other_cols]
    
    return result


# Alias for backward compatibility
run_differential_expression = marker_genes


def rank_genes_groups(
    adata: AnnData,
    groupby: str,
    *,
    groups: str | Sequence[str] = "all",
    reference: str = "rest",
    n_genes: int | None = None,
    method: str | None = None,
    use_raw: bool | None = None,
    layer: str | None = None,
    key_added: str | None = None,
    **kwargs,
) -> None:
    """
    Rank genes for characterizing groups (differential expression analysis).
    
    This function performs differential expression analysis using scanpy's rank_genes_groups,
    which supports multiple statistical methods (wilcoxon, t-test, logreg, etc.) and can
    compare multiple groups simultaneously.
    
    Parameters:
        adata: AnnData object
        groupby: Key in adata.obs that contains group labels
        groups: Subset of groups to compare. Default: "all"
        reference: Reference group for comparison. Default: "rest" (compare each group to all others)
        n_genes: Number of genes to rank. Default: None (all genes)
        method: Statistical method to use. Options: "wilcoxon" (default), "t-test", "logreg", etc.
        use_raw: Whether to use raw data. Default: None (auto-detect)
        layer: Layer to use for analysis. Default: None (use X)
        key_added: Key name for results in adata.uns. Default: None ("rank_genes_groups")
        **kwargs: Additional arguments passed to scanpy's rank_genes_groups
    
    Returns:
        None. Updates `adata.uns` with ranked genes results in `adata.uns[key_added or 'rank_genes_groups']`.
    """
    sc_tl.rank_genes_groups(
        adata,
        groupby=groupby,
        groups=groups,
        reference=reference,
        n_genes=n_genes,
        method=method,
        use_raw=use_raw,
        layer=layer,
        key_added=key_added,
        **kwargs,
    )


def filter_rank_genes_groups(
    adata: AnnData,
    *,
    key: str | None = None,
    groupby: str | None = None,
    use_raw: bool | None = None,
    key_added: str = "rank_genes_groups_filtered",
    min_in_group_fraction: float = 0.25,
    min_fold_change: float = 1,
    max_out_group_fraction: float = 0.5,
    compare_abs: bool = False,
) -> None:
    """
    Filter ranked genes groups based on expression criteria.
    
    This function filters the results from rank_genes_groups based on expression
    fraction, fold change, and other criteria to identify high-quality marker genes.
    
    Parameters:
        adata: AnnData object
        key: Key in adata.uns containing rank_genes_groups results. Default: None ("rank_genes_groups")
        groupby: Key in adata.obs that contains group labels. Default: None (auto-detect from key)
        use_raw: Whether to use raw data. Default: None (auto-detect)
        key_added: Key name for filtered results in adata.uns. Default: "rank_genes_groups_filtered"
        min_in_group_fraction: Minimum fraction of cells expressing gene in group. Default: 0.25
        min_fold_change: Minimum fold change. Default: 1
        max_out_group_fraction: Maximum fraction of cells expressing gene outside group. Default: 0.5
        compare_abs: Whether to compare absolute values. Default: False
    
    Returns:
        None. Updates `adata.uns` with filtered ranked genes results in `adata.uns[key_added]`.
    """
    sc_tl.filter_rank_genes_groups(
        adata,
        key=key,
        groupby=groupby,
        use_raw=use_raw,
        key_added=key_added,
        min_in_group_fraction=min_in_group_fraction,
        min_fold_change=min_fold_change,
        max_out_group_fraction=max_out_group_fraction,
        compare_abs=compare_abs,
    )


def marker_gene_overlap(
    adata: AnnData,
    reference_markers: dict[str, set[int] | list[int]],
    *,
    key: str = "rank_genes_groups",
    method: str = "overlap_count",
    normalize: str | None = None,
    top_n_markers: int | None = None,
    adj_pval_threshold: float | None = None,
    key_added: str = "marker_gene_overlap",
    inplace: bool = False,
) -> pd.DataFrame | None:
    """
    Analyze marker gene overlap between data and reference markers.
    
    This function compares marker genes identified in the data with reference marker
    genes (e.g., from literature or other datasets) to assess similarity and validate
    cell type annotations.
    
    Parameters:
        adata: AnnData object
        reference_markers: Dictionary mapping group names to sets/lists of marker gene names
        key: Key in adata.uns containing rank_genes_groups results. Default: "rank_genes_groups"
        method: Method for computing overlap. Options: "overlap_count", "jaccard", etc. Default: "overlap_count"
        normalize: Normalization method. Options: "reference", "data", or None. Default: None
        top_n_markers: Number of top markers to consider. Default: None (all)
        adj_pval_threshold: Adjusted p-value threshold for filtering. Default: None
        key_added: Key name for results in adata.uns. Default: "marker_gene_overlap"
        inplace: Whether to modify adata in place. Default: False
    
    Returns:
        pd.DataFrame | None: Overlap analysis results (if inplace=False)
    """
    return sc_tl.marker_gene_overlap(
        adata,
        reference_markers,
        key=key,
        method=method,
        normalize=normalize,
        top_n_markers=top_n_markers,
        adj_pval_threshold=adj_pval_threshold,
        key_added=key_added,
        inplace=inplace,
    )


__all__ = [
    "marker_genes",
    "run_differential_expression",
    "pick_groups",
    "rank_genes_groups",
    "filter_rank_genes_groups",
    "marker_gene_overlap",
]