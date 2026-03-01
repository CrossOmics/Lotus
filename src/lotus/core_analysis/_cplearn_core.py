from typing import Tuple, Optional
import numpy as np
import pandas as pd
import json
from anndata import AnnData
from lotus.lineagetracker import logged
from .cplearn.corespect import CorespectModel
from .cplearn.corespect.config import CoreSpectConfig


def _corespect_impl(
        adata: AnnData,
        *,
        use_rep: str = "X_pca",
        key_added: str = "cplearn",
        q: int = 20,
        r: int = 10,
        core_frac: float = 0.2,
        densify: bool = False,
        granularity: float = 0.5,
        resolution: float = 0.5,
        fine_grained: bool = True,
        propagate: bool = True,
        copy: bool = False,
) -> Optional[Tuple[AnnData, CorespectModel]]:
    """Run CoreSpect and write cluster labels, is_core, layer to adata.obs/uns. Internal only."""

    # Handle In-Place vs Copy Mode
    if copy:
        adata = adata.copy()

    # Prepare Data
    if use_rep not in adata.obsm:
        raise KeyError(f"Could not find representation '{use_rep}' in adata.obsm")

    X = adata.obsm[use_rep]

    # Configure Model
    cfg = CoreSpectConfig(
        q=q,
        r=r,
        core_frac=core_frac,
        densify=densify,
        granularity=granularity,
        resolution=resolution
    ).configure()

    # Initialize and Run Model
    model = CorespectModel(X, **cfg.unpack())

    # Run the core logic
    model.run(fine_grained=fine_grained, propagate=propagate)

    # Write Results back to AnnData

    # Write Cluster Labels (CoreSpect classification: one label per cell)
    labels = model.labels_.astype(str)
    labels[labels == "-1"] = "Unassigned"
    adata.obs[key_added] = pd.Categorical(labels)
    adata.obs[f"{key_added}_cluster"] = adata.obs[key_added]

    # Write Core Identity (Boolean mask) and per-cell layer index (0=core, 1, 2, ...)
    is_core = np.zeros(adata.n_obs, dtype=bool)
    layer_index = np.full(adata.n_obs, -1, dtype=np.int32)  # -1 = not in any layer
    if model.layers_ and len(model.layers_) > 0:
        core_indices = model.layers_[0]
        is_core[core_indices] = True
        for k, layer_cells in enumerate(model.layers_):
            idx = np.asarray(layer_cells).ravel()
            layer_index[idx] = k

    core_key = f"{key_added}_is_core"
    adata.obs[core_key] = is_core
    layer_key = f"{key_added}_layer"
    adata.obs[layer_key] = pd.Series(layer_index, index=adata.obs.index, dtype="int32")

    # Prepare Ragged Array for HDF5/h5py Compatibility
    # model.layers_ is a list of lists with different lengths (Ragged Array).
    # HDF5 cannot save this directly. We convert it to a JSON string.

    # Convert numpy arrays to standard python lists (for JSON serialization)
    layers_serializable = [l.tolist() if hasattr(l, "tolist") else list(l) for l in model.layers_]

    # Dump to JSON string
    layers_json_str = json.dumps(layers_serializable)

    # Write Parameters and Layers to uns (Metadata)
    adata.uns[key_added] = {
        'params': {
            'q': q, 'r': r, 'core_frac': core_frac,
            'densify': densify, 'resolution': resolution,
            'fine_grained': fine_grained, 'propagate': propagate
        },
        # Store as JSON string.
        # Downstream tools should use json.loads(adata.uns[key]['layers_indices_json'])
        'layers_indices_json': layers_json_str,
        # So we can skip CoreSpect when user only changes percentage (and did not filter cells)
        'n_obs_at_corespect': int(adata.n_obs),
    }

    # Return Logic
    return adata, model


def _run_corespect(
        adata: AnnData,
        *,
        key_added: str = "cplearn",
        force_recalc: bool = False,
        copy: bool = False,
        **kwargs,
) -> tuple[AnnData, CorespectModel | None]:
    """Run CoreSpect with cache check; used only by core_selection. Internal only."""
    if copy:
        adata = adata.copy()

    layer_key = f"{key_added}_layer"
    n_obs = adata.n_obs
    # Reuse only if layer column exists, length matches current cell count, and was computed on this n_obs (no filter since)
    uns_c = adata.uns.get(key_added) or {}
    stored_n = uns_c.get("n_obs_at_corespect")
    has_layer = layer_key in adata.obs and len(adata.obs[layer_key]) == n_obs
    can_skip = (
        not force_recalc
        and has_layer
        and stored_n is not None
        and int(stored_n) == n_obs
    )
    if can_skip:
        print(
            f"Info: CoreSpect layer found for this cell set (`adata.obs['{layer_key}']`, n_obs={n_obs}). "
            "Skipping recalculation. Use `force_recalc=True` to override."
        )
        return adata, None

    if force_recalc:
        print("Info: Force recalculation enabled. Re-running CoreSpect...")

    adata, model = _corespect_impl(
        adata,
        key_added=key_added,
        copy=True,
        **kwargs,
    )
    return adata, model

@logged
def core_selection(
        adata: AnnData,
        percentage: float,
        *,
        key_added: str = "cplearn",
        force_recalc: bool = False,
        copy: bool = False,
        **kwargs,
) -> Tuple[AnnData, Optional[CorespectModel]]:
    """
    Run cplearn (CoreSpect) and select a fraction of cells as "core" by layer order.

    This is the main user-facing API: you choose a percentage (e.g. 0.3 = 30%),
    and the function runs CoreSpect, then labels the top percentage of cells
    (when sorted by layer: layer 0 first, then 1, 2, ...) as "core", the rest as "noncore".

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix (must have representation in adata.obsm, e.g. "X_pca").
    percentage : float
        Fraction of cells to label as "core", in (0, 1]. E.g. 0.3 = 30% core, 70% noncore.
    key_added : str
        Prefix for new columns in adata.obs and adata.uns.
    force_recalc : bool
        If True, re-run CoreSpect even when results exist.
    copy : bool
        If True, copy adata before modifying.
    **kwargs
        Passed to CoreSpect (e.g. use_rep, q, r, resolution).

    Returns
    -------
    adata : AnnData
        Same object (or copy) with new columns:
        - key_added + "_core_selection": "core" | "noncore" (categorical).
        - key_added, key_added + "_cluster": CoreSpect cluster label per cell (e.g. "0", "1", "Unassigned").
        - key_added + "_is_core", key_added + "_layer" (from CoreSpect).
    model : CorespectModel or None
        The fitted model (None if results were loaded from cache). Use for coremap().

    Notes
    -----
    Core cells:  adata[adata.obs[key_added + "_core_selection"] == "core"]
    Non-core:    adata[adata.obs[key_added + "_core_selection"] == "noncore"]
    """
    if not (0.0 < percentage <= 1.0):
        raise ValueError("percentage must be in (0, 1].")

    if copy:
        adata = adata.copy()

    adata, model = _run_corespect(
        adata,
        key_added=key_added,
        force_recalc=force_recalc,
        copy=False,
        **kwargs,
    )

    layer_key = f"{key_added}_layer"
    if layer_key not in adata.obs:
        raise RuntimeError(
            f"Expected '{layer_key}' in adata.obs. "
            "Re-run with force_recalc=True or check key_added."
        )

    n_obs = adata.n_obs
    target_n = int(np.ceil(percentage * n_obs))
    target_n = min(target_n, n_obs)

    # Take whole layers only: 0, then 1, then 2, ... then -1; stop when cumulative >= target_n (never split a layer)
    layer_values = adata.obs[layer_key].values
    # Layer order: 0, 1, 2, ... then -1 (unassigned last)
    unique_layers = sorted(
        np.unique(layer_values).tolist(),
        key=lambda x: (x == -1, x),
    )
    layers_in_core = []
    cum = 0
    for L in unique_layers:
        cnt = int((layer_values == L).sum())
        layers_in_core.append(L)
        cum += cnt
        if cum >= target_n:
            break
    core_mask = np.isin(layer_values, layers_in_core)
    n_core = int(core_mask.sum())

    selection = np.where(core_mask, "core", "noncore")
    selection_key = f"{key_added}_core_selection"
    adata.obs[selection_key] = pd.Categorical(selection, categories=["core", "noncore"])

    actual_frac = n_core / n_obs
    if key_added not in adata.uns:
        adata.uns[key_added] = {}
    adata.uns[key_added]["core_selection"] = {
        "percentage_requested": percentage,
        "actual_fraction_core": actual_frac,
        "n_core": n_core,
        "n_noncore": int(n_obs - n_core),
        "layers_in_core": [int(x) for x in layers_in_core],
    }

    return adata, model
