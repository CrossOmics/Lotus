from typing import Tuple, Optional
import numpy as np
import pandas as pd
import json
from anndata import AnnData
from .cplearn.corespect import CorespectModel
from .cplearn.corespect.config import CoreSpectConfig


def corespect(
        adata: AnnData,
        *,
        use_rep: str = "X_pca",
        key_added: str = "cplearn",
        # CoreSpect Config Parameters
        q: int = 20,
        r: int = 10,
        core_frac: float = 0.2,
        densify: bool = False,
        granularity: float = 0.5,
        resolution: float = 0.5,
        # Run Parameters
        fine_grained: bool = True,
        propagate: bool = True,
        copy: bool = False,
) -> Optional[Tuple[AnnData, CorespectModel]]:
    """
    Cluster cells using the CoreSpect algorithm.
    """

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

    # Write Cluster Labels
    labels = model.labels_.astype(str)
    labels[labels == "-1"] = "Unassigned"
    adata.obs[key_added] = pd.Categorical(labels)

    # Write Core Identity (Boolean mask)
    is_core = np.zeros(adata.n_obs, dtype=bool)
    if model.layers_ and len(model.layers_) > 0:
        core_indices = model.layers_[0]
        is_core[core_indices] = True

    core_key = f"{key_added}_is_core"
    adata.obs[core_key] = is_core

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
        'layers_indices_json': layers_json_str
    }

    # Return Logic
    return adata, model

