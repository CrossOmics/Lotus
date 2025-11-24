from __future__ import annotations

from dataclasses import asdict
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from ... import cplearn as _cplearn

__all__ = [
    "corespect",
    "coremap_embedding",
    "anchored_map",
    "Coremap",
    "find_anchors",
    "cluster_core",
    "stable_core",
    "fine_grained_core",
    "propagate_from_core",
    "FlowRank",
    "FlowRankConfig",
    "StableCoreConfig",
    "FineGrainedConfig",
    "ClusterConfig",
    "PropagationConfig",
    "CorespectModel",
    "DEAnalyzer",
    "DEOptions",
    "de_from_model",
    "de_from_adata",
]

_WRAPPED_EXPORTS = {
    "corespect",
    "coremap_embedding",
    "de_from_adata",
}


def _get_array_from_adata(
    adata: AnnData,
    *,
    layer: str | None = None,
    use_rep: str | None = None,
) -> np.ndarray:
    if layer is not None and use_rep is not None:
        msg = "Specify at most one of `layer` and `use_rep`."
        raise ValueError(msg)
    if layer is not None:
        if layer not in adata.layers:
            raise KeyError(f"Layer {layer!r} not found in `adata.layers`.")
        matrix = adata.layers[layer]
    elif use_rep is not None:
        if use_rep not in adata.obsm:
            raise KeyError(f"Representation {use_rep!r} not found in `adata.obsm`.")
        matrix = adata.obsm[use_rep]
    else:
        matrix = adata.X

    if isinstance(matrix, pd.DataFrame):
        matrix = matrix.to_numpy()
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    return np.asarray(matrix, dtype=float, order="C")


def corespect(
    adata: AnnData,
    *,
    layer: str | None = None,
    use_rep: str | None = None,
    flowrank: Mapping[str, Any] | None = None,
    stable: Mapping[str, Any] | None = None,
    fine: Mapping[str, Any] | None = None,
    cluster: Mapping[str, Any] | None = None,
    propagate_cfg: Mapping[str, Any] | None = None,
    fine_grained: bool = False,
    propagate: bool = True,
    key_added: str = "cplearn",
    copy: bool = False,
) -> _cplearn.CorespectModel | tuple[AnnData, _cplearn.CorespectModel]:
    """
    Run the CPlearn CoreSpect pipeline on an AnnData object.

    Parameters
    ----------
    adata
        Input AnnData object.
    layer, use_rep
        Choose which data matrix to use. Defaults to `adata.X`.
    flowrank, stable, fine, cluster, propagate_cfg
        Optional dictionaries overriding the respective configuration dataclasses.
    fine_grained
        Whether to run the fine-grained refinement stage.
    propagate
        Whether to run diffusion / propagation after clustering.
    key_added
        Column name in `adata.obs` receiving cluster labels.
    copy
        If `True`, return a copied AnnData alongside the fitted model.

    Returns
    -------
    If ``copy=False`` (default) returns the fitted :class:`~cplearn.corespect.corespect.CorespectModel`.
    If ``copy=True`` returns ``(adata_copy, model)``.
    """

    X = _get_array_from_adata(adata, layer=layer, use_rep=use_rep)

    FlowRankConfig = getattr(_cplearn, "FlowRankConfig")
    StableCoreConfig = getattr(_cplearn, "StableCoreConfig")
    FineGrainedConfig = getattr(_cplearn, "FineGrainedConfig")
    ClusterConfig = getattr(_cplearn, "ClusterConfig")
    PropagationConfig = getattr(_cplearn, "PropagationConfig")
    CorespectModel = getattr(_cplearn, "CorespectModel")

    flowrank_cfg = FlowRankConfig(**(flowrank or {}))
    stable_cfg = StableCoreConfig(**(stable or {}))
    fine_cfg = FineGrainedConfig(**(fine or {}))
    cluster_cfg = ClusterConfig(**(cluster or {}))
    prop_cfg = PropagationConfig(**(propagate_cfg or {}))

    model = CorespectModel(
        X,
        flowrank_cfg=flowrank_cfg,
        stable_cfg=stable_cfg,
        fine_core_cfg=fine_cfg,
        cluster_cfg=cluster_cfg,
        prop_cfg=prop_cfg,
    )
    model.run(fine_grained=fine_grained, propagate=propagate)

    target = adata.copy() if copy else adata

    labels = np.asarray(model.labels_, dtype=int)
    target.obs[key_added] = pd.Categorical(labels)

    cplearn_info = {
        "layers": [np.array(layer).astype(int).tolist() for layer in (model.layers_ or [])],
        "flowrank_score": dict(model.flowrank_score_) if model.flowrank_score_ is not None else None,
        "config": {
            "flowrank": asdict(flowrank_cfg),
            "stable": asdict(stable_cfg),
            "fine": asdict(fine_cfg),
            "cluster": asdict(cluster_cfg),
            "propagation": asdict(prop_cfg),
            "fine_grained": fine_grained,
            "propagate": propagate,
        },
    }
    
    # Save core_layers_ information if available
    if hasattr(model, "core_layers_") and model.core_layers_ is not None:
        core_layers = model.core_layers_
        cplearn_info["core_layers"] = [
            np.array(layer).astype(int).tolist() 
            for layer in core_layers
        ]
        # Flatten all core layer indices
        core_indices = set()
        for core_layer in core_layers:
            core_layer_arr = np.array(core_layer).flatten()
            core_indices.update(core_layer_arr.astype(int))
        cplearn_info["core_indices"] = sorted(list(core_indices))
        cplearn_info["n_core_cells"] = len(core_indices)
        
        # Create is_core marker in obs for easy visualization
        is_core = np.zeros(target.n_obs, dtype=bool)
        is_core[list(core_indices)] = True
        target.obs[f"{key_added}_is_core"] = pd.Categorical(is_core)
    
    target.uns[f"{key_added}_cplearn"] = cplearn_info

    if copy:
        return target, model
    return model


def coremap_embedding(
    adata: AnnData,
    *,
    model: _cplearn.CorespectModel,
    use_rep: str | None = "X_umap",
    layer: str | None = None,
    key_added: str = "X_cplearn_coremap",
    fast_view: bool = True,
    anchor_finding_mode: str = "default",
    anchor_reach: float | None = None,
    copy: bool = False,
) -> AnnData | None:
    """
    Compute anchored map embeddings from a fitted :class:`CorespectModel`.

    Parameters
    ----------
    adata
        AnnData containing the data used to fit the model.
    model
        Fitted :class:`cplearn.corespect.corespect.CorespectModel` instance returned by :func:`corespect`.
    use_rep, layer
        Optional existing embedding to initialise the anchored map. Defaults to ``adata.obsm['X_umap']`` if present.
    key_added
        Key in ``adata.obsm`` to store the resulting embedding.
    fast_view, anchor_finding_mode, anchor_reach
        Arguments forwarded to :class:`cplearn.coremap.coremap.Coremap`.
    copy
        If ``True`` return a copy of ``adata`` with the embedding stored.
    """

    Coremap = getattr(_cplearn, "Coremap")
    anchored_map = getattr(_cplearn, "anchored_map")

    if use_rep is not None and use_rep in adata.obsm:
        global_umap = np.asarray(adata.obsm[use_rep])
    elif layer is not None:
        global_umap = _get_array_from_adata(adata, layer=layer)
    else:
        global_umap = None

    coremap_obj = Coremap(
        model,
        global_umap=global_umap,
        fast_view=fast_view,
        anchor_finding_mode=anchor_finding_mode,
        anchor_reach=anchor_reach,
    )
    embedding_history = anchored_map(coremap_obj)
    if not embedding_history:
        raise RuntimeError("Anchored map did not produce any embeddings.")

    layers = getattr(model, "layers_", None) or []
    if not layers:
        raise RuntimeError(
            "CorespectModel has no layers information; cannot align coremap embeddings."
        )

    ordered_indices: list[int] = [
        int(index) for layer in layers for index in layer  # preserve iteration order
    ]
    if not ordered_indices:
        raise RuntimeError(
            "CorespectModel layers are empty; anchored map cannot be aligned to AnnData."
        )

    final_round = max(embedding_history, key=lambda k: k)
    coordinates = np.asarray(embedding_history[final_round], dtype=float)
    if coordinates.shape[0] != len(ordered_indices):
        raise RuntimeError(
            "Anchored map output does not match the number of indexed observations."
        )

    full_coordinates = np.full(
        (adata.n_obs, coordinates.shape[1]),
        np.nan,
        dtype=float,
    )
    full_coordinates[np.asarray(ordered_indices, dtype=int)] = coordinates

    target = adata.copy() if copy else adata
    target.obsm[key_added] = full_coordinates
    history_payload: dict[str, dict[str, list[float] | list[int]]] = {}
    for round_id, coords in embedding_history.items():
        round_indices = [
            int(index)
            for layer in layers[: round_id + 1]
            for index in layer
        ]
        history_payload[str(round_id)] = {
            "indices": round_indices,
            "coordinates": np.asarray(coords, dtype=float).tolist(),
        }

    target.uns[f"{key_added}_history"] = history_payload
    target.uns[f"{key_added}_ordered_indices"] = ordered_indices

    return target if copy else None


def __getattr__(name: str) -> Any:
    if name in _WRAPPED_EXPORTS and name in globals():
        return globals()[name]
    return getattr(_cplearn, name)


def __dir__() -> list[str]:
    return sorted(__all__)
