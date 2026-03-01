from collections.abc import Iterable, Sequence
from typing import Any

import anndata
import celltypist
import pandas as pd
from celltypist import models
from loguru import logger

from anndata import AnnData
from lotus.lineagetracker import logged
from ._storage import (
    ANNOTATION_SCHEMA_VERSION,
    build_obs_column_name,
    sanitize_storage_key,
)


_AVAILABLE_MODELS: list[str] | None = None


def _normalize_model_names(model_names: str | Sequence[str] | None) -> list[str]:
    if model_names is None:
        return [models.get_default_model()]
    if isinstance(model_names, str):
        return [model_names]
    normalized = [name for name in model_names]
    if not normalized:
        raise ValueError("`model_names` must contain at least one model.")
    return normalized


def _sanitize_model_key(model_name: str) -> str:
    return sanitize_storage_key(model_name)


def _list_celltypist_models() -> list[str]:
    global _AVAILABLE_MODELS
    if _AVAILABLE_MODELS is not None:
        return _AVAILABLE_MODELS
    try:
        _AVAILABLE_MODELS = models.get_all_models()
    except Exception as exc:
        logger.warning("Unable to enumerate CellTypist models: {}", exc)
        _AVAILABLE_MODELS = []
    return _AVAILABLE_MODELS


def _resolve_model_names(model_names: Sequence[str]) -> list[str]:
    available_models = set(_list_celltypist_models())
    missing = [
        model_name
        for model_name in model_names
        if "/" not in model_name and model_name not in available_models
    ]
    if missing:
        raise ValueError(
            "CellTypist model(s) not available locally: "
            f"{missing}. Call `celltypist.models.models_description()` or "
            "`download_models_to_cache()` first."
        )
    return list(model_names)


def _summarize_clusters(
        adata: AnnData,
        predicted_obs: pd.DataFrame,
        prediction_column: str,
        cluster_column: str,
) -> dict[str, Any]:
    cluster_summary: dict[str, Any] = {}
    combined_obs = pd.DataFrame(
        {
            "cluster": adata.obs[cluster_column].values,
            "prediction": predicted_obs[prediction_column].values,
            "conf_score": predicted_obs["conf_score"].values if "conf_score" in predicted_obs else 0.0,
        },
        index=adata.obs.index,
    )
    for cluster_id, group in combined_obs.groupby("cluster", observed=True):
        if group.empty:
            continue
        top_label = group["prediction"].mode()[0]
        top_count = int((group["prediction"] == top_label).sum())
        total_cells = int(len(group))
        agreement_fraction = top_count / total_cells if total_cells else 0.0
        cluster_summary[str(cluster_id)] = {
            "predicted_cell_type": top_label,
            "cell_count": total_cells,
            "agreement_fraction": float(f"{agreement_fraction:.4f}"),
            "average_confidence": float(f"{group['conf_score'].mean():.4f}"),
        }
    return cluster_summary


@logged
def run_celltypist_annotation(
        adata: anndata.AnnData,
        model_names: str | Sequence[str] | None = None,
        *,
        majority_voting: bool = True,
        target_cluster_col: str | None = None,
        mode: str = "best_match",
        p_thres: float = 0.5,
        key_added: str = "celltypist",
        store_probabilities: bool = False,
        store_decision_scores: bool = False,
        inplace: bool = True,
) -> AnnData | dict[str, Any]:
    """
    Run CellTypist annotation and store results in AnnData.

    Per-cell labels and confidence scores are written to `adata.obs`.
    Run-level metadata and summaries are written to `adata.uns[key_added]`.
    """
    if not isinstance(adata, AnnData):
        raise TypeError("`adata` must be an AnnData object.")
    if target_cluster_col is not None and target_cluster_col not in adata.obs:
        raise KeyError(f"`{target_cluster_col}` not found in `adata.obs`.")

    normalized_mode = mode.replace("_", " ")
    if normalized_mode not in {"best match", "prob match"}:
        raise ValueError("`mode` must be one of {'best_match', 'prob_match'}.")

    model_name_list = _resolve_model_names(_normalize_model_names(model_names))
    result_adata = adata.copy() if not inplace else adata
    uns_bucket: dict[str, Any] = result_adata.uns.setdefault(key_added, {})
    uns_bucket["schema_version"] = ANNOTATION_SCHEMA_VERSION
    uns_bucket["provider"] = "celltypist"
    uns_bucket["analysis_type"] = "cell_type_annotation"
    model_results: dict[str, Any] = {}

    for model_name in model_name_list:
        model_key = _sanitize_model_key(model_name)
        loaded_model = models.Model.load(model=model_name)
        prediction = celltypist.annotate(
            result_adata,
            model=loaded_model,
            majority_voting=majority_voting,
            mode=normalized_mode,
            p_thres=p_thres,
        )
        has_majority_voting = "majority_voting" in prediction.predicted_labels.columns
        confidence_source = "majority_voting" if majority_voting and has_majority_voting else "predicted_labels"

        prediction_obs = prediction.to_adata(
            insert_labels=True,
            insert_conf=True,
            insert_conf_by=confidence_source,
            insert_decision=store_decision_scores,
            insert_prob=store_probabilities,
            prefix=f"{key_added}_{model_key}_",
        ).obs

        predicted_label_col = build_obs_column_name(key_added, model_key, "predicted_labels")
        majority_voting_col = build_obs_column_name(key_added, model_key, "majority_voting")
        conf_score_col = build_obs_column_name(key_added, model_key, "conf_score")

        result_adata.obs[predicted_label_col] = prediction_obs[f"{key_added}_{model_key}_predicted_labels"]
        if has_majority_voting and f"{key_added}_{model_key}_majority_voting" in prediction_obs:
            result_adata.obs[majority_voting_col] = prediction_obs[majority_voting_col]
        if conf_score_col in prediction_obs:
            result_adata.obs[conf_score_col] = prediction_obs[conf_score_col]

        cluster_summary = None
        if target_cluster_col is not None:
            summary_input = pd.DataFrame(index=result_adata.obs.index)
            summary_prediction_column = predicted_label_col
            summary_input[predicted_label_col] = result_adata.obs[predicted_label_col]
            if majority_voting and majority_voting_col in result_adata.obs:
                summary_input[majority_voting_col] = result_adata.obs[majority_voting_col]
                summary_prediction_column = majority_voting_col
            if conf_score_col in result_adata.obs:
                summary_input["conf_score"] = result_adata.obs[conf_score_col]
            cluster_summary = _summarize_clusters(
                result_adata,
                predicted_obs=summary_input,
                prediction_column=summary_prediction_column,
                cluster_column=target_cluster_col,
            )

        model_results[model_key] = {
            "model_name": model_name,
            "model_key": model_key,
            "cell_count": prediction.cell_count,
            "parameters": {
                "majority_voting_requested": majority_voting,
                "majority_voting_applied": has_majority_voting,
                "target_cluster_col": target_cluster_col,
                "mode": normalized_mode,
                "p_thres": p_thres,
                "store_probabilities": store_probabilities,
                "store_decision_scores": store_decision_scores,
                "confidence_source": confidence_source,
            },
            "provider": "celltypist",
            "model_description": dict(loaded_model.description),
            "obs_columns": {
                "predicted_labels": predicted_label_col,
                "majority_voting": majority_voting_col if has_majority_voting and majority_voting_col in result_adata.obs else None,
                "conf_score": conf_score_col if conf_score_col in result_adata.obs else None,
            },
            "summary_frequency_predicted_labels": prediction.summary_frequency("predicted_labels").to_dict("records"),
            "summary_frequency_majority_voting": (
                prediction.summary_frequency("majority_voting").to_dict("records")
                if has_majority_voting
                else None
            ),
            "cluster_summary": cluster_summary,
            "available_cell_types": loaded_model.cell_types.tolist(),
            "feature_count": int(len(loaded_model.features)),
        }

        if store_probabilities:
            model_results[model_key]["probability_columns"] = [
                f"{key_added}_{model_key}_{column}"
                for column in prediction.probability_matrix.columns
            ]
        if store_decision_scores and not store_probabilities:
            model_results[model_key]["decision_columns"] = [
                f"{key_added}_{model_key}_{column}"
                for column in prediction.decision_matrix.columns
            ]

    uns_bucket["models"] = model_results
    uns_bucket["model_order"] = list(model_results)
    return result_adata if not inplace else uns_bucket


def download_models_to_cache(model_names: str | Iterable[str] | None = None) -> None:
    """
    Pre-download CellTypist model(s) to the local cache.
    """
    if model_names is None:
        models.download_models()
        return
    if isinstance(model_names, str):
        model_names = [model_names]
    for model_name in model_names:
        models.download_model(model_name)
