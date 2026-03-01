from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import gseapy as gp
import pandas as pd
from anndata import AnnData
from loguru import logger

from lotus.lineagetracker import logged
from ._storage import ANNOTATION_SCHEMA_VERSION, _sanitize_storage_key


_AVAILABLE_LIBRARIES: dict[str, list[str]] = {}


def _normalize_gene_sets(gene_sets: str | Sequence[str] | Mapping[str, Sequence[str]]) -> list[str] | Mapping[str, Sequence[str]]:
    if isinstance(gene_sets, Mapping):
        return gene_sets
    if isinstance(gene_sets, str):
        return [gene_sets]
    normalized = [gene_set for gene_set in gene_sets]
    if not normalized:
        raise ValueError("`gene_sets` must contain at least one library or custom gene set.")
    return normalized


def _extract_cell_type_from_term(term_name: str, p_value: float) -> str | None:
    if not term_name or pd.isna(term_name):
        return None
    if pd.isna(p_value) or float(p_value) > 0.05:
        return None
    term = str(term_name).strip()
    if ":" in term:
        return term.split(":", 1)[0].strip()
    return None


def _list_enrichr_libraries(organism: str) -> list[str]:
    organism_key = organism.lower()
    if organism_key in _AVAILABLE_LIBRARIES:
        return _AVAILABLE_LIBRARIES[organism_key]
    try:
        _AVAILABLE_LIBRARIES[organism_key] = gp.get_library_name(organism=organism)
    except Exception as exc:
        logger.warning("Unable to enumerate Enrichr libraries for {}: {}", organism, exc)
        _AVAILABLE_LIBRARIES[organism_key] = []
    return _AVAILABLE_LIBRARIES[organism_key]


def _resolve_gene_sets(
        gene_sets: list[str] | Mapping[str, Sequence[str]],
        organism: str,
) -> tuple[list[str], list[str], Mapping[str, Sequence[str]] | None]:
    if isinstance(gene_sets, Mapping):
        return [], [], gene_sets

    available_libraries = set(_list_enrichr_libraries(organism))
    missing_libraries: list[str] = []
    enrichr_libraries: list[str] = []
    for gene_set in gene_sets:
        if Path(gene_set).exists():
            enrichr_libraries.append(gene_set)
            continue
        if available_libraries and gene_set not in available_libraries:
            missing_libraries.append(gene_set)
            continue
        enrichr_libraries.append(gene_set)
    return enrichr_libraries, missing_libraries, None


@logged
def run_enrichr_analysis(
        adata: AnnData,
        gene_list: Sequence[str],
        gene_sets: str | Sequence[str] | Mapping[str, Sequence[str]],
        *,
        organism: str = "human",
        limit: int | None = None,
        pval_cutoff: float = 0.05,
        key_added: str = "gseapy",
        analysis_label: str | None = None,
        background: list[str] | int | str | None = None,
        inplace: bool = True,
) -> AnnData | dict[str, Any]:
    """
    Run Enrichr analysis and store enrichment results in `adata.uns`.

    GSEApy/Enrichr does not produce per-cell outputs, so all results are stored in `uns`.
    """
    if not isinstance(adata, AnnData):
        raise TypeError("`adata` must be an AnnData object.")
    normalized_gene_list = [gene for gene in gene_list if gene]
    if not normalized_gene_list:
        raise ValueError("`gene_list` must contain at least one non-empty gene symbol.")

    normalized_gene_sets = _normalize_gene_sets(gene_sets)
    resolved_gene_sets, missing_libraries, custom_gene_sets = _resolve_gene_sets(
        normalized_gene_sets,
        organism=organism,
    )
    if missing_libraries:
        raise ValueError(
            "Enrichr libraries not available for the selected organism: "
            f"{missing_libraries}. Call `gseapy.get_library_name()` to inspect valid libraries."
        )

    result_adata = adata.copy() if not inplace else adata
    uns_bucket: dict[str, Any] = result_adata.uns.setdefault(key_added, {})
    metadata_keys = {"schema_version", "provider", "analysis_type"}
    uns_bucket["schema_version"] = ANNOTATION_SCHEMA_VERSION
    uns_bucket["provider"] = "gseapy"
    uns_bucket["analysis_type"] = "enrichment"
    result_key = _sanitize_storage_key(analysis_label) if analysis_label else f"analysis_{len([k for k in uns_bucket if k not in metadata_keys]) + 1}"
    libraries_to_run: list[str] | Mapping[str, Sequence[str]]
    libraries_to_run = custom_gene_sets if custom_gene_sets is not None else resolved_gene_sets

    enrichment = gp.enrichr(
        gene_list=normalized_gene_list,
        gene_sets=libraries_to_run,
        organism=organism,
        background=background,
        outdir=None,
        no_plot=True,
        verbose=False,
    )

    results_by_library: dict[str, Any] = {}
    if custom_gene_sets is not None:
        library_names = list(custom_gene_sets)
    else:
        library_names = resolved_gene_sets

    result_frame = enrichment.results.copy() if enrichment.results is not None else pd.DataFrame()
    for library_name in library_names:
        if "Gene_set" in result_frame.columns:
            library_frame = result_frame[result_frame["Gene_set"] == library_name].copy()
        else:
            library_frame = result_frame.copy()
        if (
            custom_gene_sets is not None
            and (library_frame is None or library_frame.empty)
            and "Term" in result_frame.columns
        ):
            library_frame = result_frame[result_frame["Term"] == library_name].copy()
        if library_frame is None or library_frame.empty:
            results_by_library[library_name] = {
                "top_terms": [],
                "term_count": 0,
                "inferred_cell_types": [],
                "ranking_metric": None,
            }
            continue

        filtered_frame = library_frame.copy()
        if "Adjusted P-value" in filtered_frame.columns:
            filtered_frame = filtered_frame[filtered_frame["Adjusted P-value"] < pval_cutoff]
        elif "P-value" in filtered_frame.columns:
            filtered_frame = filtered_frame[filtered_frame["P-value"] < pval_cutoff]

        top_frame = filtered_frame.head(limit) if limit is not None else filtered_frame
        inferred_cell_types: list[str] = []
        if "CellMarker" in library_name or "cellmarker" in library_name.lower():
            pval_column = "Adjusted P-value" if "Adjusted P-value" in top_frame.columns else "P-value"
            term_column = "Term" if "Term" in top_frame.columns else top_frame.columns[0]
            inferred_cell_types = sorted(
                {
                    cell_type
                    for _, row in top_frame.iterrows()
                    if (cell_type := _extract_cell_type_from_term(row.get(term_column), row.get(pval_column)))
                }
            )

        results_by_library[library_name] = {
            "top_terms": top_frame.to_dict("records"),
            "term_count": int(len(top_frame)),
            "inferred_cell_types": inferred_cell_types,
            "ranking_metric": (
                "Combined Score"
                if "Combined Score" in top_frame.columns
                else "Odds Ratio" if "Odds Ratio" in top_frame.columns else None
            ),
        }

    uns_bucket[result_key] = {
        "provider": "gseapy",
        "analysis_label": result_key,
        "gene_list": normalized_gene_list,
        "gene_list_size": len(normalized_gene_list),
        "gene_sets": library_names,
        "organism": organism,
        "parameters": {
            "limit": limit,
            "pval_cutoff": pval_cutoff,
            "background": background,
        },
        "results": results_by_library,
    }
    return result_adata if not inplace else uns_bucket[result_key]
