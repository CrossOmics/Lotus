from importlib import import_module
from typing import Any

__all__ = [
    "aggregate",
    "annotate_cell_types",
    "filter_rank_genes_groups",
    "ingest",
    "marker_gene_overlap",
    "obs_df",
    "rank_genes_groups",
    "rank_genes_groups_df",
    "run_celltypist_annotation",
    "run_enrichr_analysis",
    "run_enrichment",
    "score_genes",
    "score_genes_cell_cycle",
    "var_df",
]

_MODULE_EXPORTS = {
    "deg_analysis_core": {
        "rank_genes_groups",
        "filter_rank_genes_groups",
        "marker_gene_overlap",
        "score_genes",
        "score_genes_cell_cycle",
        "rank_genes_groups_df",
    },
    "_get_core": {
        "aggregate",
        "obs_df",
        "var_df",
    },
    "celltypist_core": {
        "run_celltypist_annotation",
    },
    "gseapy_core": {
        "run_enrichr_analysis",
    },
    "_annotation_workflows": {
        "annotate_cell_types",
        "run_enrichment",
    },
    "ingest_core": {
        "ingest",
    },
}


def __getattr__(name: str) -> Any:
    for module_name, exports in _MODULE_EXPORTS.items():
        if name in exports:
            module = import_module(f"lotus.annotation_analysis.{module_name}")
            value = getattr(module, name)
            globals()[name] = value
            return value
    raise AttributeError(f"module 'lotus.annotation_analysis' has no attribute {name!r}")
