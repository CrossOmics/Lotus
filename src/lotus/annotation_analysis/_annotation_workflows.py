from collections.abc import Mapping, Sequence
from typing import Any

from anndata import AnnData

from .celltypist_core import run_celltypist_annotation
from .gseapy_core import run_enrichr_analysis


def annotate_cell_types(
        adata: AnnData,
        model_names: str | Sequence[str] | None = None,
        **kwargs,
) -> AnnData | dict[str, Any]:
    """
    Unified entrypoint for cell type annotation workflows.

    This currently delegates to the CellTypist-backed implementation.
    """
    return run_celltypist_annotation(
        adata,
        model_names=model_names,
        **kwargs,
    )


def run_enrichment(
        adata: AnnData,
        gene_list: Sequence[str],
        gene_sets: str | Sequence[str] | Mapping[str, Sequence[str]],
        **kwargs,
) -> AnnData | dict[str, Any]:
    """
    Unified entrypoint for enrichment workflows.

    This currently delegates to the Enrichr-backed GSEApy implementation.
    """
    return run_enrichr_analysis(
        adata,
        gene_list=gene_list,
        gene_sets=gene_sets,
        **kwargs,
    )
