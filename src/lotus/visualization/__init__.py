from importlib import import_module
from typing import Any

__all__ = [
    "coremap",
    "scatter",
    "heatmap",
    "dotplot",
    "tracksplot",
    "violin",
    "stacked_violin",
    "matrixplot",
    "clustermap",
    "ranking",
    "dendrogram",
    "sim",
    "tsne",
    "umap",
    "diffmap",
    "draw_graph",
    "spatial",
    "embedding",
    "embedding_density",
    "pca",
    "pca_loadings",
    "pca_variance_ratio",
    "pca_overview",
    "scrublet_score_distribution",
    "highly_variable_genes",
    "highest_expr_genes",
    "dpt_groups_pseudotime",
    "dpt_timeseries",
    "paga",
    "paga_path",
    "paga_compare",
    "correlation_matrix",
    "rank_genes_groups",
    "rank_genes_groups_violin",
    "rank_genes_groups_stacked_violin",
    "rank_genes_groups_heatmap",
    "rank_genes_groups_dotplot",
    "rank_genes_groups_matrixplot",
    "rank_genes_groups_tracksplot",
    "DotPlot",
    "MatrixPlot",
    "StackedViolin",
]

_MODULE_EXPORTS = {
    "_cplearn_visualization_core": {"coremap"},
    "_generic_visualization": {
        "scatter",
        "heatmap",
        "dotplot",
        "tracksplot",
        "violin",
        "stacked_violin",
        "matrixplot",
        "clustermap",
        "ranking",
        "dendrogram",
        "sim",
    },
    "_embedding_visualization_core": {
        "tsne",
        "umap",
        "diffmap",
        "draw_graph",
        "spatial",
        "embedding",
        "embedding_density",
    },
    "_pca_visualization_core": {
        "pca",
        "pca_loadings",
        "pca_variance_ratio",
        "pca_overview",
    },
    "_preprocessing_visualization_core": {
        "scrublet_score_distribution",
        "highly_variable_genes",
        "highest_expr_genes",
    },
    "_trajectory_visualization_core": {
        "dpt_groups_pseudotime",
        "dpt_timeseries",
        "paga",
        "paga_path",
        "paga_compare",
        "correlation_matrix",
    },
    "_differential_expression_core": {
        "rank_genes_groups",
        "rank_genes_groups_violin",
        "rank_genes_groups_stacked_violin",
        "rank_genes_groups_heatmap",
        "rank_genes_groups_dotplot",
        "rank_genes_groups_matrixplot",
        "rank_genes_groups_tracksplot",
    },
    "_advanced_classes_visualization_core": {
        "DotPlot",
        "MatrixPlot",
        "StackedViolin",
    },
}


def __getattr__(name: str) -> Any:
    for module_name, exports in _MODULE_EXPORTS.items():
        if name in exports:
            module = import_module(f"lotus.visualization.{module_name}")
            value = getattr(module, name)
            globals()[name] = value
            return value
    raise AttributeError(f"module 'lotus.visualization' has no attribute {name!r}")
