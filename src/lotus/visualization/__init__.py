from ._cplearn_visualization_core import coremap
from ._generic_visualization import (
    scatter,
    heatmap,
    dotplot,
    tracksplot,
    violin,
    stacked_violin,
    matrixplot,
    clustermap,
    ranking,
    dendrogram,
)
from ._embedding_visualization_core import (
    # Embedding Plots
    tsne,
    umap,
    diffmap,
    draw_graph,
    spatial,
    embedding,
    embedding_density,
)

from ._pca_visualization_core import (
    pca,
    pca_loadings,
    pca_variance_ratio,
    pca_overview,
)

from ._preprocessing_visualization_core import (
    scrublet_score_distribution,
    highly_variable_genes,
    highest_expr_genes,
)

from ._trajectory_visualization_core import (
    dpt_groups_pseudotime,
    dpt_timeseries,
    paga,
    paga_path,
    paga_compare,
    correlation_matrix,
)

from ._differential_expression_core import (
    rank_genes_groups,
    rank_genes_groups_violin,
    rank_genes_groups_stacked_violin,
    rank_genes_groups_heatmap,
    rank_genes_groups_dotplot,
    rank_genes_groups_matrixplot,
    rank_genes_groups_tracksplot,
)

from ._advanced_classes_visualization_core import (
    DotPlot,
    MatrixPlot,
    StackedViolin,
)

__all__ = [
    # Cplearn Plots
    "coremap",
    # Generic Plots
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
    "coremap",

    # Embedding Plots
    "tsne",
    "umap",
    "diffmap",
    "draw_graph",
    "spatial",
    "embedding",
    "embedding_density",
    # PCA Plots
    "pca",
    "pca_loadings",
    "pca_variance_ratio",
    "pca_overview",
    # Preprocessing Plots
    "scrublet_score_distribution",
    "highly_variable_genes",
    "highest_expr_genes",
    # Trajectory
    "dpt_groups_pseudotime",
    "dpt_timeseries",
    "paga",
    "paga_path",
    "paga_compare",
    "correlation_matrix",
    #  Differential Expression / Maker genes
    "rank_genes_groups",
    "rank_genes_groups_violin",
    "rank_genes_groups_stacked_violin",
    "rank_genes_groups_heatmap",
    "rank_genes_groups_dotplot",
    "rank_genes_groups_matrixplot",
    "rank_genes_groups_tracksplot",
    # advanced classes
    "DotPlot",
    "MatrixPlot",
    "StackedViolin",
]
