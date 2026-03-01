from collections.abc import Iterable

import scanpy as sc
from anndata import AnnData

from lotus.lineagetracker import logged


@logged
def ingest(
        adata: AnnData,
        adata_ref: AnnData,
        *,
        obs: str | Iterable[str] | None = None,
        embedding_method: str | Iterable[str] = ('umap', 'pca'),
        labeling_method: str = 'knn',
        neighbors_key: str | None = None,
        inplace: bool = True,
        **kwargs,
) -> AnnData | None:
    """
    Map labels and embeddings from a reference dataset to new data.

    Wraps scanpy.tl.ingest.
    """
    return sc.tl.ingest(
        adata,
        adata_ref,
        obs=obs,
        embedding_method=embedding_method,
        labeling_method=labeling_method,
        neighbors_key=neighbors_key,
        inplace=inplace,
        **kwargs,
    )
