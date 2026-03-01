from collections.abc import Collection, Iterable
from typing import Literal

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from lotus.lineagetracker import logged


@logged
def obs_df(
        adata: AnnData,
        keys: Collection[str] = (),
        obsm_keys: Iterable[tuple[str, int]] = (),
        *,
        layer: str | None = None,
        gene_symbols: str | None = None,
        use_raw: bool = False,
) -> pd.DataFrame:
    """
    Return observation-aligned values from an AnnData object.

    Wraps scanpy.get.obs_df.
    """
    return sc.get.obs_df(
        adata,
        keys=keys,
        obsm_keys=obsm_keys,
        layer=layer,
        gene_symbols=gene_symbols,
        use_raw=use_raw,
    )


@logged
def var_df(
        adata: AnnData,
        keys: Collection[str] = (),
        varm_keys: Iterable[tuple[str, int]] = (),
        *,
        layer: str | None = None,
) -> pd.DataFrame:
    """
    Return variable-aligned values from an AnnData object.

    Wraps scanpy.get.var_df.
    """
    return sc.get.var_df(
        adata,
        keys=keys,
        varm_keys=varm_keys,
        layer=layer,
    )


@logged
def aggregate(
        adata: AnnData,
        by: str | Collection[str],
        func: Literal['count_nonzero', 'mean', 'sum', 'var', 'median'] | Iterable[
            Literal['count_nonzero', 'mean', 'sum', 'var', 'median']
        ],
        *,
        axis: Literal['obs', 0, 'var', 1] | None = None,
        mask: np.ndarray | str | None = None,
        dof: int = 1,
        layer: str | None = None,
        obsm: str | None = None,
        varm: str | None = None,
) -> AnnData:
    """
    Aggregate expression values over categorical groups.

    Wraps scanpy.get.aggregate.
    """
    return sc.get.aggregate(
        adata,
        by=by,
        func=func,
        axis=axis,
        mask=mask,
        dof=dof,
        layer=layer,
        obsm=obsm,
        varm=varm,
    )
