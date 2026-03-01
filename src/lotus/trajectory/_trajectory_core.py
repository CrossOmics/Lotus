
from pathlib import Path
from typing import Iterable, Literal

import numpy as np
import scanpy as sc
from anndata import AnnData
from scipy.sparse import csc_matrix, csr_matrix

from lotus.lineagetracker import logged


@logged
def diffmap(
        adata: AnnData,
        n_comps: int = 15,
        *,
        neighbors_key: str | None = None,
        random_state: int | np.random.RandomState | None = 0,
        copy: bool = False,
) -> AnnData | None:
    """
    Compute a diffusion map embedding.

    Wraps scanpy.tl.diffmap.
    """
    return sc.tl.diffmap(
        adata,
        n_comps=n_comps,
        neighbors_key=neighbors_key,
        random_state=random_state,
        copy=copy,
    )


@logged
def draw_graph(
        adata: AnnData,
        layout: Literal['fr', 'drl', 'kk', 'grid_fr', 'lgl', 'rt', 'rt_circular', 'fa'] = 'fa',
        *,
        init_pos: str | bool | None = None,
        root: int | None = None,
        random_state: int | np.random.RandomState | None = 0,
        n_jobs: int | None = None,
        adjacency: csr_matrix | csc_matrix | None = None,
        key_added_ext: str | None = None,
        neighbors_key: str | None = None,
        obsp: str | None = None,
        copy: bool = False,
        **kwds,
) -> AnnData | None:
    """
    Compute a force-directed graph layout.

    Wraps scanpy.tl.draw_graph.
    """
    return sc.tl.draw_graph(
        adata,
        layout=layout,
        init_pos=init_pos,
        root=root,
        random_state=random_state,
        n_jobs=n_jobs,
        adjacency=adjacency,
        key_added_ext=key_added_ext,
        neighbors_key=neighbors_key,
        obsp=obsp,
        copy=copy,
        **kwds,
    )


@logged
def dpt(
        adata: AnnData,
        n_dcs: int = 10,
        *,
        n_branchings: int = 0,
        min_group_size: float = 0.01,
        allow_kendall_tau_shift: bool = True,
        neighbors_key: str | None = None,
        copy: bool = False,
) -> AnnData | None:
    """
    Infer diffusion pseudotime.

    Wraps scanpy.tl.dpt.
    """
    return sc.tl.dpt(
        adata,
        n_dcs=n_dcs,
        n_branchings=n_branchings,
        min_group_size=min_group_size,
        allow_kendall_tau_shift=allow_kendall_tau_shift,
        neighbors_key=neighbors_key,
        copy=copy,
    )


@logged
def paga(
        adata: AnnData,
        groups: str | None = None,
        *,
        use_rna_velocity: bool = False,
        model: Literal['v1.2', 'v1.0'] = 'v1.2',
        neighbors_key: str | None = None,
        copy: bool = False,
) -> AnnData | None:
    """
    Compute a PAGA graph over groups.

    Wraps scanpy.tl.paga.
    """
    return sc.tl.paga(
        adata,
        groups=groups,
        use_rna_velocity=use_rna_velocity,
        model=model,
        neighbors_key=neighbors_key,
        copy=copy,
    )


@logged
def embedding_density(
        adata: AnnData,
        basis: str = 'umap',
        *,
        groupby: str | None = None,
        key_added: str | None = None,
        components: str | Iterable[str] | None = None,
) -> None:
    """
    Calculate density values over an embedding.

    Wraps scanpy.tl.embedding_density.
    """
    sc.tl.embedding_density(
        adata,
        basis=basis,
        groupby=groupby,
        key_added=key_added,
        components=components,
    )


@logged
def sim(
        model: Literal['krumsiek11', 'toggleswitch'],
        *,
        params_file: bool = True,
        tmax: int | None = None,
        branching: bool | None = None,
        nrRealizations: int | None = None,
        noiseObs: float | None = None,
        noiseDyn: float | None = None,
        step: int | None = None,
        seed: int | None = None,
        writedir: Path | str | None = None,
) -> AnnData:
    """
    Simulate dynamic gene expression data.

    Wraps scanpy.tl.sim.
    """
    return sc.tl.sim(
        model,
        params_file=params_file,
        tmax=tmax,
        branching=branching,
        nrRealizations=nrRealizations,
        noiseObs=noiseObs,
        noiseDyn=noiseDyn,
        step=step,
        seed=seed,
        writedir=writedir,
    )
