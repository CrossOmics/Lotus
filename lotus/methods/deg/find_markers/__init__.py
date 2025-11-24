from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import stats as _scipy_stats
from scipy.stats import mannwhitneyu, rankdata
from statsmodels.stats.multitest import multipletests
from scipy import sparse
from tabulate import tabulate

ArrayLike = Union[np.ndarray, Sequence[int]]
Label = Union[int, str]

__all__ = ["DEOptions", "DEAnalyzer", "de_from_model", "de_from_adata"]


@dataclass
class DEOptions:
    pseudocount: float = 1.0
    min_cells_per_group: int = 3
    min_detect_pct: float = 0.05
    use_exact: bool = False
    alternative: str = "two-sided"
    gene_filter_min_nonzero: int = 0


class DEAnalyzer:
    def __init__(
        self,
        corespect_obj: Any,
        genes: Optional[Sequence[str]] = None,
        options: Optional[DEOptions] = None,
    ) -> None:
        if corespect_obj.count_mat is None:
            cell_x_feat = corespect_obj.X
        else:
            cell_x_feat = corespect_obj.count_mat

        labels = corespect_obj.labels_
        layers = {i + 1: np.array(layer, dtype=int) for i, layer in enumerate(corespect_obj.layers_ or [])}

        self.layers_ = corespect_obj.layers_

        cell_x_feat = np.asarray(cell_x_feat)
        if cell_x_feat.ndim != 2:
            raise ValueError("count matrix must be a 2D array (n_cells, n_features)")
        self.X = cell_x_feat
        self.n_cells, self.n_genes = cell_x_feat.shape

        self.labels = np.asarray(labels)
        if self.labels.shape[0] != self.n_cells:
            raise ValueError("labels length must match number of cells:")

        self.layers = {k: np.asarray(v, dtype=int) for k, v in layers.items()}
        self._validate_layers_disjoint()

        if genes is None:
            self.genes = np.array([f"g{i}" for i in range(self.n_genes)], dtype=object)
        else:
            if len(genes) != self.n_genes:
                raise ValueError("genes length must match number of columns in X")
            self.genes = np.array(genes, dtype=object)

        self.options = options or DEOptions()

        if self.options.gene_filter_min_nonzero > 0:
            nnz = np.count_nonzero(self.X, axis=0)
            self._gene_keep_mask = nnz >= self.options.gene_filter_min_nonzero
        else:
            self._gene_keep_mask = np.ones(self.n_genes, dtype=bool)

        self.X_masked = np.ascontiguousarray(self.X[:, self._gene_keep_mask])
        self.genes_masked = self.genes[self._gene_keep_mask]

    def _validate_layers_disjoint(self) -> None:
        seen = set()
        for idx in self.layers.values():
            if idx.size == 0:
                continue
            if np.any(idx < 0) or np.any(idx >= self.n_cells):
                raise ValueError("layer has out-of-range indices")
            overlap = seen.intersection(idx.tolist())
            if overlap:
                raise ValueError("layers are not disjoint; overlap detected")
            seen.update(idx.tolist())

    @staticmethod
    def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
        _, p_adj, _, _ = multipletests(pvals, method="fdr_bh")
        return p_adj

    def _test_two_groups(
        self,
        idx_a: np.ndarray,
        idx_b: np.ndarray,
        meta: Mapping[str, Union[str, int, float]],
    ) -> pd.DataFrame:
        opts = self.options

        if idx_a.size < opts.min_cells_per_group or idx_b.size < opts.min_cells_per_group:
            return pd.DataFrame(
                columns=[
                    "gene",
                    "statistic",
                    "pvalue",
                    "p_adj",
                    "log2fc",
                    "mean_a",
                    "mean_b",
                    "pct_expr_a",
                    "pct_expr_b",
                    "z_score",
                    *meta.keys(),
                ]
            )

        X_a = self.X_masked.take(idx_a, axis=0)
        X_b = self.X_masked.take(idx_b, axis=0)
        genes_kept = self.genes_masked

        n_g = genes_kept.size
        stats = np.empty(n_g, dtype=float)
        pvals = np.empty(n_g, dtype=float)

        mean_a = X_a.mean(axis=0)
        mean_b = X_b.mean(axis=0)
        pct_a = (X_a > 0).mean(axis=0)
        pct_b = (X_b > 0).mean(axis=0)

        if opts.min_detect_pct > 0.0:
            keep2 = np.logical_or(pct_a > opts.min_detect_pct, pct_b > opts.min_detect_pct)
            if not np.all(keep2):
                X_a = X_a[:, keep2]
                X_b = X_b[:, keep2]
                genes_kept = genes_kept[keep2]
                mean_a = mean_a[keep2]
                mean_b = mean_b[keep2]
                pct_a = pct_a[keep2]
                pct_b = pct_b[keep2]
                n_g = genes_kept.size
                stats = np.empty(n_g, dtype=float)
                pvals = np.empty(n_g, dtype=float)

        for g in range(n_g):
            try:
                res = mannwhitneyu(
                    X_a[:, g],
                    X_b[:, g],
                    alternative=opts.alternative,
                    method="exact" if opts.use_exact else "asymptotic",
                )
                stats[g] = res.statistic
                pvals[g] = res.pvalue
            except Exception:
                stats[g] = np.nan
                pvals[g] = 1.0

        X_all = np.vstack([X_a, X_b])
        n1, n2 = X_a.shape[0], X_b.shape[0]
        n = n1 + n2
        ranks = np.apply_along_axis(rankdata, 0, X_all)
        R1 = ranks[:n1, :].sum(axis=0)
        U = R1 - n1 * (n1 + 1) / 2.0
        stats = U
        mean_U = n1 * n2 / 2.0
        std_U = np.sqrt(n1 * n2 * (n + 1) / 12.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            z = (U - mean_U) / std_U
        pvals = 2 * _scipy_stats.norm.sf(np.abs(z))

        with np.errstate(divide="ignore", invalid="ignore"):
            r_effect = 1.0 - (2.0 * stats) / (n1 * n2)
        pc = opts.pseudocount
        log2fc = np.log2((mean_a + pc) / (mean_b + pc))

        n_a, n_b = len(idx_a), len(idx_b)
        mean_diff = mean_a - mean_b
        std_a = X_a.std(axis=0, ddof=1)
        std_b = X_b.std(axis=0, ddof=1)
        se = np.sqrt((std_a**2) / n_a + (std_b**2) / n_b)
        z_score = np.divide(mean_diff, se, out=np.zeros_like(mean_diff), where=se > 0)

        df = pd.DataFrame(
            {
                "gene": genes_kept,
                "statistic": stats,
                "pvalue": pvals,
                "log2fc": log2fc,
                "effect_r": r_effect,
                "mean_a": mean_a,
                "mean_b": mean_b,
                "pct_expr_a": pct_a,
                "pct_expr_b": pct_b,
                "z_score": z_score,
            }
        )
        df["p_adj"] = self._bh_fdr(df["pvalue"].values)
        for key, value in meta.items():
            df[key] = value
        df = self._sort_with_pvalue_buckets(df)
        return df

    def run_default(
        self,
        set1: Iterable[Label],
        set2: Optional[Iterable[Label]] = None,
        ignore_label: Label = -1,
        name1: Optional[str] = None,
        name2: Optional[str] = None,
        starting_layer: Optional[int] = None,
        ending_layer: Optional[int] = None,
    ) -> pd.DataFrame:
        set1 = set(set1)
        if set2 is None:
            others = set(np.unique(self.labels)) - set1 - {ignore_label}
            set2 = others
        else:
            set2 = set(set2)

        if starting_layer is not None or ending_layer is not None:
            start = starting_layer or 0
            end = ending_layer or len(self.layers_)
            selected_indices = np.concatenate(self.layers_[start:end])
        else:
            selected_indices = np.arange(len(self.labels))

        labels_subset = np.asarray(self.labels)[selected_indices]
        idx_a_local = np.isin(labels_subset, list(set1))
        idx_b_local = np.isin(labels_subset, list(set2))

        idx_a = selected_indices[idx_a_local]
        idx_b = selected_indices[idx_b_local]

        meta: Dict[str, Union[str, int, float, None]] = {
            "mode": "default",
            "group_a": name1 if name1 is not None else f"{sorted(list(set1))}",
            "group_b": name2 if name2 is not None else f"{sorted(list(set2))}",
            "starting_layer": starting_layer,
            "ending_layer": ending_layer,
        }
        return self._test_two_groups(idx_a, idx_b, meta)

    @staticmethod
    def _sort_with_pvalue_buckets(df: pd.DataFrame) -> pd.DataFrame:
        bins = [0, 0.001, 0.01, 0.05, np.inf]
        labels = [1, 2, 3, 4]
        df = df.copy()
        df["abs_log2fc"] = df["log2fc"].abs()
        df["p_bucket"] = pd.cut(df["p_adj"], bins=bins, labels=labels, include_lowest=True)
        df["p_bucket"] = df["p_bucket"].astype(int)
        return df.sort_values(["p_bucket", "abs_log2fc"], ascending=[True, False], ignore_index=True)

    @staticmethod
    def observe_de(
        df: pd.DataFrame,
        n_rows: int = 10,
        fmt: str = "github",
        genes: Optional[Sequence[str]] = None,
    ) -> None:
        if df.empty:
            print("(no results)")
            return
        if genes is not None:
            df = df[df["gene"].isin(genes)]
            if df.empty:
                print(f"(no results for selected genes: {genes})")
                return
        cols = [
            "gene",
            "log2fc",
            "z_score",
            "mean_a",
            "mean_b",
            "pct_expr_a",
            "pct_expr_b",
            "pvalue",
        ]
        cols = [c for c in cols if c in df.columns]
        df_show = df[cols].copy()
        df_show = df_show.round(
            {
                "log2fc": 3,
                "z_score": 3,
                "mean_a": 3,
                "mean_b": 3,
                "pct_expr_a": 3,
                "pct_expr_b": 3,
                "pvalue": 3,
            }
        )
        print(tabulate(df_show.head(n_rows), headers="keys", tablefmt=fmt, showindex=False))


def de_from_model(
    corespect_model: Any,
    *,
    genes: Optional[Sequence[str]] = None,
    options: Optional[DEOptions] = None,
    groups_a: Iterable[Label],
    groups_b: Optional[Iterable[Label]] = None,
    **kwargs: Any,
) -> pd.DataFrame:
    analyzer = DEAnalyzer(corespect_model, genes=genes, options=options)
    return analyzer.run_default(groups_a, groups_b, **kwargs)


def _get_array_from_adata(
    adata: AnnData,
    *,
    layer: str | None = None,
    use_rep: str | None = None,
) -> np.ndarray:
    if layer is not None and use_rep is not None:
        raise ValueError("Specify at most one of `layer` and `use_rep`.")
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


class _AnnDataCorespectProxy:
    def __init__(self, matrix: np.ndarray, labels: np.ndarray, layers: Sequence[Sequence[int]]) -> None:
        self.count_mat = matrix
        self.X = matrix
        self.labels_ = labels
        self.layers_ = [np.asarray(layer, dtype=int) for layer in layers]


def de_from_adata(
    adata: AnnData,
    *,
    key: str = "cplearn",
    layer: str | None = None,
    use_rep: str | None = None,
    layers_key: str | None = None,
    genes: Optional[Sequence[str]] = None,
    options: Optional[DEOptions] = None,
    groups_a: Iterable[Label],
    groups_b: Optional[Iterable[Label]] = None,
    **kwargs: Any,
) -> pd.DataFrame:
    """
    Run differential expression using AnnData annotated with CoreSpect results.

    Parameters
    ----------
    adata
        Annotated data matrix.
    key
        Column in ``adata.obs`` containing cluster labels (default ``"cplearn"``).
    layer, use_rep
        Optional data sources to override ``adata.X``.
    layers_key
        Key inside ``adata.uns`` that stores CoreSpect layer indices. Defaults to ``f\"{key}_cplearn\"``.
    genes
        Optional explicit list of gene names. Defaults to ``adata.var_names`` when available.
    options
        Differential expression options.
    groups_a, groups_b
        Groups to compare. These mirror :meth:`DEAnalyzer.run_default`.
    kwargs
        Passed through to :meth:`DEAnalyzer.run_default`.
    """

    if key not in adata.obs:
        raise KeyError(f"{key!r} not found in adata.obs.")
    labels_series = adata.obs[key]
    if isinstance(labels_series, pd.Series):
        labels = labels_series.to_numpy()
    else:
        labels = np.asarray(labels_series)
    labels = labels.astype(object)

    matrix = _get_array_from_adata(adata, layer=layer, use_rep=use_rep)

    layers_meta_key = layers_key or f"{key}_cplearn"
    if layers_meta_key not in adata.uns:
        raise KeyError(f"{layers_meta_key!r} not found in adata.uns.")
    layers_meta = adata.uns[layers_meta_key]
    raw_layers = layers_meta.get("layers")
    if raw_layers is None:
        raise KeyError(f"Layer information missing under adata.uns[{layers_meta_key!r}]['layers'].")
    layers_list = [np.asarray(layer, dtype=int) for layer in raw_layers]
    if not layers_list:
        layers_list = [np.arange(matrix.shape[0], dtype=int)]

    proxy = _AnnDataCorespectProxy(matrix, labels, layers_list)

    if genes is None and adata.var_names is not None:
        genes = list(map(str, adata.var_names))

    analyzer = DEAnalyzer(proxy, genes=genes, options=options)
    return analyzer.run_default(groups_a, groups_b, **kwargs)

