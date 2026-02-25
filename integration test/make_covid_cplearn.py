"""
One-off script: load covid.h5ad → ensure X_pca → core_selection(0.30) → save to covid_cplearn.h5ad.

Run from repo root (requires covid.h5ad in repo root):
    python "integration test/make_covid_cplearn.py"

Output: /Users/liuzhenke/Desktop/Lotus/covid_cplearn.h5ad
"""
import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
src_root = project_root / "src"
for p in (str(project_root), str(src_root)):
    if p not in sys.path:
        sys.path.insert(0, p)

COVID_PATH = project_root / "covid.h5ad"
OUT_PATH = project_root / "covid_cplearn.h5ad"
PERCENTAGE = 0.30

def _ensure_pca(adata):
    if "X_pca" in adata.obsm:
        return adata
    import lotus.preprocessing as pp
    pp.normalize_total(adata, target_sum=1e4)
    pp.log1p(adata)
    pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
    adata = adata[:, adata.var["highly_variable"]].copy()
    pp.scale(adata, max_value=10)
    pp.pca(adata, n_comps=50)
    return adata


def main():
    import lotus.io as io
    from lotus.core_analysis import core_selection

    if not COVID_PATH.exists():
        print(f"Error: {COVID_PATH} not found. Put covid.h5ad in repo root and re-run.")
        sys.exit(1)

    print(f"Loading {COVID_PATH.name} ...")
    adata = io.standardize_load(COVID_PATH)
    print(f"  {adata.n_obs} cells × {adata.n_vars} genes")
    adata = _ensure_pca(adata)
    print(f"Running core_selection(percentage={PERCENTAGE}) ...")
    adata, _ = core_selection(adata, PERCENTAGE)
    print(f"Writing {OUT_PATH} ...")
    io.write(OUT_PATH, adata)
    print(f"Done. {OUT_PATH}")


if __name__ == "__main__":
    main()
