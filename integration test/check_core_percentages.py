"""
Load covid_cplearn.h5ad and run core_selection at 10%, 30%, 50%.
Report for each: n_core, expected by proportion, and which layers are included.

Run from repo root:
    python "integration test/check_core_percentages.py"
"""
import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
src_root = project_root / "src"
for p in (str(project_root), str(src_root)):
    if p not in sys.path:
        sys.path.insert(0, p)

DATA_PATH = project_root / "covid_cplearn.h5ad"
PERCENTAGES = (0.10, 0.30, 0.50)


def main():
    import lotus.io as io
    from lotus.core_analysis import core_selection
    import numpy as np

    if not DATA_PATH.exists():
        print(f"Error: {DATA_PATH} not found. Run make_covid_cplearn.py first.")
        sys.exit(1)

    print(f"Loading {DATA_PATH.name} ...")
    adata = io.standardize_load(DATA_PATH)
    n_obs = adata.n_obs
    layer_key = "cplearn_layer"
    print(f"  n_obs = {n_obs:,}\n")

    print("Percentage | n_core (actual) | min target (by %) | layers in core (whole layers only)")
    print("-" * 75)

    for pct in PERCENTAGES:
        adata_run, _ = core_selection(adata, pct)
        sel = adata_run.obs["cplearn_core_selection"] == "core"
        n_core = int(sel.sum())
        target_min = int(np.ceil(pct * n_obs))
        core_layers = adata_run.obs.loc[sel, layer_key]
        layer_counts = core_layers.value_counts().sort_index()
        layer_str = ", ".join(f"L{k}({v})" for k, v in layer_counts.items())

        print(f"  {pct*100:5.0f}%   | {n_core:>14,} | {target_min:>16,} | {layer_str}")

    print("-" * 75)
    print("(Whole-layer policy: actual >= min target; we never split a layer.)")


if __name__ == "__main__":
    main()
