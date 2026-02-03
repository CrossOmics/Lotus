import gseapy as gp
import pandas as pd
from typing import List, Dict, Optional, Any

# Define databases to query
# mapping user-friendly names to Enrichr library names
GSEAPY_DATABASES = {
    "cellmarker": ["CellMarker_Augmented_2021"],
    "functional": ["GO_Biological_Process_2023", "KEGG_2021_Human"],
    "msigdb": ["MSigDB_Hallmark_2020"],
    "reactome": ["Reactome_2022"],
    "wikipathway": ["WikiPathway_2021_Human"],
    "tissue": ["Human_Gene_Atlas"],
    "tf": ["Enrichr_Submissions_TF-Gene_Coocurrence"],
}


def _extract_cell_type_from_term(term_name: str, p_value: float) -> Optional[str]:
    """
    Helper: Parses CellMarker terms (usually 'CellType:Tissue').
    Returns the cell type if p-value is significant (< 0.05).
    """
    if not term_name or pd.isna(term_name):
        return None
    if not p_value or pd.isna(p_value) or float(p_value) > 0.05:
        return None

    term = str(term_name).strip()
    # If format is "CellType:Tissue", take the left part
    if ":" in term:
        return term.split(":", 1)[0].strip()
    return None


def run_enrichr_analysis(gene_list: List[str], categories: List[str] = None, limit=None) -> Dict[str, Any]:
    """
    Core function to run GSEApy Enrichr on a list of genes.

    Args:
        gene_list: List of gene symbols (e.g., top 100 markers).
        categories: List of categories to query (e.g., ['cellmarker', 'functional']).
                   If None, queries all defined databases.
        limit: Keep top N results

    Returns:
        Dictionary containing results for each requested category.
    """
    results = {}

    # If no categories specified, use all
    target_cats = categories if categories else GSEAPY_DATABASES.keys()

    for category in target_cats:
        if category not in GSEAPY_DATABASES:
            continue

        db_libraries = GSEAPY_DATABASES[category]

        try:
            # Call GSEApy Enrichr API
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=db_libraries,
                organism="human",
                outdir=None,  # Don't save to disk
                verbose=False
            )

            # Extract results dataframe
            if enr is None or not hasattr(enr, "res2d") or enr.res2d is None:
                results[category] = None
                continue

            df = enr.res2d.copy()

            # Filter by P-value < 0.05 for relevance
            if "Adjusted P-value" in df.columns:
                df = df[df["Adjusted P-value"] < 0.05]
            elif "P-value" in df.columns:
                df = df[df["P-value"] < 0.05]


            # Keep top "limit" results to keep response light
            if limit:
                df_top = df.head(limit)
            else:
                df_top = df
            # Convert to dictionary for JSON serialization
            category_result = {
                "top_terms": df_top.to_dict("records")
            }

            # Special logic for CellMarker: Infer cell types
            if category == "cellmarker":
                inferred_types = set()
                term_col = "Term" if "Term" in df_top.columns else df_top.columns[0]
                pval_col = "Adjusted P-value" if "Adjusted P-value" in df_top.columns else "P-value"

                for _, row in df_top.iterrows():
                    ct = _extract_cell_type_from_term(row.get(term_col), row.get(pval_col))
                    if ct:
                        inferred_types.add(ct)

                category_result["inferred_cell_types"] = sorted(list(inferred_types))

            results[category] = category_result

        except Exception as e:
            print(f"Error running Enrichr for {category}: {e}")
            results[category] = {"error": str(e)}

    return results
