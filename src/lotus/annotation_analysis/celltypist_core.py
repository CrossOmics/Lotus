import celltypist
from celltypist import models
import pandas as pd
import anndata
from typing import Dict, Any, Optional, List

# Default model list
AVAILABLE_MODELS = [
    "Immune_All_Low.pkl",
    "Immune_All_High.pkl",
    "Human_Lung_Atlas.pkl",
    "Healthy_COVID19_PBMC.pkl"
]


def run_celltypist_annotation(
        adata: anndata.AnnData,
        model_name: str = "Immune_All_Low.pkl",
        majority_voting: bool = True,
        target_cluster_col: Optional[str] = None
) -> Dict[str, Any]:
    """
    Core function to run CellTypist annotation on an AnnData object.

    Args:
        adata: AnnData object containing the expression matrix.
               CellTypist handles normalization internally (expects log1p or raw counts).
        model_name: Name of the pre-trained model to use (default: Immune_All_Low).
        majority_voting: Whether to refine predictions using cluster-based majority voting.
                         Requires 'target_cluster_col' to be present in adata.obs.
        target_cluster_col: The column name in adata.obs representing clusters (e.g., 'leiden').
                            Required if majority_voting is True.

    Returns:
        Dictionary containing:
            - "model_used": Name of the model
            - "cluster_summary": Dictionary summarizing predictions per cluster
            - "cell_summary": (Optional) Detailed per-cell predictions (can be large)
    """

    results = {
        "model_used": model_name,
        "parameters": {
            "majority_voting": majority_voting,
            "target_cluster_col": target_cluster_col
        }
    }

    # Validation & Model Loading
    try:
        # Check if model exists in cache or download it
        # models.models_path contains the local cache directory
        model = models.Model.load(model=model_name)
    except Exception as e:
        return {"error": f"Failed to load CellTypist model '{model_name}': {str(e)}"}

    # Run Annotation
    # annotate() returns an AnnotationResult object
    try:
        # Note: 'mode' can be 'best match' or 'prob match'. 'best match' is standard.
        predictions = celltypist.annotate(
            adata,
            model=model,
            majority_voting=majority_voting,
            mode='best match'
        )
    except Exception as e:
        return {"error": f"CellTypist annotation process failed: {str(e)}"}

    # Process Results
    # Convert result to AnnData to easily manipulate obs
    adata_pred = predictions.to_adata()

    # Determine which column holds the final prediction
    # If majority_voting is True, CellTypist adds a 'majority_voting' column
    # Otherwise, the primary column is 'predicted_labels'
    pred_col = "majority_voting" if majority_voting and "majority_voting" in adata_pred.obs else "predicted_labels"

    # Generate Cluster-Level Summary (Aggregation)
    # This is critical for the frontend to show "Cluster 0 is T cells"
    if target_cluster_col and target_cluster_col in adata.obs:
        cluster_summary = {}

        # We need to map predictions back to the original cells groupings
        # Create a temp dataframe aligning the cluster info with prediction info
        combined_obs = pd.DataFrame({
            "cluster": adata.obs[target_cluster_col].values,
            "prediction": adata_pred.obs[pred_col].values,
            # Confidence score is usually stored in 'conf_score'
            "conf_score": adata_pred.obs["conf_score"].values if "conf_score" in adata_pred.obs else 0.0
        }, index=adata.obs.index)

        # Group by Cluster
        for cluster_id, group in combined_obs.groupby("cluster"):
            # A. Find dominant cell type (Mode)
            if group.empty:
                continue

            top_label = group["prediction"].mode()[0]

            # B. Calculate statistics
            total_cells = len(group)
            top_count = len(group[group["prediction"] == top_label])
            fraction = top_count / total_cells  # "Purity" or "Agreement"

            # C. Average Confidence (of the assigned label)
            avg_conf = group["conf_score"].mean()

            cluster_summary[str(cluster_id)] = {
                "predicted_cell_type": top_label,
                "cell_count": int(total_cells),
                "agreement_fraction": float(f"{fraction:.4f}"),  # e.g. 0.95
                "average_confidence": float(f"{avg_conf:.4f}")
            }

        results["cluster_summary"] = cluster_summary

    # Cell-level detailed results
    # not return this in the JSON if dataset is huge (100k cells)
    results["cell_predictions"] = adata_pred.obs[pred_col].to_dict()

    return results


def download_models_to_cache(model_list: List[str] = None):
    """
    Helper utility to pre-download models to the server cache on startup.
    This prevents timeouts on the first user request.
    """
    targets = model_list if model_list else AVAILABLE_MODELS
    print(f"Checking/Downloading CellTypist models: {targets}")
    for m in targets:
        models.download_model(m)
