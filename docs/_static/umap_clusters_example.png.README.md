# UMAP Visualization Example Image

This directory should contain example visualization images for the documentation.

To generate the actual UMAP visualization image:

1. Run the example script:
   ```bash
   python examples/lotus_workflow.py --clusters 3 --cells-per-cluster 60
   ```

2. Copy the generated image from the output directory:
   ```bash
   cp result_*/umap_clusters.png docs/_static/umap_clusters_example.png
   ```

The image should show a UMAP plot with cells colored by cluster labels, demonstrating the clustering results from the Lotus workflow.

