# UMAP Visualization Example Image

This directory should contain example visualization images for the documentation.

## To generate the actual UMAP visualization image:

### Option 1: Run the example script
```bash
python examples/lotus_workflow.py --clusters 3 --cells-per-cluster 60
cp result_*/umap_clusters.png docs/_static/umap_clusters_example.png
```

### Option 2: Use the placeholder script
```bash
python docs/_static/create_image_placeholder.py
```

The image should show a UMAP plot with cells colored by cluster labels, demonstrating the clustering results from the Lotus workflow.

**Note**: If the image doesn't exist, the documentation will display a helpful message with instructions on how to generate it.

