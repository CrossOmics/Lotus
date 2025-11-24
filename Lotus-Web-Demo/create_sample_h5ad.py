"""
Create a sample .h5ad file for testing the Lotus Web Demo.

.h5ad is the file format used by AnnData (Annotated Data), which is the standard
data structure for single-cell RNA sequencing data in Python. It stores:
- Expression matrix (cells × genes)
- Cell metadata (observations)
- Gene metadata (variables)
- Dimensionality reductions (PCA, UMAP, etc.)
- Layers (raw counts, normalized data, etc.)
"""

import numpy as np
import pandas as pd
from anndata import AnnData
from pathlib import Path

def create_sample_h5ad(output_path='web_demo/sample_data.h5ad', n_cells=500, n_genes=2000):
    """
    Create a sample .h5ad file with synthetic single-cell data.
    
    Parameters:
    -----------
    output_path : str
        Path where to save the .h5ad file
    n_cells : int
        Number of cells to generate
    n_genes : int
        Number of genes to generate
    """
    print(f"Creating sample .h5ad file with {n_cells} cells and {n_genes} genes...")
    
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Generate synthetic count matrix
    # Using negative binomial distribution to simulate real scRNA-seq data
    counts = np.random.negative_binomial(
        n=5, 
        p=0.3, 
        size=(n_cells, n_genes)
    ).astype(np.float32)
    
    # Add cluster-specific expression patterns
    n_clusters = 5
    cells_per_cluster = n_cells // n_clusters
    
    # Define cluster-specific marker genes
    cluster_markers = {
        0: list(range(0, 50)),      # Cluster 0: genes 0-49
        1: list(range(50, 100)),    # Cluster 1: genes 50-99
        2: list(range(100, 150)),   # Cluster 2: genes 100-149
        3: list(range(150, 200)),   # Cluster 3: genes 150-199
        4: list(range(200, 250))    # Cluster 4: genes 200-249
    }
    
    # Create cluster labels
    cluster_labels = []
    cell_types = []
    
    for cluster_id in range(n_clusters):
        n_in_cluster = cells_per_cluster if cluster_id < n_clusters - 1 else n_cells - cluster_id * cells_per_cluster
        
        for i in range(n_in_cluster):
            cluster_labels.append(cluster_id)
            
            # Assign cell types
            if cluster_id == 0:
                cell_types.append('T_cell')
            elif cluster_id == 1:
                cell_types.append('B_cell')
            elif cluster_id == 2:
                cell_types.append('NK_cell')
            elif cluster_id == 3:
                cell_types.append('Monocyte')
            else:
                cell_types.append('Dendritic_cell')
            
            # Enhance marker gene expression for this cluster
            if cluster_id in cluster_markers:
                for gene_idx in cluster_markers[cluster_id]:
                    counts[len(cluster_labels) - 1, gene_idx] *= np.random.uniform(2.0, 5.0)
    
    # Create gene names
    gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]
    
    # Create cell names
    cell_names = [f"Cell_{i:04d}" for i in range(n_cells)]
    
    # Create AnnData object
    adata = AnnData(counts)
    adata.var_names = gene_names
    adata.obs_names = cell_names
    
    # Add cell metadata (observations)
    adata.obs['cluster'] = pd.Categorical(cluster_labels)
    adata.obs['cell_type'] = pd.Categorical(cell_types)
    adata.obs['n_genes'] = (counts > 0).sum(axis=1)  # Number of expressed genes per cell
    adata.obs['total_counts'] = counts.sum(axis=1)    # Total UMI counts per cell
    
    # Add gene metadata (variables)
    adata.var['gene_id'] = gene_names
    adata.var['highly_variable'] = np.random.choice([True, False], size=n_genes, p=[0.1, 0.9])
    
    # Simulate some preprocessing results
    # PCA (first 20 components)
    n_pcs = 20
    pca_coords = np.random.randn(n_cells, n_pcs).astype(np.float32)
    adata.obsm['X_pca'] = pca_coords
    
    # Create a latent representation (combining PCA with some noise)
    latent = np.hstack([pca_coords, np.random.randn(n_cells, 5).astype(np.float32)])
    adata.obsm['X_latent'] = latent
    
    # Save raw counts in a layer
    adata.layers['raw_counts'] = counts.copy()
    
    # Save the file
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    adata.write(output_path)
    
    print(f"✓ Successfully created sample .h5ad file: {output_path}")
    print(f"  - Cells: {adata.n_obs}")
    print(f"  - Genes: {adata.n_vars}")
    print(f"  - Clusters: {len(adata.obs['cluster'].unique())}")
    print(f"  - Cell types: {', '.join(adata.obs['cell_type'].unique())}")
    print(f"  - File size: {output_path.stat().st_size / 1024 / 1024:.2f} MB")
    print(f"\n💡 You can now load this file in the web demo!")
    
    return adata

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Create a sample .h5ad file for testing')
    parser.add_argument('--output', '-o', default='web_demo/sample_data.h5ad',
                       help='Output path for .h5ad file (default: web_demo/sample_data.h5ad)')
    parser.add_argument('--n-cells', type=int, default=500,
                       help='Number of cells (default: 500)')
    parser.add_argument('--n-genes', type=int, default=2000,
                       help='Number of genes (default: 2000)')
    
    args = parser.parse_args()
    
    create_sample_h5ad(
        output_path=args.output,
        n_cells=args.n_cells,
        n_genes=args.n_genes
    )

