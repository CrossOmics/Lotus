"""
Create sample data files in .h5, .csv, and .tsv formats for testing.

This script generates:
- .h5: 10x Genomics H5 format (HDF5)
- .csv: CSV matrix format (genes as rows, cells as columns)
- .tsv: TSV matrix format (genes as rows, cells as columns)
"""

import numpy as np
import pandas as pd
from pathlib import Path
import h5py
import sys

# Add parent directory to path to import lotus if needed
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    import scanpy as sc
    SCANPY_AVAILABLE = True
except ImportError:
    SCANPY_AVAILABLE = False
    print("Warning: scanpy not available. Some formats may not be generated.")


def create_sample_data(output_dir='data', n_cells=500, n_genes=2000):
    """
    Create sample data files in multiple formats.
    
    Parameters:
    -----------
    output_dir : str
        Directory where to save the data files
    n_cells : int
        Number of cells to generate
    n_genes : int
        Number of genes to generate
    """
    print(f"Creating sample data files with {n_cells} cells and {n_genes} genes...")
    
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
    
    # Enhance marker gene expression for each cluster
    for cluster_id in range(n_clusters):
        start_idx = cluster_id * cells_per_cluster
        end_idx = start_idx + cells_per_cluster if cluster_id < n_clusters - 1 else n_cells
        
        if cluster_id in cluster_markers:
            for gene_idx in cluster_markers[cluster_id]:
                counts[start_idx:end_idx, gene_idx] *= np.random.uniform(2.0, 5.0, size=end_idx-start_idx)
    
    # Create gene names
    gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]
    
    # Create cell names
    cell_names = [f"Cell_{i:04d}" for i in range(n_cells)]
    
    # Ensure output directory exists
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Create CSV file (genes as rows, cells as columns)
    print("\n1. Creating CSV file...")
    csv_path = output_dir / 'sample_data.csv'
    df_csv = pd.DataFrame(counts.T, index=gene_names, columns=cell_names)
    df_csv.to_csv(csv_path)
    print(f"   ✓ Created: {csv_path}")
    print(f"   - Shape: {df_csv.shape} (genes × cells)")
    print(f"   - Size: {csv_path.stat().st_size / 1024 / 1024:.2f} MB")
    
    # 2. Create TSV file (genes as rows, cells as columns)
    print("\n2. Creating TSV file...")
    tsv_path = output_dir / 'sample_data.tsv'
    df_tsv = pd.DataFrame(counts.T, index=gene_names, columns=cell_names)
    df_tsv.to_csv(tsv_path, sep='\t')
    print(f"   ✓ Created: {tsv_path}")
    print(f"   - Shape: {df_tsv.shape} (genes × cells)")
    print(f"   - Size: {tsv_path.stat().st_size / 1024 / 1024:.2f} MB")
    
    # 3. Create H5 file (10x Genomics format)
    print("\n3. Creating H5 file (10x Genomics format)...")
    h5_path = output_dir / 'sample_data.h5'
    
    if SCANPY_AVAILABLE:
        # Use scanpy to create proper 10x H5 format
        try:
            # Create AnnData object
            from anndata import AnnData
            adata = AnnData(counts)
            adata.var_names = gene_names
            adata.obs_names = cell_names
            
            # Convert to sparse matrix for efficiency
            import scipy.sparse as sp
            adata.X = sp.csr_matrix(adata.X)
            
            # Save as 10x H5 format
            # Note: scanpy's write_10x_h5 requires specific structure
            # We'll create a simplified H5 format that can be read
            with h5py.File(h5_path, 'w') as f:
                # Create matrix group (10x format structure)
                matrix_group = f.create_group('matrix')
                
                # Store data, indices, indptr (CSR format)
                matrix_group.create_dataset('data', data=adata.X.data, compression='gzip')
                matrix_group.create_dataset('indices', data=adata.X.indices, compression='gzip')
                matrix_group.create_dataset('indptr', data=adata.X.indptr, compression='gzip')
                matrix_group.attrs['shape'] = adata.X.shape
                
                # Store gene names
                gene_group = matrix_group.create_group('features')
                gene_group.create_dataset('id', data=[g.encode() for g in gene_names])
                gene_group.create_dataset('name', data=[g.encode() for g in gene_names])
                gene_group.create_dataset('feature_type', data=[b'Gene Expression'] * n_genes)
                
                # Store barcodes (cell names)
                matrix_group.create_dataset('barcodes', data=[c.encode() for c in cell_names])
            
            print(f"   ✓ Created: {h5_path}")
            print(f"   - Shape: {counts.shape} (cells × genes)")
            print(f"   - Size: {h5_path.stat().st_size / 1024 / 1024:.2f} MB")
            print(f"   - Format: 10x Genomics H5 (simplified)")
            
        except Exception as e:
            print(f"   ✗ Error creating H5 file with scanpy: {e}")
            print(f"   Creating basic H5 format...")
            # Fallback: create a simple H5 file
            with h5py.File(h5_path, 'w') as f:
                f.create_dataset('counts', data=counts, compression='gzip')
                f.create_dataset('gene_names', data=[g.encode() for g in gene_names])
                f.create_dataset('cell_names', data=[c.encode() for c in cell_names])
            print(f"   ✓ Created basic H5: {h5_path}")
    else:
        # Create basic H5 format without scanpy
        print("   Creating basic H5 format (scanpy not available)...")
        with h5py.File(h5_path, 'w') as f:
            f.create_dataset('counts', data=counts, compression='gzip')
            f.create_dataset('gene_names', data=[g.encode() for g in gene_names])
            f.create_dataset('cell_names', data=[c.encode() for c in cell_names])
        print(f"   ✓ Created basic H5: {h5_path}")
        print(f"   - Shape: {counts.shape} (cells × genes)")
        print(f"   - Size: {h5_path.stat().st_size / 1024 / 1024:.2f} MB")
        print(f"   - Note: This is a simplified format. For full 10x compatibility, install scanpy.")
    
    print("\n" + "="*60)
    print("✓ All sample data files created successfully!")
    print("="*60)
    print(f"\nFiles created in: {output_dir.absolute()}")
    print(f"  - {csv_path.name} (CSV format)")
    print(f"  - {tsv_path.name} (TSV format)")
    print(f"  - {h5_path.name} (H5 format)")
    print(f"\n💡 You can now use these files in the web demo!")
    
    return {
        'csv': csv_path,
        'tsv': tsv_path,
        'h5': h5_path
    }


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Create sample data files in .h5, .csv, and .tsv formats')
    parser.add_argument('--output-dir', '-o', default='data',
                       help='Output directory for data files (default: data)')
    parser.add_argument('--n-cells', type=int, default=500,
                       help='Number of cells (default: 500)')
    parser.add_argument('--n-genes', type=int, default=2000,
                       help='Number of genes (default: 2000)')
    
    args = parser.parse_args()
    
    create_sample_data(
        output_dir=args.output_dir,
        n_cells=args.n_cells,
        n_genes=args.n_genes
    )

