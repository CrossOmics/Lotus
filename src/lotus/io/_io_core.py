import warnings
from pathlib import Path
from typing import Literal, Sequence, Union, Mapping, Iterable, Iterator, Any

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from numba.cuda.cudadrv.nvvm import logger
from scipy.sparse import csr_matrix, csc_matrix, issparse


def standardize_load(
        file_path: Union[str, Path],
        *,
        source_type: Literal[
            'auto', '10x_mtx', '10x_h5', 'h5ad', 'csv', 'tsv', 'txt',
            'loom', 'excel', 'visium', 'hdf', 'mtx'
        ] = 'auto',
        backed: Literal['r', 'r+'] | None = None,
        make_unique: bool = True,
        convert_to_sparse: bool = True,
        var_names: Literal['gene_symbols', 'gene_ids'] = 'gene_symbols',
        **kwargs
) -> AnnData:
    """
    Load single-cell data from various formats into a standardized AnnData object.

    - Auto-detects format if source_type='auto'
    - Ensures unique obs/var names
    - Optionally converts dense matrices to sparse
    - Performs basic integrity checks
    """

    # Path validation
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Path does not exist: {path}")

    # Format detection
    if source_type == 'auto':
        if path.is_dir():
            if (
                    (path / 'matrix.mtx').exists()
                    or (path / 'matrix.mtx.gz').exists()
                    or (path / 'features.tsv').exists()
                    or (path / 'genes.tsv').exists()
            ):
                source_type = '10x_mtx'
            elif (path / 'spatial').exists():
                source_type = 'visium'
            else:
                raise ValueError(f"Unrecognized directory format: {path}")
        else:
            suffix = path.suffix.lower()
            suffix_map = {
                '.h5ad': 'h5ad',
                '.h5': '10x_h5',
                '.csv': 'csv',
                '.tsv': 'tsv',
                '.txt': 'txt',
                '.loom': 'loom',
                '.xlsx': 'excel',
                '.xls': 'excel',
                '.mtx': 'mtx',
                '.hdf': 'hdf',
                '.hdf5': 'hdf',
            }
            if suffix in suffix_map:
                source_type = suffix_map[suffix]
            elif path.suffixes[-2:] == ['.mtx', '.gz']:
                source_type = 'mtx'
            else:
                raise ValueError(f"Cannot auto-detect format: {path}")

    # Load data
    try:
        if source_type == '10x_mtx':
            adata = sc.read_10x_mtx(
                path,
                var_names=var_names,
                make_unique=kwargs.get('make_unique', True),
                cache=kwargs.get('cache', False),
                gex_only=kwargs.get('gex_only', True),
            )

        elif source_type == '10x_h5':
            adata = sc.read_10x_h5(
                path,
                genome=kwargs.get('genome'),
                gex_only=kwargs.get('gex_only', True),
            )

        elif source_type == 'h5ad':
            adata = sc.read_h5ad(
                path,
                backed=backed,
                as_sparse=kwargs.get('as_sparse', ()),
                chunk_size=kwargs.get('chunk_size', 6000),
            )

        elif source_type in {'csv', 'tsv'}:
            adata = sc.read_csv(
                path,
                delimiter=',' if source_type == 'csv' else '\t',
                first_column_names=kwargs.get('first_column_names'),
                dtype=kwargs.get('dtype', 'float32'),
            )

        elif source_type == 'txt':
            adata = sc.read_text(
                path,
                delimiter=kwargs.get('delimiter'),
                first_column_names=kwargs.get('first_column_names'),
                dtype=kwargs.get('dtype', 'float32'),
            )

        elif source_type == 'loom':
            adata = sc.read_loom(
                path,
                sparse=kwargs.get('sparse', True),
                cleanup=kwargs.get('cleanup', False),
                X_name=kwargs.get('X_name', 'spliced'),
                obs_names=kwargs.get('obs_names', 'CellID'),
                var_names=kwargs.get('var_names_loom', 'Gene'),
                dtype=kwargs.get('dtype', 'float32'),
            )

        elif source_type == 'excel':
            adata = read_excel(
                path,
                sheet=kwargs.get('sheet', 0),
                dtype=kwargs.get('dtype', 'float32'),
            )


        elif source_type == 'visium':
            adata = sc.read_visium(
                path,
                genome=kwargs.get('genome'),
                count_file=kwargs.get('count_file', 'filtered_feature_bc_matrix.h5'),
                library_id=kwargs.get('library_id'),
                load_images=kwargs.get('load_images', True),
            )

        elif source_type == 'mtx':
            adata = sc.read_mtx(path, dtype=kwargs.get('dtype', 'float32'))

        elif source_type == 'hdf':
            adata = sc.read_hdf(path, key=kwargs.get('key', 'data'))

        else:
            raise ValueError(f"Unsupported source_type: {source_type}")

    except Exception as e:
        raise RuntimeError(f"Failed to load {path} ({source_type}): {e}")

    # Standardization
    if make_unique:
        if not adata.var_names.is_unique:
            warnings.warn("Non-unique var names; making unique.")
            adata.var_names_make_unique()
        if not adata.obs_names.is_unique:
            warnings.warn("Non-unique obs names; making unique.")
            adata.obs_names_make_unique()

    if convert_to_sparse and adata.X is not None and not issparse(adata.X):
        sparsity = 1 - (adata.X != 0).sum() / adata.X.size
        if sparsity > 0.5:
            adata.X = csr_matrix(adata.X)

    # Validation
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError("AnnData contains no cells or no genes.")

    if adata.X is not None and adata.X.shape != (adata.n_obs, adata.n_vars):
        raise ValueError("X shape does not match obs/var dimensions.")

    if adata.obs_names.isna().any() or adata.var_names.isna().any():
        raise ValueError("Found NaN values in obs or var names.")

    return adata


def write(
        filename: str | Path,
        adata: AnnData,
        *,
        ext: Literal['h5', 'csv', 'txt', 'npz'] | None = None,
        compression: Literal['gzip', 'lzf'] | None = 'gzip',
        compression_opts: int | None = None,
) -> None:
    """
    Write AnnData objects to file.

    Parameters:
        filename: If the filename has no file extension, it is interpreted as a key for
                  generating a filename via `sc.settings.writedir`.
        adata: Annotated data matrix.
        ext: File extension from which to infer file format. If None, defaults to
             `sc.settings.file_format_data`.
        compression: Compression strategy. See h5py documentation.
        compression_opts: Compression options. See h5py documentation.

    Returns:
        None.
    """
    sc.write(
        filename,
        adata,
        ext=ext,
        compression=compression,
        compression_opts=compression_opts
    )


def read(
        filename: str | Path,
        backed: Literal['r', 'r+'] | None = None,
        *,
        sheet: str | None = None,
        ext: str | None = None,
        delimiter: str | None = None,
        first_column_names: bool = False,
        backup_url: str | None = None,
        cache: bool = False,
        cache_compression: Literal['gzip', 'lzf'] | None = None,
        **kwargs,
) -> AnnData:
    """
    Read file and return AnnData object.

    This is a generic read function that infers the file format from the extension.

    Parameters:
        filename: File name or path.
        backed: If 'r', load AnnData in backed mode (memory efficient).
        sheet: Name of sheet/table in hdf5 or Excel file.
        ext: Extension that indicates the file type.
        delimiter: Delimiter that separates data within text file.
        first_column_names: Assume the first column stores row names.
        backup_url: Retrieve the file from an URL if not present on disk.
        cache: If True, read from fast 'h5ad' cache.
        cache_compression: Compression for the cache file.
        **kwargs: Parameters passed to read_loom().

    Returns:
        AnnData object.
    """
    return sc.read(
        filename,
        backed=backed,
        sheet=sheet,
        ext=ext,
        delimiter=delimiter,
        first_column_names=first_column_names,
        backup_url=backup_url,
        cache=cache,
        cache_compression=cache_compression,
        **kwargs
    )


def read_10x_h5(
        filename: str | Path,
        *,
        genome: str | None = None,
        gex_only: bool = True,
        backup_url: str | None = None,
) -> AnnData:
    """
    Read 10x-Genomics-formatted hdf5 file.

    Parameters:
        filename: Path to a 10x hdf5 file.
        genome: Filter expression to genes within this genome.
        gex_only: Only keep 'Gene Expression' data and ignore other feature types.
        backup_url: Retrieve the file from an URL if not present on disk.

    Returns:
        AnnData object.
    """
    return sc.read_10x_h5(
        filename,
        genome=genome,
        gex_only=gex_only,
        backup_url=backup_url
    )


def read_10x_mtx(
        path: str | Path,
        *,
        var_names: Literal['gene_symbols', 'gene_ids'] = 'gene_symbols',
        make_unique: bool = True,
        cache: bool = False,
        cache_compression: Literal['gzip', 'lzf'] | None = None,
        gex_only: bool = True,
        prefix: str | None = None,
) -> AnnData:
    """
    Read 10x-Genomics-formatted mtx directory.

    Parameters:
        path: Path to directory containing .mtx and .tsv files.
        var_names: The variables index ('gene_symbols' or 'gene_ids').
        make_unique: Whether to make the variables index unique.
        cache: If True, read from fast 'h5ad' cache.
        cache_compression: Compression for the cache file.
        gex_only: Only keep 'Gene Expression' data.
        prefix: Any prefix before matrix.mtx, genes.tsv and barcodes.tsv.

    Returns:
        AnnData object.
    """
    return sc.read_10x_mtx(
        path,
        var_names=var_names,
        make_unique=make_unique,
        cache=cache,
        cache_compression=cache_compression,
        gex_only=gex_only,
        prefix=prefix
    )


def read_visium(
        path: str | Path,
        genome: str | None = None,
        *,
        count_file: str = 'filtered_feature_bc_matrix.h5',
        library_id: str | None = None,
        load_images: bool = True,
        source_image_path: str | Path | None = None,
) -> AnnData:
    """
    Read 10x-Genomics-formatted Visium dataset.

    Note: As of Scanpy 1.11.0, use of squidpy.read.visium() is recommended.

    Parameters:
        path: Path to directory for visium datafiles.
        genome: Filter expression to genes within this genome.
        count_file: Which file in the passed directory to use as the count file.
        library_id: Identifier for the visium library.
        load_images: Whether to load images.
        source_image_path: Path to the high-resolution tissue image.

    Returns:
        AnnData object containing spatial data.
    """
    return sc.read_visium(
        path,
        genome=genome,
        count_file=count_file,
        library_id=library_id,
        load_images=load_images,
        source_image_path=source_image_path
    )


def read_h5ad(
        filename: str | Path,
        backed: Literal['r', 'r+'] | bool | None = None,
        *,
        as_sparse: Sequence[str] = (),
        as_sparse_fmt: type[csr_matrix | csc_matrix] = csr_matrix,
        chunk_size: int = 6000,
) -> AnnData:
    """
    Read .h5ad-formatted hdf5 file.

    Parameters:
        filename: File name of data file.
        backed: If 'r', load AnnData in backed mode.
        as_sparse: If an array was saved as dense, read it as sparse.
        as_sparse_fmt: Sparse format class to read elements from as_sparse in as.
        chunk_size: Chunk size for loading sparse dataset stored as dense.

    Returns:
        AnnData object.
    """
    return sc.read_h5ad(
        filename,
        backed=backed,
        as_sparse=as_sparse,
        as_sparse_fmt=as_sparse_fmt,
        chunk_size=chunk_size
    )


def read_csv(
        filename: str | Path | Iterator[str],
        delimiter: str | None = ',',
        *,
        first_column_names: bool | None = None,
        dtype: str = 'float32',
) -> AnnData:
    """
    Read .csv file.

    Parameters:
        filename: Data file.
        delimiter: Delimiter that separates data.
        first_column_names: Assume the first column stores row names.
        dtype: Numpy data type.

    Returns:
        AnnData object.
    """
    return sc.read_csv(
        filename,
        delimiter=delimiter,
        first_column_names=first_column_names,
        dtype=dtype
    )


def read_excel(
        filename: str | Path,
        sheet: str | int,
        dtype: str = 'float32',
) -> AnnData:
    """
    Read .xlsx (Excel) file.

    Parameters:
        filename: File name to read from.
        sheet: Name of sheet in Excel file.
        dtype: Numpy data type.

    Returns:
        AnnData object.
    """
    converted_data = sc.read_excel(
        filename,
        sheet,
        dtype=dtype
    )
    # fix the bug of original scanpy
    # cast Excel data to float32 to prevent object dtype errors
    converted_data.X = converted_data.X.astype('float32')
    return converted_data


def read_hdf(
        filename: str | Path,
        key: str,
) -> AnnData:
    """
    Read .h5 (hdf5) file.

    Parameters:
        filename: Filename of data file.
        key: Name of dataset in the file.

    Returns:
        AnnData object.
    """
    return sc.read_hdf(
        filename,
        key
    )


def read_loom(
        filename: str | Path,
        *,
        sparse: bool = True,
        cleanup: bool = False,
        X_name: str = 'spliced',
        obs_names: str = 'CellID',
        obsm_names: str | None = None,
        var_names: str = 'Gene',
        varm_names: str | None = None,
        dtype: str = 'float32',
        obsm_mapping: Mapping[str, Iterable[str]] = dict(),
        varm_mapping: Mapping[str, Iterable[str]] = dict(),
        **kwargs,
) -> AnnData:
    """
    Read .loom-formatted hdf5 file.

    Parameters:
        filename: The filename.
        sparse: Whether to read the data matrix as sparse.
        cleanup: Whether to collapse unique obs/var fields.
        X_name: Loompy key for data matrix X.
        obs_names: Loompy key for observation/cell names.
        obsm_names: Deprecated/Loompy specific.
        var_names: Loompy key for variable/gene names.
        varm_names: Deprecated/Loompy specific.
        dtype: Data type.
        obsm_mapping: Mapping for observation matrices.
        varm_mapping: Mapping for variable matrices.
        **kwargs: Arguments to loompy.connect.

    Returns:
        AnnData object.
    """
    return sc.read_loom(
        filename,
        sparse=sparse,
        cleanup=cleanup,
        X_name=X_name,
        obs_names=obs_names,
        obsm_names=obsm_names,
        var_names=var_names,
        varm_names=varm_names,
        dtype=dtype,
        obsm_mapping=obsm_mapping,
        varm_mapping=varm_mapping,
        **kwargs
    )


def read_mtx(
        filename: str | Path,
        dtype: str = 'float32',
) -> AnnData:
    """
    Read .mtx file.

    Parameters:
        filename: The filename.
        dtype: Numpy data type.

    Returns:
        AnnData object.
    """
    return sc.read_mtx(
        filename,
        dtype=dtype
    )


def read_text(
        filename: str | Path | Iterator[str],
        delimiter: str | None = None,
        *,
        first_column_names: bool | None = None,
        dtype: str = 'float32',
) -> AnnData:
    """
    Read .txt, .tab, .data (text) file.

    Parameters:
        filename: Data file, filename or stream.
        delimiter: Delimiter that separates data.
        first_column_names: Assume the first column stores row names.
        dtype: Numpy data type.

    Returns:
        AnnData object.
    """
    return sc.read_text(
        filename,
        delimiter=delimiter,
        first_column_names=first_column_names,
        dtype=dtype
    )


def read_umi_tools(
        filename: str | Path,
        dtype: str | None = None,
) -> AnnData:
    """
    Read a gzipped condensed count matrix from umi_tools.

    Parameters:
        filename: File name to read from.
        dtype: Data type.

    Returns:
        AnnData object.
    """
    return sc.read_umi_tools(
        filename,
        dtype=dtype
    )
