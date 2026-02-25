# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Lotus is a Python SDK for single-cell RNA-seq analysis, built as a high-level wrapper around Scanpy with custom algorithms (CoreSpect/CoreMap from `cplearn`). It provides a unified API for the full scRNA-seq workflow: IO, preprocessing, clustering, annotation, and visualization.

## Build & Development

```bash
# Install with uv (project uses hatchling build system)
uv sync

# Install in editable mode
uv pip install -e .

# Run all tests
python -m pytest test/

# Run a single test
python -m pytest test/test_lotus/test_preprocess.py::test_run_preprocessing_mock_data -v
```

## Architecture

All source code lives under `src/lotus/`. The package is organized by analysis stage:

- **`io/`** - Data loading/writing. Wraps `scanpy.read_*` / `scanpy.write` with a unified `standardize_load()` that auto-detects format, ensures unique names, and converts to sparse.
- **`preprocessing/`** - QC, filtering, normalization, HVG selection, PCA, neighbors. `run_preprocessing()` is the all-in-one pipeline. All functions delegate to `scanpy.preprocessing`.
- **`clustering/`** - Leiden, Louvain, UMAP (via `scanpy.tl`), plus `corespect_clustering()` which wraps the custom CoreSpect algorithm.
- **`core_analysis/`** - The CoreSpect/CoreMap algorithms from `cplearn`. `_cplearn_core.py` is the entry point that configures `CorespectModel`, runs it, and writes results back to `adata.obs`/`adata.uns`.
- **`annotation_analysis/`** - DEG analysis (`rank_genes_groups`), pathway enrichment (GSEApy), cell type annotation (CellTypist).
- **`visualization/`** - Plot functions organized by category (embedding, PCA, preprocessing, trajectory, DE, generic). All wrap Scanpy plot functions.
- **`utils/metadata.py`** - `LoggingMeta` metaclass and `@logged` decorator for automatic method logging and `adata.X` snapshotting into layers.

### Key Patterns

- **Scanpy delegation**: Most functions are thin wrappers that pass through to Scanpy, preserving the exact same parameter signatures. This is intentional — Lotus adds orchestration, not reimplementation.
- **In-place by default**: Following Scanpy convention, most functions modify `adata` in-place and return `None`. Functions with `copy=True` return a new AnnData.
- **Module `__init__.py` as public API**: Each module's `__init__.py` explicitly imports and re-exports its public functions via `__all__`. This is the user-facing API surface.
- **CoreSpect results storage**: Cluster labels go to `adata.obs[key_added]`, core identity to `adata.obs[f"{key_added}_is_core"]`, and ragged layer indices are JSON-serialized into `adata.uns[key_added]['layers_indices_json']`.

## Language

The primary developer communicates in Chinese. Respond in the same language the user uses.
