

# Lotus

A single-cell RNA sequencing analysis package for preprocessing, clustering, visualization, and differential expression analysis.

## Features

- **Preprocessing**: Quality control, normalization, and feature selection
- **Clustering**: Cell type identification and clustering workflows
- **Visualization**: UMAP embedding and interactive visualizations
- **Differential Expression**: Marker gene identification and DEG analysis
- **Core Analysis**: Integrated analysis workflows

## Installation

Install the lotus package in development mode:

```bash
pip install -e .
```

This will make the `lotus` package available system-wide without needing to set `PYTHONPATH`.

## Quick Start

See `examples/lotus_workflow.py` for a complete workflow example.

## Lotus Embedding Projector

We developed an interactive analysis tool, available at `Interactive-Lotus/`. See the [Interactive-Lotus README](Interactive-Lotus/README.md) for setup instructions. We also support a [lightweight web demo](https://lotus.zeabur.app/) for demonstration. 

**Note:** The web demo has computational resource limitations and only supports very small datasets.

## Documentation

Full documentation is available at https://crossomics.github.io/Lotus/.
