# Lotus

Lotus is a single-cell RNA-seq analysis pipeline built on top of Scanpy and AnnData.

## Installation

Install the Python dependencies with `uv`:

```bash
uv sync
```

## Lineage tracker prerequisite

The lineage tracker writes JSON data by itself, but PNG rendering also requires the Graphviz `dot` executable.

- `uv` / `pip` install the Python packages only
- they do **not** install the native Graphviz binary required by `diagrams`

Install Graphviz with your system package manager before expecting lineage PNG output:

### macOS

```bash
brew install graphviz
```

### Ubuntu / Debian

```bash
sudo apt-get install graphviz
```

### Windows

```powershell
winget install Graphviz.Graphviz
```

### Verify

```bash
dot -V
```

If Graphviz is missing, Lotus will still save `lineage.json` and will skip PNG rendering with a warning.