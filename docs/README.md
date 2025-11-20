# Lotus Documentation Build Instructions

This document explains how to build the Lotus API reference documentation.

## Prerequisites

First, install Sphinx and related dependencies:

```bash
pip install sphinx furo
```

Or install all dependencies from the project root (including documentation build dependencies):

```bash
pip install -e .
```

## Building the Documentation

### Using Makefile (Recommended for Linux/macOS)

```bash
cd docs
make html
```

The generated HTML documentation will be in the `docs/_build/html/` directory.

### Using sphinx-build (All platforms)

```bash
cd docs
sphinx-build -b html . _build/html
```

### Windows

```bash
cd docs
make.bat html
```

## Viewing the Documentation

After building, open `docs/_build/html/index.html` in your browser to view the documentation.

## Cleaning Build Files

```bash
cd docs
make clean
```

Or

```bash
cd docs
rm -rf _build
```

## Other Build Options

- `make html` - Build HTML documentation
- `make clean` - Clean the build directory
- `make help` - View all available options

## Documentation Structure

- `index.rst` - Documentation homepage
- `getting_started.rst` - Quick start guide
- `user_guide/` - User guide
- `api/` - API reference documentation

## Development

The documentation uses Sphinx to automatically extract docstrings from code. Main configuration is in `conf.py`.

To modify the documentation:
1. Modify the corresponding `.rst` files
2. Re-run `make html`
3. View results in `_build/html/`

## Notes

- If you see "No module named 'scanpy'" warnings during build, this is normal. These warnings occur because Sphinx cannot find scanpy when importing modules, but they do not affect documentation generation. In actual usage environments, lotus will correctly install scanpy as a dependency.
- The documentation automatically extracts docstrings from code, so after modifying docstrings in code, rebuild to update the documentation.
