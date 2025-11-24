# Lotus Web Demo

Web application for Lotus Embedding Projector.

## Features

- Upload and visualize single-cell RNA-seq data
- Interactive embedding visualization
- Clustering and preprocessing workflows

## Quick Start

### Method 1: Using run.sh (Recommended)

```bash
cd web_demo
./run.sh
```

### Method 2: Manual Setup

```bash
cd web_demo

# Create virtual environment
python3 -m venv venv

# Activate virtual environment
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Install Lotus (in development mode from parent directory)
pip install -e .. --no-deps

# Run the application
python3 app.py
```

## Access the Application

Open your browser and navigate to: http://localhost:5000

## Restart Server

If the server is already running, you may need to stop it first:

```bash
# Find and stop the process
lsof -ti:5000 | xargs kill -9

# Then restart
cd web_demo
./run.sh
```

## Requirements

See `requirements.txt` for full list of dependencies. Main requirements:
- Flask >= 2.3.0
- numpy >= 1.24.0
- pandas >= 2.0.0
- anndata >= 0.9.0
- scanpy >= 1.9.0
- lotus (from parent directory)

## License

See parent repository for license information.

