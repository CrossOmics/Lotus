# Lotus Web Demo

Web application for Lotus Embedding Projector and provide a [lightweight web demo](https://huggingface.co/spaces/zzq1zh/Lotus-hf). Due to computational resource limitations, this web demo only supports small datasets. We encourage you to deploy locally using the following method.


## Features

- Upload and visualize single-cell RNA-seq data
- Interactive embedding visualization
- Clustering and preprocessing workflows

## Quick Start

```bash
cd Lotus-Web-Demo

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
cd Lotus-Web-Demo
source venv/bin/activate
python3 app.py

