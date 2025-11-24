#!/bin/bash
# Run script for Lotus Web Demo

echo "=========================================="
echo "Lotus Embedding Projector Web Application"
echo "=========================================="
echo ""

# Check if virtual environment exists
if [ ! -d "venv" ]; then
    echo "Creating virtual environment..."
    python3 -m venv venv
fi

# Activate virtual environment
echo "Activating virtual environment..."
source venv/bin/activate

# Install/upgrade dependencies
echo "Installing dependencies..."
pip install --upgrade pip
pip install -r requirements.txt

# Check if lotus is available and install if needed
echo ""
echo "Checking Lotus availability..."
python3 -c "from lotus.workflows import preprocess, clustering; from lotus import read; print('✓ Lotus is available and working')" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "⚠ Lotus not found. Installing Lotus in development mode..."
    pip install -e .. --no-deps
    echo "Verifying installation..."
    python3 -c "from lotus.workflows import preprocess, clustering; from lotus import read; print('✓ Lotus installed successfully')" || echo "⚠ Warning: Lotus installation may have issues"
fi

# Run the application
echo ""
echo "Starting Flask server..."
echo "Open http://localhost:5000 in your browser"
echo "Press Ctrl+C to stop"
echo ""

python3 app.py

