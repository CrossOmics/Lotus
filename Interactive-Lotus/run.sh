#!/bin/bash
# Run script for Lotus Web Demo

echo "=========================================="
echo "Lotus Web Application"
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
python3 -c "from lotus.workflows import preprocess, cluster; from lotus import read; print('✓ Lotus is available and working')" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "⚠ Lotus not found. Installing Lotus in development mode..."
    pip install -e .. --no-deps
    echo "Verifying installation..."
    python3 -c "from lotus.workflows import preprocess, cluster; from lotus import read; print('✓ Lotus installed successfully')" || echo "⚠ Warning: Lotus installation may have issues"
fi

# Run the application
echo ""
echo "Starting Flask server..."
<<<<<<< HEAD
echo "Open http://localhost:5000 in your browser"
echo "Press Ctrl+C to stop"
echo ""

python3 app.py
=======
# Default port is 5259 (5000 is often used by AirPlay on macOS)
PORT=${PORT:-5259}
echo "Open http://localhost:${PORT} in your browser"
echo "Press Ctrl+C to stop"
echo ""

PORT=${PORT} python3 app.py
>>>>>>> 8b02f6559a88ee0942e4fa185dfd86d6e1447c40

