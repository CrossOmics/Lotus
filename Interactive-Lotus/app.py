"""
Flask backend for Lotus Web Application
Uses the modular API structure from lotus.api
"""

import sys
from pathlib import Path

# Add parent directory to path to import lotus
sys.path.insert(0, str(Path(__file__).parent.parent))

# Import the main function from lotus.api
from lotus.api.app import main

if __name__ == '__main__':
    main()
