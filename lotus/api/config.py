"""
Configuration for Lotus API
"""

from pathlib import Path
import tempfile
import warnings

warnings.filterwarnings('ignore', category=UserWarning)

# Try to import lotus and scanpy
# Note: We're inside lotus package, so we need to import from parent or use absolute imports
import sys
from pathlib import Path

# Add parent directory to path to import lotus
_project_root = Path(__file__).parent.parent.parent.resolve()
if str(_project_root) not in sys.path:
    sys.path.insert(0, str(_project_root))

try:
    from lotus.workflows import preprocess, clustering, umap as lotus_umap
    from lotus import read, read_10x_h5, read_10x_mtx
    import lotus as lt
    LOTUS_AVAILABLE = True
    SCANPY_AVAILABLE = True
    import scanpy as sc
except ImportError as e:
    LOTUS_AVAILABLE = False
    preprocess = None
    clustering = None
    read = None
    read_10x_h5 = None
    read_10x_mtx = None
    lt = None
    try:
        import scanpy as sc
        SCANPY_AVAILABLE = True
    except ImportError:
        SCANPY_AVAILABLE = False
        sc = None
    print(f"Warning: Lotus not available ({e}). Some features will be disabled.")

# Upload folder configuration
UPLOAD_FOLDER = Path(tempfile.gettempdir()) / 'lotus_web_demo'
UPLOAD_FOLDER.mkdir(exist_ok=True, parents=True)

# Max file size: 500MB
MAX_CONTENT_LENGTH = 500 * 1024 * 1024

# Memory optimization settings for 4GB server
# Reduce default parameters to save memory
DEFAULT_N_PCS = 15  # Reduced from 20
DEFAULT_N_NEIGHBORS = 10  # Reduced from 15
DEFAULT_N_TOP_GENES = 1000  # Reduced from 2000

