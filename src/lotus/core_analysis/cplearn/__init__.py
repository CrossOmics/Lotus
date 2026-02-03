""" CPLEARN Library """
import warnings
from numba.core.errors import NumbaPendingDeprecationWarning


__version__ = "0.2.1"

# Silence specific warnings
warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message=r"Graph is not fully connected*",
    module=r"sklearn\.manifold\._spectral_embedding"
)

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=NumbaPendingDeprecationWarning)

__all__ = ["corespect"]