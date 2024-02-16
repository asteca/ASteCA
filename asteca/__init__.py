import importlib.metadata
from .asteca import isochrones
from .asteca import synth_clusters
from .asteca import cluster
from .asteca import likelihood

# Pull version from pyproject.toml
__version__ = importlib.metadata.version(__package__ or __name__)
