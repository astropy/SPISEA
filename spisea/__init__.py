import importlib.metadata

try:
    __version__ = importlib.metadata.version(__package__)
except importlib.metadata.PackageNotFoundError:
    # Fallback for development mode if not installed
    __version__ = "2.4.1"
