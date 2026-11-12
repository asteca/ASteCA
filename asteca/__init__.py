import importlib.metadata
from contextlib import suppress
from pathlib import Path

from .cluster import Cluster
from .isochrones import Isochrones
from .membership import Membership
from .sampler import Sampler
from .synthetic import Synthetic

# The order of this list sets the order of the classes in the API documentation, so
# it is not sorted alphabetically
__all__ = [  # noqa: RUF022
    "Cluster",
    "Isochrones",
    "Membership",
    "Synthetic",
    "Sampler",
]


def extract_version() -> str:
    """Returns either the version of the installed package or the one
    found in nearby pyproject.toml
    """
    with suppress(FileNotFoundError, StopIteration):
        pp_file = Path(__file__).parent.parent / "pyproject.toml"
        with open(pp_file, encoding="utf-8") as pyproject_toml:
            version = (
                next(line for line in pyproject_toml if line.startswith("version"))
                .split("=")[1]
                .strip("'\"\n ")
            )
            return f"{version}"
    return importlib.metadata.version(__package__ or __name__)


__version__ = extract_version()
