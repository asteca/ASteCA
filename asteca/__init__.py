from .cluster import Cluster as cluster
from .membership import Membership as membership
from .isochrones import Isochrones as isochrones
from .synthetic import Synthetic as synthetic
from .likelihood import Likelihood as likelihood
from . import plot as plot

from contextlib import suppress
import importlib.metadata
from pathlib import Path


__all__ = ["cluster", "membership", "isochrones", "synthetic", "likelihood", "plot"]


def extract_version() -> str:
    """Returns either the version of the installed package or the one
    found in nearby pyproject.toml

    https://stackoverflow.com/a/76206192/1391441

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
    # return importlib.metadata.version(__package__
    #                                   or __name__.split(".", maxsplit=1)[0])
    return importlib.metadata.version(__package__ or __name__)


__version__ = extract_version()
