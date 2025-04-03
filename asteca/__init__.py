import importlib.metadata
from contextlib import suppress
from pathlib import Path

from .cluster import Cluster
from .isochrones import Isochrones
from .likelihood import Likelihood
from .membership import Membership
from .synthetic import Synthetic

__all__ = [
    "Cluster",
    "Membership",
    "Isochrones",
    "Synthetic",
    "Likelihood",
]


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

# # Check if an updated version exists.
# # If this file exists, skip update check
# if os.path.isfile("asteca_disable_check.txt") is False:
#     # Get the latest version from PyPI
#     pypi_url = "https://pypi.org/pypi/asteca/json"
#     pypi_response = requests.get(pypi_url, timeout=3)
#     pypi_data = pypi_response.json()
#     new_v = pypi_data["info"]["version"]

#     # Parse and compare versions
#     if new_v != __version__:
#         print("\n--------------------------------------------------")
#         print(f"New version of ASteCA is available: {__version__} -> {new_v}")
#         print("   Update with: pip install --upgrade asteca\n")
#         print("   See what's new at: http://asteca.github.io/")
#         print("--------------------------------------------------\n")
