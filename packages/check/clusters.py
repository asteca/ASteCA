
from packages.inp import input_clusters


def check(mypath, file_end):
    """
    Check the existence of cluster data files.
    """

    # Read paths and names of all clusters stored inside /input.
    cl_files = input_clusters.main(mypath, file_end)

    # Check that at least one photometric cluster file exists.
    if not cl_files:
        raise ValueError(
            "no photometric data files found in '{}'\n"
            "folder. Halting.".format('/input' + file_end))
    else:
        return cl_files
