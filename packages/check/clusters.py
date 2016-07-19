
import sys


def check(cl_files, file_end):
    """
    Check the existence of cluster data files.
    """

    # Check that at least one photometric cluster file exists.
    if not cl_files:
        sys.exit("ERROR: no photometric data files found in '{}'\n"
                 "folder. Halting.".format('/input' + file_end))
