
import sys
from os import extsep


def check(cl_files, run_mode, **kwargs):
    """
    Check that the parameters are properly written.
    """
    semi_file = 'semi_input.dat'  # HARDCODED
    # Mode is semi.
    if run_mode == 'semi':
        for clust_path in cl_files:
            clust_name = clust_path[-1].split(extsep)[0]
            # Flag to indicate if cluster was found in file.
            flag_clust_found = False
            with open(semi_file, "r") as f_cl_dt:
                for line in f_cl_dt:
                    li = line.strip()
                    # Skip comments.
                    if not li.startswith("#"):
                        reader = li.split()
                        # Prevent empty lines with spaces detected as a cluster
                        # line from crashing the code.
                        if reader:
                            # If cluster is found in file.
                            if reader[0] == clust_name:
                                # Set flag to True if the cluster was found.
                                flag_clust_found = True

            # Cluster not found.
            if not flag_clust_found:
                # Name of cluster not found in semi_input file.
                sys.exit("ERROR: 'semi' mode is set but '{}' was not found\n"
                         "in 'semi_input.dat' file".format(clust_name))
