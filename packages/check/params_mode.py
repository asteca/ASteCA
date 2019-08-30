
import sys
from os import extsep


def check(mypath, cl_files, run_mode, semi_input, **kwargs):
    """
    Check that running mode parameters are properly written.
    """

    # Check mode.
    if run_mode not in ('auto', 'semi'):
        sys.exit("ERROR: 'mode' value selected ('{}') is not valid.".format(
            run_mode))

    if run_mode == 'semi':
        for clust_path in cl_files:
            cl_split = clust_path[-1].split(extsep)
            clust_name = '.'.join(cl_split[:-1])

            # Flag to indicate if cluster was found in file.
            flag_clust_found = 0
            for line in semi_input:
                if line[0] == clust_name:
                    # Set flag to True if the cluster was found.
                    flag_clust_found += 1

            # Cluster not found.
            if flag_clust_found == 0:
                sys.exit("ERROR: 'semi' mode is set but '{}' was not found\n"
                         "in the 'Semi Input' block".format(clust_name))

            # Cluster is duplicated.
            if flag_clust_found > 1:
                sys.exit("ERROR: 'semi' mode is set and '{}' is duplicated\n"
                         "in the 'Semi Input' block".format(clust_name))
