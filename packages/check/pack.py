

import sys
import imp


def check():
    """
    Check if all the essential (and Non-essential) packages are installed.

    Source: https://stackoverflow.com/a/14050282/1391441
    """
    # Required packages.
    missing_pckg, inst_packgs_lst = [], []
    for pckg in ['numpy', 'matplotlib', 'scipy', 'astropy']:
        try:
            imp.find_module(pckg)
            inst_packgs_lst.append(pckg)
        except ImportError:
            missing_pckg.append(pckg)

    if missing_pckg:
        print("ERROR: the following required packages are missing:\n")
        for p in missing_pckg:
            print(" - {}".format(p))
        sys.exit("Missing packages")

    # Non-essential packages
    for pckg in ['rpy2']:
        try:
            imp.find_module(pckg)
            inst_packgs_lst.append(pckg)
        except ImportError:
            pass

    return inst_packgs_lst


if __name__ == "__main__":
    check()
