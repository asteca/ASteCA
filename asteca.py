# -*- coding: utf-8 -*-

from os.path import join, realpath, dirname
from os import getcwd
import traceback
from functions import __version__
from functions.checker import check
from functions.get_in_clusters import in_clusters


def main():
    '''
    Main function. Reads input data files and calls the container function.
    '''
    print '\n'
    print '-------------------------------------------'
    print '            ASteCA %s' % __version__
    print '-------------------------------------------\n'

    # Get path where the code is running
    mypath = realpath(join(getcwd(), dirname(__file__)))

    # Read paths and names of all clusters stored inside /input.
    cl_files = in_clusters(mypath)

    # Checker function to verify that things are in place before running.
    # As part of the checking process, and to save time, the isochrone
    # files are read and stored here.
    # The 'R_in_place' flag indicates that R and rpy2 are installed.
    ip_list, R_in_place = check(mypath, cl_files)

    # Import here to ensure the check has passed and all the necessary
    # packages are installed.
    from functions.func_caller import asteca_funcs as af

    # Iterate through all cluster files.
    for cl_file in cl_files:

        try:
            # Call function that calls all sub-functions sequentially.
            af(mypath, cl_file, ip_list, R_in_place)
        except Exception:
            print '\n!!! --> {}/{} could not be processed <-- !!!\n'.format(
                *cl_file)
            print traceback.format_exc()

    print 'Full iteration completed.'


if __name__ == "__main__":
    main()