# -*- coding: utf-8 -*-

import time
from os.path import join, realpath, dirname
from os import getcwd
import argparse
import traceback
from functions import __version__
from functions.checker import check
from functions._in.get_in_clusters import in_clusters


def num_exec():
    """
    Parse optional command-line argument. The integer (if) passed, will be used
    as the two last characters in the file name of a
    'params_input_XX.dat' file.
    Integer must be smaller than 99, else default 'params_input.dat' file is
    used.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-N", help="Integer. Set the params_input_N.dat"
                        " file to be used by this run.", type=int)
    args = parser.parse_args()
    file_end = ''
    if args.N:
        if args.N < 99:
            file_end = "{:0>2}".format(args.N)
            print("Will load parameters from params_input_{}.dat\n".format(
                file_end))
        else:
            print("Integer must be smaller than 99. Fall back to default"
                  " params_input.dat file.\n")

    return file_end


def main():
    """
    Main function. Reads input data files and calls the container function.
    """
    # Start timing loop.
    start = time.time()

    print('\n')
    print('-------------------------------------------')
    print('             [ASteCA {}]'.format(__version__))
    print('-------------------------------------------\n')

    file_end = num_exec()

    # Root path where the code is running. Remove 'functions' from path.
    mypath = realpath(join(getcwd(), dirname(__file__)))[:-9]

    # Read paths and names of all clusters stored inside /input.
    cl_files = in_clusters(mypath, file_end)

    # Checker function to verify that things are in place before running.
    # As part of the checking process, and to save time, the isochrone
    # files are read and stored here.
    # The 'R_in_place' flag indicates that R and rpy2 are installed.
    ip_list, R_in_place = check(mypath, file_end, cl_files)

    # Store those global variables that could be changed when processing each
    # cluster.
    # Import *after* checker function.
    import functions._in.get_in_params as g
    global mode, er_params
    mode_orig, er_params_orig = g.mode, g.er_params

    # Import here to ensure the check has passed and all the necessary
    # packages are installed.
    from functions.func_caller import asteca_funcs as af

    # Iterate through all cluster files.
    for cl_file in cl_files:

        # Set these variables to their original global values.
        g.mode, g.er_params = mode_orig, er_params_orig

        try:
            # Call function that calls all sub-functions sequentially.
            af(cl_file, ip_list, R_in_place)
        except Exception:
            print '\n!!! --> {}/{} could not be processed <-- !!!\n'.format(
                cl_file[-2], cl_file[-1])
            print traceback.format_exc()

    # End of run.
    elapsed = time.time() - start
    m, s = divmod(elapsed, 60)
    print 'Full run completed in {:.0f}m {:.0f}s.'.format(m, s)


if __name__ == "__main__":
    main()
