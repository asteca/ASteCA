
import time
from os.path import join, realpath, dirname
from os import getcwd
import argparse
import traceback
from _version import __version__
from packages.first_run import check_1strun
from packages.inp import input_clusters
from packages.checker import check_all


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
            print("Will load parameters from 'params_input_{}.dat'"
                  " file.\n".format(file_end))
        else:
            print("Integer must be smaller than 99. Fall back to\ndefault"
                  " 'params_input.dat' file.\n")

    return file_end


def main():
    """
    Reads input data files and calls the container function.
    """
    # Start timing loop.
    start = time.time()

    print('\n-------------------------------------------')
    print('             [ASteCA {}]'.format(__version__))
    print('-------------------------------------------\n')

    # Root path where the code is running. Remove 'packages' from path.
    mypath = realpath(join(getcwd(), dirname(__file__)))[:-8]

    # Check .first_run file.
    check_1strun(mypath)

    # Read command-line argument.
    file_end = num_exec()

    # Read paths and names of all clusters stored inside /input.
    cl_files = input_clusters.main(mypath, file_end)

    # Checker function to verify that things are in place before running.
    # As part of the checking process, and to save time, the isochrone
    # files are read and stored here.
    # The 'R_in_place' flag indicates that R and rpy2 are installed.
    ip_list, R_in_place = check_all(mypath, file_end, cl_files)

    # Store those global variables that could be changed when processing each
    # cluster.
    # Import *after* checker function.
    import packages.inp.input_params as g
    global mode, er_params
    mode_orig, er_params_orig = g.mode, g.er_params

    # Import here to ensure the check has passed and all the necessary
    # packages are installed.
    from packages import func_caller

    # Iterate through all cluster files.
    for cl_file in cl_files:

        # Set these variables to their original global values.
        g.mode, g.er_params = mode_orig, er_params_orig

        try:
            # Call module that calls all sub-modules sequentially.
            func_caller.main(cl_file, ip_list, R_in_place)
        except Exception:
            print '\n!!! --> {}/{} could not be processed <-- !!!\n'.format(
                cl_file[-2], cl_file[-1])
            print traceback.format_exc()

    # End of run.
    elapsed = time.time() - start
    m, s = divmod(elapsed, 60)
    print 'Full run completed in {:.0f}m {:.0f}s.'.format(m, s)
