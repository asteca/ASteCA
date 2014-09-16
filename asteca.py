# -*- coding: utf-8 -*-

from os.path import join, realpath, dirname
from os import getcwd
import traceback
# Import files with defined functions.
from functions.get_in_clusters import in_clusters
from functions.func_caller import asteca_funcs as af


def main():
    '''
    Main function. Reads input data files and calls the container function.
    '''

    __version__ = "v3.0.0-beta"

    print '-------------------------------------------'
    print '            ASteCA %s' % __version__
    print '-------------------------------------------\n'

    # Path where the code is running
    mypath = realpath(join(getcwd(), dirname(__file__)))
    # Read paths and names of all clusters stored in /input.
    cl_files = in_clusters(mypath)

    # Iterate through all cluster files.
    for cl_file in cl_files:

        try:
            # Call function that calls all sub-functions sequentially.
            af(mypath, cl_file)
        except Exception, err:
            print 'FATAL: {}/{} could not be processed.'.format(*cl_file)
            print 'Error:', str(err), '\n'
            print traceback.format_exc()

    print 'Full iteration completed.'

if __name__ == "__main__":
    main()