# -*- coding: utf-8 -*-

from os.path import join, realpath, dirname
from os import getcwd
import traceback
# Import files with defined functions.
from functions.read_paths import read_paths as rp
from functions.get_in_params import get_in_params as gip
from functions.create_out_data_file import create_out_data_file as c_o_d_f
from functions.func_caller import asteca_funcs as af


def main():
    '''
    Main function. Reads input data files and calls the container function.
    '''

    __version__ = "v2.0.1-beta"

    print '-------------------------------------------'
    print '            ASteCA %s' % __version__
    print '-------------------------------------------\n'

    # Path where the code is running
    mypath = realpath(join(getcwd(), dirname(__file__)))

    # Read input parameters from params_input.dat file.
    gip_params = gip(mypath)
    # Read input/output paths.
    input_dir, output_dir = gip_params[1][:2]
    # Read paths and names of all clusters stored in input_dir.
    dir_files = rp(input_dir)
    # Create output data file in output_dir (append if file already exists)
    out_file_name = c_o_d_f(output_dir)

    # Iterate through all cluster files.
    for f_indx, sub_dir in enumerate(dir_files[0]):

        # Store name of file in 'myfile'.
        myfile = dir_files[1][f_indx]

        # Read input parameters from params_input.dat file again to reset params
        # lists.
        gip_params = gip(mypath)

        try:
            # Call function that calls sub-functions sequentially.
            af(myfile, sub_dir, out_file_name, gip_params)
        except Exception, err:
            print 'FATAL: %s could not be processed. ' % myfile[:-4]
            print 'Error:', str(err), '\n'
            print traceback.format_exc()

    print 'Full iteration completed.'

if __name__ == "__main__":
    main()