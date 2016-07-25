
import sys
from os.path import join, isfile
import traceback


def check(mypath, file_end):
    """
    Check the existence of the 'params_input.dat' file initialize
    global variables, and check that the parameters are properly written.
    """

    # Check if params_input_XX.dat file exists.
    pars_f_name = 'params_input' + file_end + '.dat'
    pars_f_path = join(mypath, pars_f_name)
    if file_end == '':
        if not isfile(pars_f_path):
            sys.exit("ERROR: '{}' file does not exist.".format(pars_f_name))
    else:
        if not isfile(pars_f_path):
            print ("  WARNING: {} file does not exist.\n  Falling back to"
                   " 'params_input.dat' file.\n".format(pars_f_name))

            # Fall back to default file.
            pars_f_name = 'params_input.dat'
            pars_f_path = join(mypath, pars_f_name)
            if not isfile(pars_f_path):
                sys.exit("ERROR: '{}' file does not exist.".format(
                    pars_f_name))

    # Module with all input parameters loaded as global variables.
    import packages.inp.input_params as g
    # Check if params_input_XX.dat file is properly formatted.
    try:
        # Read input parameters from params_input.dat file. Initialize
        # global variables.
        g.init(mypath, pars_f_path)
    except Exception:
        # Halt code.
        print traceback.format_exc()
        sys.exit("ERROR: '{}' is badly formatted.".format(pars_f_name))

    # Return 'check for available update', and 'force backend' flags.
    return g.up_flag, g.flag_back_force
