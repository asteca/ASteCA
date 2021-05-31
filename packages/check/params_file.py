
from os.path import join, isfile
import traceback
from packages.inp import input_params


def check(mypath, file_end):
    """
    Check the existence of the 'asteca.ini' file, and check that
    the parameters are properly written.
    """

    # Check if 'asteca_XX.ini' file exists.
    pars_f_name = 'asteca' + file_end + '.ini'
    pars_f_path = join(mypath, pars_f_name)
    if file_end == '':
        if not isfile(pars_f_path):
            raise ValueError("'{}' file does not exist.".format(pars_f_name))
    else:
        if not isfile(pars_f_path):
            print("  WARNING: {} file does not exist.\n  Falling back to"
                  " 'asteca.ini' file.\n".format(pars_f_name))

            # Fall back to default file.
            pars_f_name = 'asteca.ini'
            pars_f_path = join(mypath, pars_f_name)
            if not isfile(pars_f_path):
                raise ValueError("'{}' file does not exist.".format(
                    pars_f_name))
    # Check if asteca_XX.ini file is properly formatted.
    try:
        # Read input parameters from asteca.ini file.
        pd = input_params.main(pars_f_path)
    except Exception:
        # Halt code.
        print(traceback.format_exc())
        raise ValueError("'{}' is badly formatted.".format(pars_f_name))

    # Add to parameters dictionary.
    pd['file_end'] = file_end

    if pd['best_fit_algor'] not in pd['optimz_algors']:
        raise ValueError("the selected best fit method '{}' does not match"
                         " a valid input.".format(pd['best_fit_algor']))

    return pd
