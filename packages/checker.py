
from .check import first_run
from .check import clusters
from .check import params_file
from .check import update
from .check import params_data
from .check import params_kinem
from .check import params_out
from .check import params_struct
from .check import params_decont
from .check import params_synthcl
from .check import params_match


def check_all(mypath, file_end):
    """
    Checks that the necessary files are in place, the parameters stored
    in the input file are valid, and that the ranges given for the cluster
    parameters are consistent with the isochrones available before moving
    on with the code.
    """

    # Check .first_run file.
    first_run.main(mypath)

    print('Checking input parameters...\n')

    # # DEPRECATED 05/2021
    # # Check that all the essential packages are installed.
    # inst_packgs_lst = pack.check()

    # Check if input cluster files exist.
    cl_files = clusters.check(mypath, file_end)

    # Read parameters from 'asteca.ini' file. Return a dictionary
    # containing all the parameter values.
    pd = params_file.check(mypath, file_end)

    # Check if a new version is available.
    update.check()

    # Check that the data column indexes/names were properly given, and that
    # the magnitude and color names were properly defined.
    # If they are, store also the name of the proper isochrones folders.
    pd = params_data.check(mypath, pd)

    params_kinem.check(pd)

    # Check output parameters.
    pd = params_out.check(pd)

    # Check structural parameters.
    pd = params_struct.check(pd)

    # Check decontamination algorithm parameters.
    params_decont.check(cl_files, **pd)

    # Check synthetic clusters parameters. Generate the
    # 'fundam_params_all' variable.
    pd = params_synthcl.check(mypath, pd)

    # Check the best match parameters.
    pd = params_match.check(pd)

    # Filters and colors names.
    fs = ', '.join(_[1] for _ in pd['filters'])
    cs = ', '.join('(' + _[1].replace(',', '-') + ')' for _ in pd['colors'])
    print("Filter: {}".format(fs))
    print("Color:  {}\n".format(cs))

    return cl_files, pd
