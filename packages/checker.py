
from .check import first_run
from .check import pack


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

    # Check that all the essential packages are installed.
    inst_packgs_lst = pack.check()

    # Import here after the needed packages were checked to be present.
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
    from .check import read_met_files

    # Check if input cluster files exist.
    cl_files = clusters.check(mypath, file_end)

    # Read parameters from 'params_input.dat' file. Return a dictionary
    # containing all the parameter values.
    pd = params_file.check(mypath, file_end, inst_packgs_lst)

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
    params_struct.check(**pd)

    # Check decontamination algorithm parameters.
    params_decont.check(cl_files, **pd)

    # Check synthetic clusters parameters. Generate the
    # 'fundam_params' variable.
    pd = params_synthcl.check(mypath, pd)

    # Check the best match parameters.
    pd = params_match.check(pd)

    # Filters and colors names.
    fs = ', '.join(_[1] for _ in pd['filters'])
    cs = ', '.join('(' + _[1].replace(',', '-') + ')' for _ in pd['colors'])
    print("Filter: {}".format(fs))
    print("Color:  {}\n".format(cs))

    if pd['best_fit_algor'] != 'n':
        # Check and store metallicity files.
        pd = read_met_files.check_get(pd)

    return cl_files, pd
