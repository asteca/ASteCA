
from .check import pack
from .check import first_run
from .check import update
from .check import clusters
from .check import params_file
from .check import params_mode
from .check import params_data
from .check import params_out
from .check import params_struct
from .check import params_decont


def X_is_running():
    """
    Detect if X11 is available. Source:
    https://stackoverflow.com/a/1027942/1391441
    """
    from subprocess import Popen, PIPE
    p = Popen(["xset", "-q"], stdout=PIPE, stderr=PIPE)
    p.communicate()
    return p.returncode == 0


def check_all(mypath, file_end):
    """
    Checks that the necessary files are in place, the parameters stored
    in the input file are valid, and that the ranges given for the cluster
    parameters are consistent with the isochrones available before moving
    on with the code.
    """
    print('Checking input parameters...\n')

    # Check that all the essential packages are installed.
    inst_packgs_lst = pack.check()

    # Check .first_run file.
    first_run_flag = first_run.main(mypath)

    # Import here after the needed packages were checked to be present.
    from .check import params_match
    from .check import read_met_files

    # Check if input cluster files exist.
    cl_files = clusters.check(mypath, file_end)

    # Read parameters from 'params_input.dat' file. Return a dictionary
    # containing all the parameter values.
    pd = params_file.check(mypath, file_end, inst_packgs_lst)

    # Check if a new version is available.
    update.check(**pd)

    # Check running mode. If mode is 'semi', check that all cluster
    # in '/input' folder are listed.
    params_mode.check(mypath, cl_files, **pd)

    # Check that the data column indexes/names were properly given, and that
    # the magnitude and color names were properly defined.
    # If they are, store also the name of the proper isochrones folders.
    pd = params_data.check(mypath, pd)

    # Check output parameters.
    params_out.check(**pd)

    # Check structural parameters.
    params_struct.check(**pd)

    # DEPRECATED 31/10/18
    # # Check that R and rpy2 are installed, if necessary.
    # pd = params_pval.check(pd)

    # Check decontamination algorithm parameters.
    params_decont.check(cl_files, **pd)

    # Check the best synthetic cluster match parameters.
    params_match.check(**pd)

    # Filters and colors names.
    fs = ', '.join(_[1] for _ in pd['filters'])
    cs = ', '.join('(' + _[1].replace(',', '-') + ')' for _ in pd['colors'])
    print("Filter: {}".format(fs))
    print("Color:  {}\n".format(cs))

    # Check and store metallicity files.
    pd = read_met_files.check_get(pd)

    # Force matplotlib to not use any Xwindows backend. This call prevents
    # the code from crashing when used in a computer cluster. See:
    # http://stackoverflow.com/a/3054314/1391441
    if not X_is_running():
        import matplotlib
        matplotlib.use('Agg')
        print("(Force matplotlib to not use any Xwindows backend)\n")

    print("Full check done.\n\nNumber of clusters to analyze: {}\n".format(
        len(cl_files)))

    # Change these values if this is the first run, for quick processing.
    if first_run_flag:
        pd['pvalue_runs'], pd['bayesda_runs'], pd['N_bootstrap'],\
            pd['N_pop'], pd['N_gen'] = 1, 2, 2, 50, 10

    return cl_files, pd
