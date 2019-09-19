
import platform
from .check import pack
from .check import first_run
from .check import update
from .check import clusters
from .check import params_file
from .check import params_data
from .check import params_out
from .check import params_struct
from .check import params_decont


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

    # Check output parameters.
    params_out.check(**pd)

    # Check structural parameters.
    params_struct.check(**pd)

    # Check decontamination algorithm parameters.
    params_decont.check(cl_files, **pd)

    # Check the best synthetic cluster match parameters.
    pd = params_match.check(mypath, pd)

    # Filters and colors names.
    fs = ', '.join(_[1] for _ in pd['filters'])
    cs = ', '.join('(' + _[1].replace(',', '-') + ')' for _ in pd['colors'])
    print("Filter: {}".format(fs))
    print("Color:  {}\n".format(cs))

    # Check and store metallicity files.
    pd = read_met_files.check_get(pd)

    # Force matplotlib to not use Xwindows backend. This call prevents
    # the code from crashing when used in a computer cluster. See:
    # http://stackoverflow.com/a/3054314/1391441
    if not X_is_running():
        import matplotlib
        matplotlib.use('Agg')
        print("(Force matplotlib to not use Xwindows backend)\n")

    print("Full check done.\n\nNumber of clusters to analyze: {}\n".format(
        len(cl_files)))

    return cl_files, pd


def X_is_running():
    """
    Detect if X11 is available. Source:
    https://stackoverflow.com/a/1027942/1391441
    """
    if platform.system() == 'Linux':
        from subprocess import Popen, PIPE
        p = Popen(["xset", "-q"], stdout=PIPE, stderr=PIPE)
        p.communicate()
        return p.returncode == 0
    else:
        # If this is not a Linux system, assume that it is either Mac OS or
        # Windows, and thus assume that a windows system is present.
        return True
