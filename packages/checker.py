
from check import pack
from check import update
from check import clusters
from check import params_file
from check import params_input_struct
from check import params_input_pval
from check import params_input_decont


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

    # Check if input cluster files exist.
    cl_files = clusters.check(mypath, file_end)

    # Read parameters from 'params_input.dat' file.
    pd = params_file.check(mypath, file_end)

    # Check that R and rpy2 are installed, if necessary.
    R_in_place = params_input_pval.check(inst_packgs_lst, **pd)

    # Check if a new version is available.
    if pd['up_flag']:
        update.check()

    # Check that structural parameters are properly given.
    params_input_struct.check(mypath, cl_files, **pd)

    # Define dictionary of accepted binning methods.
    bin_methods_dict = {'blocks', 'knuth', 'scott', 'freedman', 'sturges',
                        'sqrt', 'bb'}

    # Check decontamination algorithm parameters.
    params_input_decont.check(cl_files, bin_methods_dict, **pd)

    # Check the best synthetic cluster match parameters.
    # Import here after the needed packages were checked to be present, since
    # this imports numpy.
    from check import params_input_match
    params_input_match.check(bin_methods_dict, **pd)

    # Check and store metallicity files.
    from check import read_met_files
    ip_list = read_met_files.check_get(pd)

    print("Full check done. Clusters to process: {}\n".format(
        len(cl_files)))

    # Force matplotlib to not use any Xwindows backend. This call prevents
    # the code from crashing when used in a computer cluster. See:
    # http://stackoverflow.com/a/3054314/1391441
    if pd['flag_back_force']:
        import matplotlib
        matplotlib.use('Agg')
        print("(Force matplotlib to not use any Xwindows backend)\n")

    return cl_files, pd, ip_list, R_in_place
