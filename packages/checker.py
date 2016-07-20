
from check import pack
from check import update
from check import clusters
from check import params_file
from check import params_input_struct
from check import params_input_pval
from check import params_input_decont
from check import params_input_match
from check import read_met_files


def check_all(mypath, file_end, cl_files):
    """
    Checks that the necessary files are in place, the parameters stored
    in the input file are valid and that the ranges given for the cluster
    parameters are consistent with the isochrones available before moving
    on with the code.
    """
    print('Checking input parameters...\n')

    # Check that all the essential packages are installed.
    inst_packgs_lst = pack.check()

    # Check if input cluster files exist.
    clusters.check(cl_files, file_end)

    # Check if 'params_input.dat' file exists. Initialize global variables.
    flag_updt, flag_back_force = params_file.check(mypath, file_end)

    # Check if a new version is available.
    if flag_updt:
        update.check()

    # Check that structural parameters are properly given.
    params_input_struct.check(mypath, cl_files)

    # Check that R and rpy2 are installed, if necessary.
    R_in_place = params_input_pval.check(inst_packgs_lst)

    # Define dictionary of accepted binning methods.
    bin_methods_dict = {'blocks', 'knuth', 'scott', 'freedman', 'sturges',
                        'sqrt', 'bb'}

    # Check decontamination algorithm parameters.
    params_input_decont.check(cl_files, bin_methods_dict)

    # Check the best synthetic cluster match parameters.
    params_input_match.check(bin_methods_dict)

    # Check and store metallicity files.
    ip_list = read_met_files.check_get()

    print 'Full check done.\n'

    # Force matplotlib to not use any Xwindows backend. This call prevents
    # the code from crashing when used in a computer cluster. See:
    # http://stackoverflow.com/a/3054314/1391441
    if flag_back_force:
        import matplotlib
        matplotlib.use('Agg')
        print("(Force matplotlib to not use any Xwindows backend)\n")

    return ip_list, R_in_place
