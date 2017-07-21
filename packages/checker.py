
from check import pack
from check import update
from check import clusters
from check import params_file
from check import semi_list
from check import photom_names
from check import params_struct
from check import params_errors
from check import params_pval
from check import params_decont


def check_all(mypath, first_run_flag, file_end):
    """
    Checks that the necessary files are in place, the parameters stored
    in the input file are valid, and that the ranges given for the cluster
    parameters are consistent with the isochrones available before moving
    on with the code.
    """
    print('Checking input parameters...\n')

    # Check that all the essential packages are installed.
    inst_packgs_lst = pack.check()
    # Import here after the needed packages were checked to be present.
    from check import params_match
    from check import read_met_files

    # Check if input cluster files exist.
    cl_files = clusters.check(mypath, file_end)

    # Read parameters from 'params_input.dat' file. Return a dictionary
    # containing all the parameter values.
    pd = params_file.check(mypath, file_end)

    # If mode is 'semi', check that all cluster in '/input' folder are listed.
    semi_list.check(cl_files, **pd)

    # Check if a new version is available.
    update.check(**pd)

    # Check that the magnitude and color names were properly given.
    # If they are, store also the name of the proper isochrones folders.
    pd = photom_names.check(mypath, pd)

    # Check that mode, coordinates, figure, and structural parameters
    # are properly given.
    params_struct.check(mypath, **pd)

    # Check that the errors module parameters are correct.
    params_errors.check(**pd)

    # Check that R and rpy2 are installed, if necessary.
    pd = params_pval.check(inst_packgs_lst, pd)

    # Check decontamination algorithm parameters.
    params_decont.check(cl_files, **pd)

    # Check the best synthetic cluster match parameters.
    params_match.check(**pd)

    # Check and store metallicity files.
    pd = read_met_files.check_get(pd)

    # Force matplotlib to not use any Xwindows backend. This call prevents
    # the code from crashing when used in a computer cluster. See:
    # http://stackoverflow.com/a/3054314/1391441
    if pd['flag_back_force']:
        import matplotlib
        matplotlib.use('Agg')
        print("(Force matplotlib to not use any Xwindows backend)\n")

    print("Full check done.\n\nNumber of clusters to analyze: {}\n".format(
        len(cl_files)))

    # Change these values if this is the first run, for quick processing.
    if first_run_flag:
        pd['pvalue_mode'], pd['pvalue_runs'], pd['bayesda_mode'],\
            pd['bayesda_runs'], pd['N_bootstrap'], pd['N_pop'], pd['N_gen'] =\
            'manual', 2, 'manual', 2, 2, 50, 10

    return cl_files, pd
