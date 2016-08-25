
from check import pack
from check import update
from check import clusters
from check import params_file
from check import photom_names
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

    # Read parameters from 'params_input.dat' file. Return a dictionary
    # containing all the parameter values.
    pd = params_file.check(mypath, file_end)

    # Check that the magnitude and color names were properly given.
    # If they are, store also the name of the proper isochrones folders.
    pd = photom_names.check(mypath, pd)

    # Check that R and rpy2 are installed, if necessary.
    pd = params_input_pval.check(inst_packgs_lst, pd)

    # Check if a new version is available.
    update.check(**pd)

    # Check that structural parameters are properly given.
    params_input_struct.check(mypath, cl_files, **pd)

    # Check decontamination algorithm parameters.
    params_input_decont.check(cl_files, **pd)

    # Print info about tracks.
    # Map isochrones set selection to proper name.
    iso_select = pd['ps_params'][2]
    iso_print = pd['tracks_dict'].get(iso_select)
    # Extract photometric system used,m from the isochrone's folder name.
    syst = pd['ps_params'][0].split('_', 1)[1]
    print("Process {} theoretical isochrones".format(iso_print))
    print("in the '{}' photometric system.\n".format(syst))

    # Check the best synthetic cluster match parameters.
    # Import here after the needed packages were checked to be present, since
    # this imports numpy.
    from check import params_input_match
    params_input_match.check(**pd)

    # Check and store metallicity files.
    from check import read_met_files
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

    return cl_files, pd
