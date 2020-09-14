
from os.path import join, exists
from os import makedirs, extsep


def main(cl_file, **kwargs):
    '''
    Generate names and paths to be used by several functions.
    '''
    # Hard-coded output folder name.
    out_fold = 'output'

    # Extract cluster's name from file.
    clust_name = get_clust_name(cl_file)

    # Generate hardcoded file names and paths.
    data_file = join(*cl_file)
    # Path to membership probabilities file if it exists.
    memb_file = memb_file_name(cl_file)
    # Root output dir.
    output_dir = join(cl_file[0], out_fold)
    # Output subdir and 'done' dir.
    if cl_file[1][-2:].isdigit():
        # If the 'input' folder has numbers at the end.
        output_subdir = join(output_dir, cl_file[1], cl_file[2], clust_name)
    else:
        output_subdir = join(output_dir, cl_file[2], clust_name)

    # Generate output dir/subdir if it doesn't exist.
    if not exists(output_subdir):
        makedirs(output_subdir)

    memb_file_out = join(output_subdir, clust_name + '_memb.dat')
    mcmc_file_out = join(output_subdir, clust_name + '_mcmc.pickle')
    synth_file_out = join(output_subdir, clust_name + '_synth.dat')
    mass_file_out = join(output_subdir, clust_name + '_mass.dat')
    write_name = join(cl_file[2], clust_name)
    out_file_name = join(output_dir, 'asteca_output.dat')
    params_out = join(output_subdir, clust_name + '_params_input.dat')

    print("Analyzing cluster {}".format(clust_name))

    npd = {
        'clust_name': clust_name, 'data_file': data_file,
        'memb_file': memb_file, 'output_dir': output_dir,
        'out_file_name': out_file_name, 'output_subdir': output_subdir,
        'memb_file_out': memb_file_out, 'synth_file_out': synth_file_out,
        'mass_file_out': mass_file_out, 'write_name': write_name,
        'mcmc_file_out': mcmc_file_out, 'params_out': params_out}
    return npd


def get_clust_name(cl_file):
    """
    Extract cluster's name from file.
    """
    # Split cluster's file into sections separated by dots.
    cl_split = cl_file[-1].split(extsep)
    # Join all the sections except the last one (the extension) and store the
    # cluster's clean name.
    clust_name = '.'.join(cl_split[:-1])

    return clust_name


def memb_file_name(cl_file):
    """
    Name of file with membership data.
    """
    # Call function again to avoid passing 'clust_name'. This allows the
    # module 'params_input_decont' to use this function.
    clust_name = get_clust_name(cl_file)
    memb_file = join(cl_file[0], cl_file[1], cl_file[2],
                     clust_name + '_memb.dat')
    return memb_file
