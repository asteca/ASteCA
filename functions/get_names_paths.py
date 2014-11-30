
from os.path import join, exists
from os import makedirs, extsep


def names_paths(mypath, cl_file):
    '''
    Generate names and paths to be used by several functions.
    '''

    # Hardcoded in/out folder names.
    in_fold = 'input'
    out_fold = 'output'

    # Store cluster's name without extension.
    clust_name = cl_file[1].split(extsep)[0]

    # Generate hardcoded file names and paths.
    data_file = join(mypath, in_fold, *cl_file)
    memb_file = join(mypath, in_fold, cl_file[0], clust_name + '_memb.dat')
    output_dir = join(mypath, out_fold)
    output_subdir = join(output_dir, cl_file[0])

    # Generate output dir/subdir if it doesn't exist.
    if not exists(output_subdir):
        makedirs(output_subdir)

    memb_file_out = join(output_subdir, clust_name + '_memb.dat')
    write_name = join(cl_file[0], clust_name)

    return clust_name, data_file, memb_file, output_dir, output_subdir,\
    memb_file_out, write_name