
from os.path import join, exists
from os import makedirs, extsep


def main(cl_file, done_dir, mode, **kwargs):
    '''
    Generate names and paths to be used by several functions.
    '''
    # Hardcoded in/out folder names.
    out_fold = 'output'

    # Split cluster's file into sections separated by dots.
    cl_split = cl_file[-1].split(extsep)
    # Join all the sections except the last one (the extension) and store the
    # cluster's clean name.
    clust_name = '.'.join(cl_split[:-1])

    # Generate hardcoded file names and paths.
    data_file = join(*cl_file)
    # Path to membership probabilities file if it exists.
    memb_file = join(cl_file[0], cl_file[1], cl_file[2],
                     clust_name + '_memb.dat')
    # Root output dir.
    output_dir = join(cl_file[0], out_fold)
    # Output subdir and 'done' dir.
    if cl_file[1][-2:].isdigit():
        # If the 'input' folder has numbers at the end.
        output_subdir = join(output_dir, cl_file[1], cl_file[2])
        dst_dir = join(done_dir, cl_file[1], cl_file[2])
    else:
        output_subdir = join(output_dir, cl_file[2])
        dst_dir = join(done_dir, cl_file[2])

    # Generate output dir/subdir if it doesn't exist.
    if not exists(output_subdir):
        makedirs(output_subdir)

    memb_file_out = join(output_subdir, clust_name + '_memb.dat')
    synth_file_out = join(output_subdir, clust_name + '_synth.dat')
    write_name = join(cl_file[2], clust_name)
    out_file_name = join(output_dir, 'asteca_output.dat')

    print("Analyzing cluster {} ({} mode).".format(clust_name, mode))

    npd = {
        'clust_name': clust_name, 'data_file': data_file,
        'memb_file': memb_file, 'output_dir': output_dir,
        'out_file_name': out_file_name, 'output_subdir': output_subdir,
        'dst_dir': dst_dir, 'memb_file_out': memb_file_out,
        'synth_file_out': synth_file_out, 'write_name': write_name}
    return npd
