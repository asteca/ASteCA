# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 14:03:44 2014

@author: gabriel
"""

from os.path import join, exists
from os import makedirs, extsep


def names_paths(cl_file):
    '''
    Generate names and paths to be used by several functions.
    '''

    # Hardcoded in/out folder names.
    out_fold = 'output'

    # Store cluster's name without extension.
    clust_name = cl_file[-1].split(extsep)[0]

    # Generate hardcoded file names and paths.
    data_file = join(*cl_file)
    # Path to membership probabilities file if it exists.
    memb_file = join(cl_file[0], cl_file[1], cl_file[2],
        clust_name + '_memb.dat')
    # Root output dir.
    output_dir = join(cl_file[0], out_fold)
    # Output subdir.
    output_subdir = join(output_dir, cl_file[2])

    # Generate output dir/subdir if it doesn't exist.
    if not exists(output_subdir):
        makedirs(output_subdir)

    memb_file_out = join(output_subdir, clust_name + '_memb.dat')
    write_name = join(cl_file[2], clust_name)

    return clust_name, data_file, memb_file, output_dir, output_subdir,\
    memb_file_out, write_name