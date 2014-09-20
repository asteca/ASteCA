# -*- coding: utf-8 -*-
"""
Created on Fri Apr 04 2014

@author: gabriel
"""

from os.path import join
from os import listdir, walk


def in_clusters(mypath):
    '''
    Store the paths and names of all the input clusters stored in the
    input_dir folder.
    '''

    # Path where cluster data files are stored.
    input_dir = join(mypath, 'input/')

    # Store subdir names and file names in cl_files.
    cl_files = []
    for root, dirs, files in walk(input_dir):

        # Don't read this sub-folder so it can be used as a container.
        if 'dont_read' in dirs:
            dirs.remove('dont_read')

        # For files not in sub-dirs.
        if files and root == input_dir:
            for name in files:
                # Don't attempt to read membership data files.
                if not name.endswith('_memb.dat'):
                    cl_files.append(['', name])

        # For files in sub-dirs.
        if dirs:
            for subdir in dirs:
                for name in listdir(join(input_dir, subdir)):
                    # Don't attempt to read membership data files.
                    if not name.endswith('_memb.dat'):
                        cl_files.append([subdir, name])

    return cl_files