# -*- coding: utf-8 -*-
"""
Created on Fri Apr 04 2014

@author: gabriel
"""

from os.path import join
from os import walk
import re


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

        for f in files:
            # Remove input_dir from sub-dir path.
            subdir = re.sub(r'^' + re.escape(input_dir), '', root)
            # Don't attempt to read membership or .mf files.
            if not f.endswith(('_memb.dat', '.md')):
                cl_files.append([subdir, f])

    # Return sorted list by cluster file name.
    cl_files.sort(key=lambda x: x[1].lower())

    return cl_files