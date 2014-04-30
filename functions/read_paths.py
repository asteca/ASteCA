# -*- coding: utf-8 -*-
"""
Created on Fri Apr 04 2014

@author: gabriel
"""

from os.path import join
from os import listdir, walk


def read_paths(input_dir):
    '''
    Store the paths and names of all the input clusters stored in the
    input_dir folder.
    '''

    # Store subdir names [0] and file names [1] inside each subdir.
    dir_files = [[], []]
    for root, dirs, files in walk(input_dir):
        if dirs:
            for subdir in dirs:
                for name in listdir(join(input_dir, subdir)):
                    # Check to see if it's a valid data file.
                    if name.endswith(('.DAT', '.MAG', '.OUT', '.TEX', '.TXT',
                    '.dat', '.mag', '.out', '.tex', '.txt')) and \
                    not name.endswith('_memb.dat'):
                        dir_files[0].append(subdir)
                        dir_files[1].append(name)

    return dir_files