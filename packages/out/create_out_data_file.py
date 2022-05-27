
import numpy as np
from time import strftime, sleep
from os.path import isfile
from .._version import __version__


txt = """#
# [ASteCA {}]
#
# Created:  [{}]
# Modified: [{}]
#
# NAME     : Cluster's name
# c_x      : x center coordinate
# c_y      : y center coordinate
# r_cl     : radius
# r_c      : Core radius (King profile)
# r_t      : Tidal radius (King profile)
# CI       : Contamination index
# n_memb_k : Number of members obtained integrating the fitted King profile
# n_memb   : Number of members assuming a uniform background
# a_f      : Fraction of cluster's area that is present in frame
#
# Parameter values are in the sense: mean, median, mode, 16th percentile,
# 84th percentile, and STDDEV.
#
# M  : total initial mass
# z  : metallicity
# a  : log(age)
# bf : binary fraction
# B  : beta
# Av : Visual absorption
# DR : differential reddening
# Rv : ratio of total to selective absorption
# dm : distance modulus
#
"""
txt += "NAME,c_x,c_y,r_cl,r_c,rc16,rc84,r_t,rt16,rt84,CI,n_memb_k,n_memb,a_f,"
for par in ('M', 'z', 'a', 'bf', 'B', 'Av', 'DR', 'Rv', 'dm'):
    txt += "{}_mean,{}_median,{}_mode,{}_16th,{}_84th,{}_std,".format(
        *[par] * 6)
# Remove last comma, add new line
txt = txt[:-1] + "\n"


def main(npd):
    """
    Create output data file with headers. This will not overwrite any old
    output data file already in the folder.
    """
    out_file_name = npd['out_file_name']
    # Current time and date.
    now_time = strftime("%Y-%m-%d %H:%M:%S")

    # Check if file already exists. If it does append new values instead of
    # deleting/creating the file again.
    if isfile(out_file_name):

        # File already exists -> don't create a new one and replace old lines.
        while True:
            with open(out_file_name, 'r') as f:
                # Read file into data var.
                data = f.readlines()

            try:
                # Modify these two lines
                data[1] = '# [ASteCA {}]\n'.format(__version__)
                data[4] = '# Modified: [{}]\n'.format(now_time)
                break
            except IndexError:
                # Wait a random number of seconds (max 10) before reading the
                # file again. This is here to avoid several parallel runs from
                # finishing all at once and overlapping each other.
                sleep(np.random.uniform(10))

        # Write everything back.
        with open(out_file_name, 'w') as f:
            f.writelines(data)

    # File doesn't exist -> create new one.
    else:
        with open(out_file_name, 'w') as out_data_file:
            out_data_file.write(txt.format(__version__, now_time, now_time))
