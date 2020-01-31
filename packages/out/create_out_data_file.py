
import numpy as np
from time import strftime, sleep
from os.path import isfile
from .._version import __version__


def main(npd):
    '''
    Create output data file with headers. This will not overwrite any old
    output data file already in the folder.
    '''
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
            out_data_file.write("#\n\
# [ASteCA {}]\n\
#\n\
# Created:  [{}]\n\
# Modified: [{}]\n\
#\n\
# NAME: Cluster's name.\n\
# c_x: Cluster's x center coordinate.\n\
# c_y: Cluster's y center coordinate.\n\
# r_cl: Cluster's radius.\n\
# r_c: Core radius (3-P King profile).\n\
# r_t: Tidal radius (3-P King profile).\n\
#\n\
# CI: 'Contamination index' is a  measure of the contamination of field\n\
#      stars in the cluster region. The closer to 1, the more contaminated \n\
#      the cluster region is.\n\
# n_memb_k: Approximate number of cluster's members obtained integrating the\n\
#           fitted King profile (if it converged).\n\
# n_memb: Approximate number of cluster's members assuming a uniform\n\
#         background.\n\
# n_memb_da: Approximate number of cluster's members obtained via the DA\n\
#            algorithm.\n\
# memb_par: Members parameter comparing the approximate structural number of\n\
#           members ('n_memb') with the approximate photometric number of\n\
#           members ('n_memb_da').\n\
# a_f: Fraction of cluster's area that is present in frame.\n\
#\n\
# Parameters values are in the sense: mean, MAP/ML, median, mode.\n\
# Parameters uncertainties are 16th, 84th percentiles.\n\
# z: Metallicity value.\n\
# a: log(age).\n\
# E: extinction E(B-V).\n\
# d: Distance modulus.\n\
# M: Total initial mass.\n\
# b: Binary fraction.\n\
# Nt: Number of samples used to estimate the parameters values.\n\
#\n\
NAME                                  c_x        c_y       \
r_cl         16         84        r_c         16         84        \
r_t         16         84      CI   n_memb_k     n_memb  n_memb_da  \
memb_par     a_f     \
z_mean      z_MAP   z_median     z_mode       16th       84th        std    R^2     \
a_mean      a_MAP   a_median     a_mode       16th       84th        std    R^2     \
E_mean      E_MAP   E_median     E_mode       16th       84th        std    R^2     \
d_mean      d_MAP   d_median     d_mode       16th       84th        std    R^2     \
M_mean      M_MAP   M_median     M_mode       16th       84th        std    R^2     \
b_mean      b_MAP   b_median     b_mode       16th       84th        std    \
R^2\n".format(__version__, now_time, now_time))
            print("Output data file created")
