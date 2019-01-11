
from time import strftime
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
        with open(out_file_name, 'r') as f:
            # Read file into data var.
            data = f.readlines()

        # Modify these two lines
        data[1] = '# [ASteCA {}]\n'.format(__version__)
        data[4] = '# Modified: [{}]\n'.format(now_time)

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
# e_rcl: Cluster's radius error.\n\
# r_c: Core radius (3-P King profile).\n\
# e_rc: Core radius error.\n\
# r_t: Tidal radius (3-P King profile).\n\
# e_rt: Tidal radius error.\n\
# kcp: King's profile concentration parameter.\n\
#\n\
# CI: 'Contamination index' is a  measure of the contamination of field\n\
#      stars in the cluster region. The closer to 1, the more contaminated \n\
#      the cluster region is.\n\
# n_memb_k: Approximate number of cluster's members obtained integrating the\n\
#         fitted 3-P King profile (if it converged).\n\
# n_memb: Approximate number of cluster's members assuming a uniform\n\
#       background.\n\
# n_memb_da: Approximate number of cluster's members obtained via the DA\n\
#          algorithm.\n\
# memb_par: Members parameter comparing the approximate structural number of\n\
#           members ('n_memb') with the approximate photometric number of\n\
#           members ('n_memb_da').\n\
# a_f: Fraction of cluster's area that is present in frame.\n\
#\n\
# Parameters values are in the sense: mean, MAP/ML, median.\n\
# Parameters uncertainties are: 16th, 84th percentiles.\n\
# z: Metallicity value.\n\
# a: log(age).\n\
# e: extinction E(B-V).\n\
# d: Distance modulus.\n\
# M: Total initial mass.\n\
# b: Binary fraction.\n\
#\n\
# M1 Indicates that the center was set manually.\n\
# M2 Indicates that the radius was set manually.\n\
#\n\
# f1 The standard deviation for either center coordinate is larger than 10%\n\
#    of that coordinate's range.\n\
# f2 The background value is smaller than a third of the maximum radial\n\
#    density value.\n\
# f3 Not enough points found stabilized around the background value -->\n\
#    clust_rad was set to the middle value in the density profile.\n\
# f4 The delta range around the background used to attain the stable\n\
#    condition to determine the radius is greater than 10%%. This indicates\n\
#    a possible variable background.\n\
# f5 The process to fit a 3-P King profile to the density points did not\n\
#    converge or did so to a tidal radius beyond the ranges of the frame.\n\
# f6 The number of approximate structural cluster members ('n_memb') is <10.\n\
# f7 The number of approximate structural and photometric cluster members\n\
#    differ greatly --> abs(n_memb_par) > 0.33.\n\
#\n\
# FC (flags count): Sum of all the flags values. The bigger this value the\n\
#    more likely it is that there's a problem with the frame, ie: no\n\
#    cluster, more than one cluster present in the frame, variable or too\n\
#    crowded field, etc.\n\
#\n\
#NAME                 c_x      c_y     r_cl    e_rcl      r_c     e_rc      \
r_t     e_rt      kcp      CI   n_memb_k     n_memb  n_memb_da  memb_par     \
a_f     \
z_mean      z_MAP   z_median       16th       84th        std     \
a_mean      a_MAP   a_median       16th       84th        std     \
e_mean      e_MAP   e_median       16th       84th        std     \
d_mean      d_MAP   d_median       16th       84th        std     \
M_mean      M_MAP   M_median       16th       84th        std     \
b_mean      b_MAP   b_median       16th       84th        std      \
M1 M2  f1 f2 f3 f4 f5 f6 f7  FC\n".format(__version__, now_time, now_time))
            print('Output data file created.')
