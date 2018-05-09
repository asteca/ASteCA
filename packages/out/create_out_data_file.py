
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
        data[3] = '# Modified: [{}]\n'.format(now_time)

        # Write everything back.
        with open(out_file_name, 'w') as f:
            f.writelines(data)

    # File doesn't exist -> create new one.
    else:
        with open(out_file_name, 'w') as out_data_file:
            out_data_file.write("#\n\
# [ASteCA {}]\n\
#\n\
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
# prob_cl: Statistical comparison of cluster vs field KDEs. It is obtained\n\
#          as 1 minus the overlap area between the KDEs. If the KDEs are\n\
#          very similar this value will be low indicating the overdensity is\n\
#          probably not a true cluster.\n\
# int_col: Integrated color magnitude for all stars inside the cluster\n\
#          radius, except those that were rejected due to large errors.\n\
# met: Metallicity value (z).\n\
# e_m: Metallicity error.\n\
# age: log(age).\n\
# e_a: log(age) error.\n\
# E(B-V): extinction.\n\
# e_E: Extinction error.\n\
# dist: Distance modulus.\n\
# e_d: Distance error.\n\
# M_i: Total initial mass.\n\
# e_M: Mass error.\n\
# bin_fr: Binary fraction.\n\
# e_bf: Binary fraction error.\n\
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
#NAME                 c_x      c_y     r_cl    e_rcl      \
r_c     e_rc      r_t     e_rt      kcp      CI   n_memb_k     n_memb  \
n_memb_da  memb_par     a_f  prob_cl  int_col      met      e_m      age      \
e_a   E(B-V)      e_E     dist      e_d      M_i      e_M   bin_fr     \
e_bf      \
M1 M2  f1 f2 f3 f4 f5 f6 f7 f8  FC\n".format(__version__, now_time))
            print 'Output data file created.'
