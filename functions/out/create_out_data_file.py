"""
@author: gabriel
"""

from time import strftime
from os.path import exists, join
from os import mkdir
from functions import __version__


def create_out_data_file(output_dir):
    '''
    Create output data file with headers. This will not overwrite any old
    output data file already in the folder.
    '''

    # Output file name.
    out_file_name = join(output_dir, 'asteca_output.dat')

    # Check if output folder exists, if not create it.
    if not exists(output_dir):
        mkdir(output_dir)

    # Check if file already exists. If it does append new values instead of
    # deleting/creating the file again.
    try:
        # File already exists -> don't create a new one and append new lines.
        with open(out_file_name):
            pass
    # File doesn't exist -> create new one.
    except IOError:
        print 'Output data file created.'
        out_data_file = open(out_file_name, 'w')
        now_time = strftime("%Y-%m-%d %H:%M:%S")
        out_data_file.write("#\n\
# ASteCA {}\n\
#\n\
# Created: [{}]\n\
#\n\
# NAME: Cluster's name.\n\
# c_x: Cluster's x center coordinate.\n\
# e_x: Cluster's x center coordinate error.\n\
# c_y: Cluster's y center coordinate.\n\
# e_y: Cluster's x center coordinate error.\n\
# r_cl: Cluster's radius.\n\
# e_rcl: Cluster's radius error.\n\
# r_c: Core radius (3-P King profile).\n\
# e_rc: Core radius error.\n\
# r_t: Tidal radius (3-P King profile).\n\
# e_rt: Tidal radius error.\n\
#\n\
# CI: 'Contamination index' is a  measure of the contamination of field\n\
#      stars in the cluster region. The closer to 1, the more contaminated \n\
#      the cluster region is.\n\
# memb: Approximate number of cluster's members assuming a uniform\n\
#       background.\n\
# memb_k: Approximate number of cluster's members obtained integrating the\n\
#         fitted 3-P King profile (if it converged).\n\
# prob_cl: Statistical comparision of cluster vs field KDEs. It is obtained\n\
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
# M1 (flag_center_manual): Indicates that the center was set manually.\n\
# M2 (flag_radius_manual): Indicates that the radius was set manually.\n\
#\n\
# f1 (flag_center_med): Either median cluster's central coordinates (obtained\n\
#    using several standard deviation smoothing values) is more than 10%% \n\
#    away from the central position assigned.\n\
# f2 (flag_center_std): The standard deviation for either center coordinate\n\
#    is larger than 10%% of the coordinate's value.\n\
# f3 (flag_delta_total): The background value is smaller than a third of the\n\
#    maximum radial density value.\n\
# f4 (flag_not_stable): Not enough points found stabilized around the\n\
#    background value -> r = middle value of density profile.\n\
# f5 (flag_delta): The delta range around the background used to attain the\n\
#    stable condition to determine the radius is greater than 10%%. This\n\
#    indicates a possible variable background.\n\
# f6 (flag_king_no_conver): The process to fit a 3-P King profile to the\n\
#    density points did not converge or did so to a tidal radius beyond the\n\
#    ranges of the frame.\n\
# f7 (err_all_fallback): no error rejection was possible.\n\
# f8 (err_max_fallback): the function had to fall back to the 'e_max'-based\n\
#    rejection method since the selected one failed.\n\
# f9 (flag_num_memb_low): The number of approximate cluster members is < 10.\n\
#\n\
# FC (flags count): Sum of all the flags values. The bigger this value the\n\
#    more likely it is that there's a problem with the frame, ie: no cluster,\n\
#    more than one cluster present in the frame, variable or too crowded\n\
#    field, etc.\n\
#\n\
#NAME                 c_x      e_x      c_y      e_y     r_cl    e_rcl      \
r_c     e_rc      r_t     e_rt       CI     memb  memb_k   prob_cl  \
int_col      met      e_m      age      e_a   E(B-V)      e_E     dist      \
e_d      M_i      e_M   bin_fr     e_bf      M1 M2  f1 f2 f3 f4 f5 f6 f7 f8 \
f9  FC\n".format(__version__, now_time))
        out_data_file.close()

    return out_file_name