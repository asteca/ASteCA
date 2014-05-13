"""
@author: gabriel
"""

from time import strftime
from os.path import exists
from os import mkdir


def create_out_data_file(output_dir):
    '''
    Create output data file with headers. This will overwrite any old output
    data file already in the folder.
    '''

    out_file_name = 'ocaat_output.dat'

    # Check if output folder exists, if not create it.
    if not exists(output_dir):
        mkdir(output_dir)

    # Check if file already exists. If it does append new values instead of
    # deleting/creating the file again.
    try:
        # File already exists -> don't create a new one and append new lines.
        with open(output_dir + out_file_name):
            pass
        print 'Output data file already exists.'
    # File doesn't exist -> create new one.
    except IOError:
        print 'Output data file created.'
        out_data_file = open(output_dir + out_file_name, 'w')
        now_time = strftime("%Y-%m-%d %H:%M:%S")
        out_data_file.write("#\n\
# [%s]\n\
#\n\
# NAME: Cluster's name.\n\
# c_x[px]: Cluster's x center coordinate in pixels.\n\
# c_y[px]: Cluster's y center coordinate in pixels.\n\
# r_cl[px]: Cluster's radius in pixels.\n\
# r_c[px]: Core radius (3-P King profile) in pixels.\n\
# r_t[px]: Tidal radius (3-P King profile) in pixels.\n\
#\n\
# cont_ind: 'Contamination index' is a  measure of the contamination of field\n\
#           stars in the cluster region. The closer to 1, the more\n\
#           contaminated the cluster region is. E.g.: a value of 0.5 means\n\
#           one should expect to find the same number of field stars and\n\
#           cluster members inside the cluster region. A value of 1 means\n\
#           that *all* of the stars in the cluster region are expected to\n\
#           be field stars.\n\
# memb: Approximate number of cluster's members assuming a uniform\n\
#       background.\n\
# memb_k: Approximate number of cluster's members obtained integrating the\n\
#         fitted 3-P King profile (if it converged).\n\
# prob_cl: Statistical comparision of cluster vs field KDEs. It is obtained\n\
#          as 1 minus the overlap area between the KDEs. If the KDEs are\n\
#          very similar this value will be low indicating the overdensity is\n\
#          probably not a true cluster.\n\
# CCC: Concordance correlation coefficient (inverted). Measures the agreement\n\
#      between the quantiles and the identity line and gives an idea of how \n\
#      similar the shapes of the KDEs are. A value close to 1 means bad\n\
#      agreement, ie: the shapes of the KDEs are very different. Low values\n\
#      of prob_cl and CCC imply that the overdensity has little chance of\n\
#      being a true cluster.\n\
# mag_int: Integrated color magnitude for all stars inside the cluster\n\
#          radius, except those that were rejected due to large errors.\n\
# met: Metallicity value (z) obtained via synthetic cluster fitting.\n\
# e_m: Metallicity error.\n\
# age: log(age) for age in Gyr, idem metallicity.\n\
# e_a: log(age) error.\n\
# E(B-V): extinction, idem metallicity.\n\
# e_E: Extinction error.\n\
# dist: Distance modulus, idem metallicity.\n\
# e_d: Distance error.\n\
#\n\
# M1 (flag_center_manual): Indicates that the center was set manually.\n\
# M2 (flag_radius_manual): Indicates that the radius was set manually.\n\
# M3 (rjct_errors_fit): Indicates that all stars with errors < e_max were\n\
#    accepted, meaning that the automatic error rejecting fit was poor.\n\
#\n\
# f1 (flag_center_med): Either median cluster's central coordinates (obtained\n\
#    using all bin widths) is more than 10%% away from the values obtained \n\
#    with the min bin width.\n\
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
# f7 (flag_num_memb_low): The number of approximate cluster members is < 10.\n\
#\n\
# FC (flags count): Sum of all the flags values. The bigger this value the\n\
#    more likely it is that there's a problem with the frame, ie: no cluster,\n\
#    more than one cluster present in the frame, variable or too crowded\n\
#    field, etc.\n\
#\n\
#NAME            c_x[px] c_y[px]   e_c[px] r_cl[px]  e_r[px] r_c[px] e_rc[px] \
r_t[px] e_rt[px] cont_ind memb memb_k prob_cl   CCC mag_int     met     e_m   \
age   e_a  E(B-V)   e_E   dist   e_d  M1 M2 M3  f1 f2 f3 f4 f5 f6 f7  \
FC\n" % now_time)
        out_data_file.close()

    return out_file_name