"""
@author: gabriel
"""

from time import strftime

def create_out_data_file(output_dir):
    '''
    Create output data file with headers. This will overwrite any old output
    data file already in the folder.
    '''
    
    # Check if file already exists. If it does append new values instead of
    # deleting/creating the file again.
    try:
        # File already exists -> don't create a new one and append new lines.
        with open(output_dir+'data_output'): pass
        print('\nOutput data file already exists.')
    # File doesn't exist -> create new one.
    except IOError:
        out_data_file = open(output_dir+'data_output','w')
        now_time = strftime("%Y-%m-%d %H:%M:%S")
        out_data_file.write("#\n\
# [%s]\n\
#\n\
# NAME: Cluster's name.\n\
# c_x[px]: Cluster's x center coordinate in pixels.\n\
# c_y[px]: Cluster's y center coordinate in pixels.\n\
# R_cl[px]: Cluster's radius in pixels.\n\
# R_c[px]: Core radius (3-P King profile) in pixels.\n\
# R_t[px]: Tidal radius (3-P King profile) in pixels.\n\
#\n\
# cont_ind: 'Contamination index' is a  measure of the contamination of field\n\
#           stars in the cluster region. The closer to 1, the more contaminated\n\
#           the cluster region is. E.g.: a value of 0.5 means I should expect\n\
#           to find the same number of field stars and cluster members inside the\n\
#           cluster region.\n\
# memb: Approximate number of cluster's members assuming a uniform background.\n\
# memb_k: Approximate number of cluster's members obtained integrating the fitted\n\
#         3-P King profile (if it converged).\n\
# p_value: Statistical comparision of cluster vs field KDEs. A value over 0.05\n\
#          (5%%) means there is little chance this is an actual cluster (ie: the\n\
#          null hypothesis that both CMD realizations come from the same distribution\n\
#          can't be rejected).\n\
# mag_int: Integrated magnitude value for all stars inside the cluster radius, except\n\
#          those that were rejected due to large errors.\n\
#\n\
# M1 (flag_center_manual): Indicates that the center was set manually.\n\
# M2 (flag_radius_manual): Indicates that the radius was set manually.\n\
# M3 (flag_errors_manual): Indicates that all stars with errors < 0.3 were\n\
#    accepted, meaning that the automatic error rejecting fit was poor.\n\
#\n\
# f1 (flag_center): Either the x or y coordinate assigned as center deviates\n\
#    more than 50px from the mean value obtained using all bin widths.\n\
# f2 (flag_std_dev): Standard deviation of either x or y of the mean center\n\
#    coordinates is > 50px.\n\
# f3 (flag_bin_count): When calculating the background value, the inner limit\n\
#    pushed all the bins outside the frame. This could be indicative of a small\n\
#    field given the size of the cluster.\n\
# f4 (flag_delta_total): The background value is smaller than a third of the\n\
#    maximum radial density value.\n\
# f5 (flag_not_stable): Not enough points found stabilized around the background\n\
#    value -> r = 500px assigned.\n\
# f6 (flag_rad_500): Radius obtained was bigger than 500px -> trimmed to 500px.\n\
# f7 (flag_delta): The delta range around the background used to attain the\n\
#    stable condition to determine the radius is greater than 10%%. This\n\
#    indicates a possible variable background.\n\
# f8 (flag_delta_points): The number of points that fall outside the delta\n\
#    range is higher than the number of points to the left of the radius. This\n\
#    indicates a possible variable background.\n\
# f9 (flag_king_no_conver): The process to fit a 3-P King profile to the density\n\
#    points did not converge.\n\
# f10 (flag_num_memb_low): The number of approximate cluster members is < 10.\n\
# f11 (flag_no_memb): The number of approximate cluster members found was < 0 so\n\
#     a value of 1/10th the total number of stars with r<R_cl or 9, whichever\n\
#     number was larger, was set.\n\
#\n\
# FC (flags count): Sum of all the flags values. The bigger this value the more\n\
#    likely it is that there's a problem with the frame, ie: no cluster, more \n\
#    than one cluster present in the frame, variable or too crowded field, etc.\n\
#\n\
#NAME            c_x[px] c_y[px] R_cl[px] R_c[px] R_t[px] cont_ind memb memb_k \
p_value mag_int  M1 M2 M3  f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11  FC\n" % now_time)
        out_data_file.close()
