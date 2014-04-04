#
# [2014-04-03 16:01:13]
#
# NAME: Cluster's name.
# c_x[px]: Cluster's x center coordinate in pixels.
# c_y[px]: Cluster's y center coordinate in pixels.
# r_cl[px]: Cluster's radius in pixels.
# r_c[px]: Core radius (3-P King profile) in pixels.
# r_t[px]: Tidal radius (3-P King profile) in pixels.
#
# cont_ind: 'Contamination index' is a  measure of the contamination of field
#           stars in the cluster region. The closer to 1, the more
#           contaminated the cluster region is. E.g.: a value of 0.5 means
#           one should expect to find the same number of field stars and
#           cluster members inside the cluster region. A value of 1 means
#           that *all* of the stars in the cluster region are expected to
#           be field stars.
# memb: Approximate number of cluster's members assuming a uniform
#       background.
# memb_k: Approximate number of cluster's members obtained integrating the
#         fitted 3-P King profile (if it converged).
# prob_cl: Statistical comparision of cluster vs field KDEs. It is obtained
#          as 1 minus the overlap area between the KDEs. If the KDEs are
#          very similar this value will be low indicating the overdensity is
#          probably not a true cluster.
# CCC: Concordance correlation coefficient (inverted). Measures the agreement
#      between the quantiles and the identity line and gives an idea of how 
#      similar the shapes of the KDEs are. A value close to 1 means bad
#      agreement, ie: the shapes of the KDEs are very different. Low values
#      of prob_cl and CCC imply that the overdensity has little chance of
#      being a true cluster.
# mag_int: Integrated magnitude value for all stars inside the cluster
#          radius, except those that were rejected due to large errors.
# met: Metallicity value (z) obtained via synthetic cluster fitting.
# e_m: Metallicity error.
# age: log(age) for age in Gyr, idem metallicity.
# e_a: log(age) error.
# E(B-V): extinction, idem metallicity.
# e_E: Extinction error.
# dist: Distance modulus, idem metallicity.
# e_d: Distance error.
#
# M1 (flag_center_manual): Indicates that the center was set manually.
# M2 (flag_radius_manual): Indicates that the radius was set manually.
# M3 (rjct_errors_fit): Indicates that all stars with errors < e_max were
#    accepted, meaning that the automatic error rejecting fit was poor.
#
# f1 (flag_center_med): Either median cluster's central coordinates (obtained
#    using all bin widths) is more than 10% away from the values obtained 
#    with the min bin width.
# f2 (flag_center_std): The standard deviation for either center coordinate
#    is larger than 10% of the coordinate's value.
# f3 (flag_delta_total): The background value is smaller than a third of the
#    maximum radial density value.
# f4 (flag_not_stable): Not enough points found stabilized around the
#    background value -> r = middle value of density profile.
# f5 (flag_delta): The delta range around the background used to attain the
#    stable condition to determine the radius is greater than 10%. This
#    indicates a possible variable background.
# f6 (flag_king_no_conver): The process to fit a 3-P King profile to the
#    density points did not converge or did so to a tidal radius beyond the
#    ranges of the frame.
# f7 (flag_num_memb_low): The number of approximate cluster members is < 10.
#
# FC (flags count): Sum of all the flags values. The bigger this value the
#    more likely it is that there's a problem with the frame, ie: no cluster,
#    more than one cluster present in the frame, variable or too crowded
#    field, etc.
#
#NAME            c_x[px] c_y[px] r_cl[px] r_c[px] e_rc[px] r_t[px] e_rt[px] cont_ind memb memb_k prob_cl   CCC mag_int     met     e_m   age   e_a  E(B-V)   e_E   dist   e_d  M1 M2 M3  f1 f2 f3 f4 f5 f6 f7  FC
