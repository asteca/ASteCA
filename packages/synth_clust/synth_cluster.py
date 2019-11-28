
from . import move_isochrone
from . import cut_max_mag
from . import mass_distribution
from . import mass_interp
from . import binarity
from . import completeness_rm
from . import add_errors


def main(
    isochrone, model_proper, err_max, err_lst, completeness, max_mag_syn,
        st_dist_mass, R_V, ext_coefs, N_fc, m_ini, cmpl_rnd, err_rnd):
    """
    Takes an isochrone and returns a synthetic cluster created according to
    a certain mass distribution.

    The synthetic cluster returned has the shape:

    synth_clust = [photometry, binary_idxs + extra_pars]

    photometry = [photom, errors]
    photom = [f1, f2, ..., fF, c1, c2, ..., cC]
    (where F and C are the total number of filters and colors defined)

    errors = [ef1, ef2, ..., efF, ec1, ec2, ..., ecC]
    (photometric errrors for each photometric dimension defined)

    Correct indexes of binary systems after completeness removal.
    binary_idxs = [i1, i2, ..., iN]

    Six lists containing the theoretical tracks extra parameters.
    extra_pars = [l1, l2, ..., l6]
    """

    # Unpack synthetic cluster parameters. The first two elements are the
    # metallicity and the age, which are already incorporated in the selected
    # isochrone.
    e, d, M_total, bin_frac = model_proper

    # import time
    # t1, t2, t3, t4, t5, t6, t7 = 0., 0., 0., 0., 0., 0., 0.

    # s = time.clock()
    # Move theoretical isochrone using the values 'e' and 'd'.
    isoch_moved = move_isochrone.main(isochrone, e, d, R_V, ext_coefs, N_fc)
    # t1 = time.clock() - s

    ##############################################################
    # # To generate a synthetic cluster with the full isochrone length,
    # # un-comment this line.
    # # This takes the max magnitude from the isochrone itself instead of using
    # # the input cluster file.
    # print "\nCluster's log(age): {:0.2f}".format(synth_cl_params[1])
    # print 'Fixed total mass: {:0.2f}'.format(M_total)
    # max_mag = max(isoch_moved[1]) + 0.5
    ##############################################################

    # s = time.clock()
    # Get isochrone minus those stars beyond the magnitude cut.
    isoch_cut = cut_max_mag.main(isoch_moved, max_mag_syn)
    # t2 = time.clock() - s

    # Empty list to pass if at some point no stars are left.
    synth_clust = []
    # Check for an empty array.
    if isoch_cut.any():

        # s = time.clock()
        # Mass distribution to produce a synthetic cluster based on
        # a given IMF and total mass.
        mass_dist = mass_distribution.main(st_dist_mass, M_total)
        # t3 = time.clock() - s

        # s = time.clock()
        # Interpolate masses in mass_dist into the isochrone rejecting those
        # masses that fall outside of the isochrone's mass range.
        isoch_mass = mass_interp.main(isoch_cut, mass_dist, m_ini)
        # t4 = time.clock() - s

        if isoch_mass.any():

            ##############################################################
            # # For plotting purposes: store a copy of this list before
            # # adding binaries since the list gets overwritten.
            # from copy import deepcopy
            # isoch_mass0 = deepcopy(isoch_mass)
            ##############################################################

            # s = time.clock()
            # Assignment of binarity.
            isoch_binar = binarity.main(isoch_mass, bin_frac, m_ini, N_fc)
            # t5 = time.clock() - s

            # s = time.clock()
            # Completeness limit removal of stars.
            isoch_compl = completeness_rm.main(
                isoch_binar, completeness, cmpl_rnd)
            # t6 = time.clock() - s

            ##############################################################
            # # Use when producing synthetic clusters from isochrones.
            # # Comment the line above.
            # isoch_compl = compl_func2(isoch_binar)
            ##############################################################

            if isoch_compl.any():

                # s = time.clock()
                # Get errors according to errors distribution.
                synth_clust = add_errors.main(
                    isoch_compl, err_lst, err_max, m_ini, err_rnd)
                # t7 = time.clock() - s

    ################################################################
    # # Plot synthetic cluster.
    # from synth_plot import synth_clust_plot as s_c_p
    # m, a = synth_cl_params[:2]
    # print(m, a, M_total)
    # out_name = str(m).split('.')[1] + '_' + str(a)
    # # out_name = 'synth_clust'
    # out_folder = '/home/gabriel/Descargas/'
    # path = out_folder + out_name + '.png'
    # s_c_p(N_fc, mass_dist, isochrone, synth_cl_params, isoch_moved,
    #       isoch_cut, isoch_mass0, isoch_binar, binar_idx0, isoch_compl,
    #       binar_idx, synth_clust, path)
    ################################################################

    # return np.array([t1, t2, t3, t4, t5, t6, t7])
    return synth_clust
