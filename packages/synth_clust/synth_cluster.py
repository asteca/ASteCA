
from . import zaWAverage
from . import move_isochrone
from . import cut_max_mag
from . import mass_distribution
from . import mass_interp
from . import binarity
from . import completeness_rm
from . import add_errors


def main(
    fundam_params, varIdxs, model, theor_tracks, err_max, err_lst,
    completeness, max_mag_syn, st_dist_mass, R_V, ext_coefs, N_fc, m_ini,
        cmpl_rnd, err_rnd):
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

    Lists containing the theoretical tracks extra parameters.
    extra_pars = [l1, l2, ..., l6]
    """

    # Generate a weighted average isochrone from the (z, log(age)) values in
    # the 'model'. Return proper values for fixed parameters.
    isochrone, model_proper = zaWAverage.main(
        theor_tracks, fundam_params, varIdxs, model)
    e, d, M_total, bin_frac = model_proper

    # Move theoretical isochrone using the values 'e' and 'd'.
    isoch_moved = move_isochrone.main(isochrone, e, d, R_V, ext_coefs, N_fc)

    ##############################################################
    # # To generate a synthetic cluster with the full isochrone length,
    # # un-comment this line.
    # # This takes the max magnitude from the isochrone itself instead of using
    # # the input cluster file.
    # print "\nCluster's log(age): {:0.2f}".format(synth_cl_params[1])
    # print 'Fixed total mass: {:0.2f}'.format(M_total)
    # max_mag = max(isoch_moved[1]) + 0.5
    ##############################################################

    # Get isochrone minus those stars beyond the magnitude cut.
    isoch_cut = cut_max_mag.main(isoch_moved, max_mag_syn)

    # Empty list to pass if at some point no stars are left.
    synth_clust = []
    # Check for an empty array.
    if isoch_cut.any():
        # Mass distribution to produce a synthetic cluster based on
        # a given IMF and total mass.
        mass_dist = mass_distribution.main(st_dist_mass, M_total)

        # Interpolate masses in mass_dist into the isochrone rejecting those
        # masses that fall outside of the isochrone's mass range.
        # This destroys the order by magnitude.
        isoch_mass = mass_interp.main(isoch_cut, mass_dist, m_ini)

        if isoch_mass.any():
            # Assignment of binarity.
            isoch_binar = binarity.main(isoch_mass, bin_frac, m_ini, N_fc)

            # Completeness limit removal of stars.
            isoch_compl = completeness_rm.main(
                isoch_binar, completeness, cmpl_rnd)

            if isoch_compl.any():
                # Get errors according to errors distribution.
                synth_clust = add_errors.main(
                    isoch_compl, err_lst, err_max, m_ini, err_rnd)

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

    return synth_clust
