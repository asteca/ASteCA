
import move_isochrone
import cut_max_mag
import mass_distribution
import mass_interp
import binarity
import completeness_rm
import add_errors


def main(err_max, err_lst, completeness, max_mag_syn, st_dist_mass,
         isochrone, R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd, synth_cl_params):
    """
    Takes an isochrone and returns a synthetic cluster created according to
    a certain mass distribution.

    The synthetic cluster returned has the shape:

    synth_clust = [photometry, binary_idxs, extra_pars]

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
    e, d, M_total, bin_frac = synth_cl_params[2:]

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
        # Index of m_ini (theoretical initial mass), stored in the theoretical
        # isochrones.
        m_ini = 2 * N_fc[0] + 2 * N_fc[1] + 2
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


if __name__ == "__main__":

    import pickle
    import numpy as np
    from synth_clust import imf
    import sys
    mpath = sys.path[0].replace('synth_clust', '').replace('packages/', '')

    # Data for RUP42 file (2.2 Mb) with BV, UB colors read
    print("Reading data")
    bin_mr = .7
    with open(mpath + 'theor_tracks.pickle', 'rb') as f:
        theor_tracks = pickle.load(f)
    err_max = .3
    err_lst = [
        np.array([4.82905373e-17, 1.60540044e+00, 3.04244535e-02]),
        np.array([2.91192261e-10, 8.99472478e-01, 4.53498125e-02]),
        np.array([5.83560740e-06, 4.79509209e-01, 3.82828812e-02])]
    completeness = [
        np.array([8.695, 8.96248, 9.22996, 9.49744, 9.76492, 10.0324,
                  10.29988, 10.56736, 10.83484, 11.10232, 11.3698, 11.63728,
                  11.90476, 12.17224, 12.43972, 12.7072, 12.97468, 13.24216,
                  13.50964, 13.77712, 14.0446, 14.31208, 14.57956, 14.84704,
                  15.11452, 15.382, 15.64948, 15.91696, 16.18444, 16.45192,
                  16.7194, 16.98688, 17.25436, 17.52184, 17.78932, 18.0568,
                  18.32428, 18.59176, 18.85924, 19.12672, 19.3942, 19.66168,
                  19.92916, 20.19664, 20.46412, 20.7316, 20.99908, 21.26656,
                  21.53404, 21.80152, 22.069]), 37,
        np.array([1., 0.98684211, 0.90131579, 0.82017544, 0.63048246,
                  0.36074561, 0.18311404, 0.10087719, 0.04166667, 0.01644737,
                  0.01754386, 0.00438596, 0.00109649])]
    max_mag_syn = 20.99
    ext_coefs = [
        (0.99974902186052628, -0.0046292182005786527),
        [(0.9999253931841896, 0.94553192328962365),
         (0.99974902186052628, -0.0046292182005786527)],
        [(0.96420188342499813, 1.784213363585738),
         (0.9999253931841896, 0.94553192328962365)]]
    N_fc = [1, 2]
    R_V = 3.1
    print("Data read.")

    Nm, Na = np.shape(theor_tracks)[0], np.shape(theor_tracks)[1]

    M_max_all = [1000., 5000., 10000., 25000., 50000., 100000., 250000.]
    print("Running")
    times_all = 0.
    for M_max in M_max_all:

        # set step for 100 masses
        M_step = M_max / 100.
        M_min = 50.
        masses = np.arange(M_min, M_max, M_step)
        m_sample_flag = False
        st_dist_mass = imf.main('kroupa_2002', 150., m_sample_flag, masses)

        cmpl_rnd = np.random.uniform(0., 1., 1000000)
        err_rnd = np.random.normal(0, 1, 1000000)

        N = 10000
        times = np.array([0., 0., 0., 0., 0., 0., 0.])
        for _ in range(N):
            ext = np.random.uniform(0., 2.)
            dm = np.random.uniform(8., 20.)
            M_total = np.random.choice(masses)
            bf = np.random.uniform(0., 1.)
            synth_cl_params = (np.nan, np.nan, ext, dm, M_total, bf)

            mi, ai = np.random.choice(Nm), np.random.choice(Na)
            isochrone = theor_tracks[mi][ai]

            times = times + main(
                err_max, err_lst, completeness, max_mag_syn, st_dist_mass,
                isochrone, R_V, ext_coefs, N_fc, cmpl_rnd, err_rnd,
                synth_cl_params)

        times_perc = np.round(100. * times / times.sum(), 1)
        print("    {:2.0f}-{:6.0f} {:7.2f}    {}".format(
            M_min, M_max, times.sum(), "    ".join(map(str, times_perc))))

        times_all += times.sum()
    print("    Total time: ~{:.3f} s --> ~xx% faster".format(times_all))
