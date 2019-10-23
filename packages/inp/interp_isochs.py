
import numpy as np
from ..synth_clust import binarity
from .. import update_progress


def main(
    mags_theor, cols_theor, mags_cols_theor, extra_pars, binar_fracs,
        bin_mr):
    """
    Interpolate extra points into all the filters, colors, filters of colors,
    and extra parameters (masses, etc). This allows the later IMF sampled
    masses to be more accurately interpolated into the theoretical
    isochrones.

    Return list structured as:

    theor_tracks = [m1, m2, .., mN]
    mX = [age1, age2, ..., ageM]
    ageX = [f1,.., c1, c2,.., f1b,.., c1b, c2b,.., bp, mb, m_ini,.., m_bol]
    where:
    fX:  individual filters (mags)
    cX:  colors
    fXb: filters with binary data
    cXb: colors with the binary data
    bp:  binary probabilities
    mb:  binary masses
    Mini,...: extra parameters.

    """

    N_mass_interp = interpPoints(mags_theor)
    print("Interpolating extra points ({}) into the isochrones".format(
        N_mass_interp))
    mags_intp, cols_intp, mags_cols_intp, extra_pars_intp = interp_isoch_data(
        (mags_theor, cols_theor, mags_cols_theor, extra_pars), N_mass_interp)

    # # Size of arrays in memory
    # sz = 0.
    # for arr in [mags_intp, cols_intp, mags_cols_intp, extra_pars_intp]:
    #     sz += np.array(arr).size * np.array(arr).itemsize
    # print("{:.3f} Mbs".format(sz / (1024.**2)))

    # The magnitudes for each defined color ('mags_cols_intp') are used here
    # and discarded after the colors (and magnitudes) with binarity assignment
    # are obtained.
    mags_binar, cols_binar, probs_binar, mass_binar = binarity.binarGen(
        binar_fracs, N_mass_interp, mags_intp, cols_intp, mags_cols_intp,
        extra_pars_intp, bin_mr)

    # Combine all data into a single array of shape:
    # (N_z, N_age, N_data, N_IMF_interp), where 'N_data' is the number of
    # sub-arrays in each array.
    theor_tracks = np.concatenate(
        (mags_intp, cols_intp, mags_binar, cols_binar, probs_binar, mass_binar,
            extra_pars_intp), axis=2)

    # DEPRECATED 02-10-2019
    #
    # This block is not needed anymore since there is no sorting anymore,
    # and the (z, a) final values can now be averaged from grid values.
    #
    # The mag sorting prevented us from being able to interpolate new (z, a)
    # values (when the MCMC or GA samples them) because the mass order is
    # not preserved. So, in order to be able to generate that interpolation
    # of values (which greatly improves the performance), we're back to the
    # 'old' cut_max_mag() function with no magnitude ordering.

    # # Sort all isochrones according to the main magnitude (min to max).
    # # This is necessary so that the cut_max_mag() function does not need
    # # to do this every time a new synthetic cluster is generated.
    # theor_tracks = [[[] for _ in a[0]] for _ in a]
    # for i, mx in enumerate(comb_data):
    #     for j, ax in enumerate(mx):
    #         # Ordering the data according to magnitude is necessary
    #         # if the cut_max_mag() function expects the data to be sorted
    #         # this way.
    #         theor_tracks[i][j] = ax[:, ax[0].argsort(kind='mergesort')]
    #
    # # The above sorting destroys the original order of the isochrones. This
    # # results in messy plots for the "best fit isochrone" at the end.
    # # (see: https://stackoverflow.com/q/35606712/1391441,
    # #       https://stackoverflow.com/q/37742358/1391441)
    # # To avoid this, we also store the not-interpolated, not sorted original
    # # values for the magnitudes and colors; just for the purpose of plotting
    # # the final isochrones.
    # plot_isoch_data = np.concatenate((mags_theor, cols_theor), axis=2)
    #
    # DEPRECATED 02-10-2019

    # # In place for #415
    # filename = 'temp.memmap'
    # sp = np.array(theor_tracks).shape
    # theor_tracks_mp = np.memmap(
    #     filename, dtype='float64', mode='w+', shape=sp)
    # theor_tracks_mp[:] = theor_tracks[:]

    return theor_tracks


def interpPoints(mags_theor):
    """
    Find the maximum number of points in all the read ages for all the
    metallicities.
    """
    N_pts_max = 0
    for z in mags_theor:
        for a in z:
            N_pts_max = max(N_pts_max, len(a[0]))
    # Fixed value. Will only change if some isochrone contains more than
    # this number of stars (not likely)
    N_mass_interp = 1500
    if N_mass_interp < N_pts_max:
        N_mass_interp = N_pts_max + 100

    return N_mass_interp


def interp_isoch_data(data_all, N):
    """
    Interpolate extra values for all the parameters in the theoretic
    isochrones.
    """
    all_interp_data = []
    for i, data in enumerate(data_all):
        interp_data = []
        # For each metallicity value list.
        for met in data:
            m = []
            # For each age value list.
            for age in met:
                a = []
                # For each filter/color/extra parameter in list.
                for fce in age:
                    t, xp = np.linspace(0., 1., N), np.linspace(0, 1, len(fce))
                    a.append(np.interp(t, xp, fce))
                m.append(a)
            interp_data.append(m)

        all_interp_data.append(interp_data)
        update_progress.updt(len(data_all), i + 1)

    return all_interp_data
