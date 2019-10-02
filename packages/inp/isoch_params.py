
import sys
import numpy as np
from scipy.interpolate import interp1d
from . import read_isochs
from ..synth_clust import binarity
from .. import update_progress


def main(met_f_filter, age_values, cmd_evol_tracks, evol_track, bin_mr,
         all_syst_filters, cmd_systs, filters, colors, fundam_params,
         N_mass_interp, **kwargs):
    '''
    Read isochrones and parameters if best fit function is set to run.
    '''
    # Print info about tracks.
    nt = '' if len(all_syst_filters) == 1 else 's'
    print("Processing {} theoretical isochrones\n"
          "in the photometric system{}:".format(
              cmd_evol_tracks[evol_track][1], nt))
    for syst in all_syst_filters:
        print(" * {}".format(cmd_systs[syst[0]][0]))

    # Get isochrones and their extra parameters (mass, etc.).
    isoch_list, extra_pars = read_isochs.main(
        met_f_filter, age_values, cmd_evol_tracks, evol_track,
        all_syst_filters)

    # Take the synthetic data from the unique filters read, create the
    # necessary colors, and position the magnitudes and colors in the
    # same sense they are read from the cluster's data file.
    # The mags_cols_theor list contains the magnitudes used to create the
    # defined colors. This is necessary to properly add binarity to the
    # synthetic clusters below.
    mags_theor, cols_theor, mags_cols_theor = arrange_filters(
        isoch_list, all_syst_filters, filters, colors)

    # Interpolate extra points into all the filters, colors, filters of colors,
    # and extra parameters (masses, etc). This allows the later IMF sampled
    # masses to be more accurately interpolated into the theoretical
    # isochrones.

    # # Find the maximum number of points in all the read ages for all the
    # # metallicities.
    # N_pts_max, m_range = 0, 0.
    # for z in extra_pars:
    #     for j, m_ini in enumerate(z):
    #         m_min, m_max = m_ini[0][0], m_ini[0][-1]
    #         N_pts_max = max(N_pts_max, len(m_ini[0]))
    #         m_range = max(m_range, m_max - m_min)
    #         # print("{:.3f} {}   {:.2f}   {:.3f} {:.3f}".format(
    #         #     age_values[j], len(m_ini[0]), m_max-m_min, m_min, m_max))
    # print(N_pts_max, m_range)

    print("Interpolating masses into the isochrones")
    interp_data = interp_isoch_data(
        mags_theor, cols_theor, mags_cols_theor, extra_pars)

    # interp_data = []
    # for i, data in enumerate([
    #         mags_theor, cols_theor, mags_cols_theor, extra_pars]):
    #     interp_data.append(interp_isoch_data(data, N_mass_interp))
    #     update_progress.updt(4, i + 1)

    # a, b, c, d = mags_theor, cols_theor, mags_cols_theor, extra_pars
    a, b, c, d = interp_data

    # # Size of arrays in memory
    # sz = 0.
    # for arr in [a, b, c, d]:
    #     sz += np.array(arr).size * np.array(arr).itemsize
    # print("{:.3f} Mbs".format(sz / (1024.**2)))

    # The magnitudes for each defined color ('c') are used here and
    # discarded after the colors (and magnitudes) with binarity assignment
    # are obtained.
    mags_binar, cols_binar, probs_binar, mass_binar = binarity.binarGen(
        fundam_params[5], a, b, c, d, bin_mr)

    # Create list structured as:
    # theor_tracks = [m1, m2, .., mN]
    # mX = [age1, age2, ..., ageM]
    # ageX = [f1,.., c1, c2,.., f1b,.., c1b, c2b,.., bp, mb, m_ini,.., m_bol]
    # where:
    # fX:  individual filters (mags)
    # cX:  colors
    # fXb: filters with binary data
    # cXb: colors with the binary data
    # bp:  binary probabilities
    # mb:  binary masses
    # Mini,...: extra parameters.
    # theor_tracks = [[[] for _ in a[0]] for _ in a]
    theor_tracks = []
    for i, mx in enumerate(a):
        m = []
        for j, ax in enumerate(mx):
            m.append([ax, b[i][j], mags_binar[i][j], cols_binar[i][j],
                      probs_binar[i][j], mass_binar[i][j], d[i][j]])
        theor_tracks.append(m)

    dd = []
    for m in theor_tracks:
        m1 = []
        for a in m:
            m1.append(np.array(a)[:,0,:])
        dd.append(m1)

    theor_tracks = dd    

    # Combine all data into a single array of shape:
    # (N_z, N_age, N_data, N_IMF_interp), where 'N_data' is the number of
    # sub-arrays in each array.
    # comb_data = np.concatenate(
    #     (a, b, mags_binar, cols_binar, probs_binar, mass_binar, d), axis=2)

    # # Sort all isochrones according to the main magnitude (min to max).
    # # This is necessary so that the cut_max_mag() function does not need
    # # to do this every time a new synthetic cluster is generated.
    # theor_tracks = [[[] for _ in a[0]] for _ in a]
    # for i, mx in enumerate(comb_data):
    #     for j, ax in enumerate(mx):
    #         # TODO ordering the data according to magnitude is necessary
    #         # if the cut_max_mag() function expects the data to be sorted
    #         # this way. This however, prevents us from being able to
    #         # interpolate new (z, a) values when the MCMC samples them
    #         # because the mass order is not preserved anymore. So, in order
    #         # to be able to generate that interpolation of values (which
    #         # greatly improves the ptemcee performance), we're back to the
    #         # 'old' cut_max_mag() function and no magnitude ordering.
    #         # theor_tracks[i][j] = ax[:, ax[0].argsort(kind='mergesort')]
    #         theor_tracks[i][j] = ax

    # The above sorting destroys the original order of the isochrones. This
    # results in messy plots for the "best fit isochrone" at the end.
    # (see: https://stackoverflow.com/q/35606712/1391441,
    #       https://stackoverflow.com/q/37742358/1391441)
    # To avoid this, we also store the not-interpolated, not sorted original
    # values for the magnitudes and colors; just for the purpose of plotting
    # the final isochrones.
    plot_isoch_data = np.concatenate((mags_theor, cols_theor), axis=2)

    print("\nGrid values: {} [z], {} [log(age)]".format(
        len(fundam_params[0]), len(fundam_params[1])))

    # # In place for #415
    # filename = 'temp.memmap'
    # sp = np.array(theor_tracks).shape
    # theor_tracks_mp = np.memmap(filename, dtype='float64', mode='w+', shape=sp)
    # theor_tracks_mp[:] = theor_tracks[:]

    return theor_tracks, plot_isoch_data


def arrange_filters(isoch_list, all_syst_filters, filters, colors):
    """
    Take the list of filters stored, create the necessary colors, and arrange
    all magnitudes and colors according to the order given to the photometric
    data read from file.
    """
    # Extract names of all read filters in the order in which they are stored
    # in 'isoch_list'.
    all_filts = []
    for ps in all_syst_filters:
        all_filts = all_filts + list(ps[1:])

    # Store the index of each filter read from data, as they are stored in
    # 'isoch_list'.
    fi = []
    for f in filters:
        fi.append(all_filts.index(f[1]))
    # Create list of theoretical magnitudes, in the same orders as they are
    # read from the cluster's data file.
    mags_theor = []
    for met in isoch_list:
        m = []
        for age in met:
            a = []
            for i in fi:
                a.append(np.array(age[i]))
            m.append(a)
        mags_theor.append(m)

    # Store the index of each filter for each color read from data, as they
    # are stored in 'isoch_list'.
    fci = []
    for c in colors:
        ci = []
        for f in c[1].split(','):
            ci.append(all_filts.index(f))
        fci.append(ci)
    # Create list of theoretical colors, in the same orders as they are
    # read from the cluster's data file.
    cols_theor = []
    for met in isoch_list:
        m = []
        for age in met:
            a = []
            for ic in fci:
                # Generate color in the sense it was given in
                # 'params_input.dat'.
                a.append(np.array(age[ic[0]]) - np.array(age[ic[1]]))
            m.append(a)
        cols_theor.append(m)

    # Create list of theoretical colors, in the same orders as they are
    # read from the cluster's data file.
    # mags_cols_theor = [met1, met2, ..., metN]
    # metX = [age1, age2, ..., age_M]
    # ageX = [filter1, filter2, filter3, filter4, ..., filterQ]
    # such that: color1 = filter1 - filter2, color1 = filter3 - filter4, ...
    mags_cols_theor = []
    for met in isoch_list:
        m = []
        for age in met:
            a = []
            # For each color defined.
            for ic in fci:
                # For each filter of this color.
                a.append(age[ic[0]])
                a.append(age[ic[1]])
            m.append(a)
        mags_cols_theor.append(m)

    return mags_theor, cols_theor, mags_cols_theor


def interp_isoch_data(mags_theor, cols_theor, mags_cols_theor, extra_pars):
    """
    Interpolate extra values for all the parameters in the theoretic
    isochrones.
    """
    # mags_theor, cols_theor, mags_cols, extra_pars --> list
    # list = [m0, m1, .., mN] <-- N metallicities
    # mx = [a0, a1, .., aM]   <-- M log(ages)
    # ax = [[s1, s2, .., sQ]] <-- Q stars
    # list[i][j][0] --> i metallicity, j age

    mass_step = 0.01
    mass_vals1 = np.arange(.01, 5., mass_step)
    mass_step = 0.001
    mass_vals2 = np.arange(5., 18., mass_step)
    mass_step = 0.0005
    mass_vals3 = np.arange(18., 150., mass_step)
    mass_vals = np.array(
        mass_vals1.tolist() + mass_vals2.tolist() + mass_vals3.tolist())

    all_interp_data = []
    for data in (mags_theor, cols_theor, mags_cols_theor, extra_pars):
        interp_data = []
        # For each metallicity value.
        for i, met in enumerate(data):
            m = []
            # For each age value.
            for j, age in enumerate(met):
                a = []
                # For each filter/color/extra parameter in list.
                for fce in age:
                    # Initial masses for this  isochrone.
                    mass_ini = np.array(extra_pars[i][j][0])
                    # mass values must be in the range of the initial masses
                    mmin = np.searchsorted(mass_vals, mass_ini.min())
                    mmax = np.searchsorted(mass_vals, mass_ini.max())
                    mass_vals_cut = mass_vals[mmin:mmax]

                    y = np.array(data[i][j][0])
                    mfunc = interp1d(mass_ini, y)
                    ynew = mfunc(mass_vals_cut)
                    a.append(ynew)

                m.append(a)
            interp_data.append(m)
        all_interp_data.append(interp_data)

    # np.shape(all_interp_data): (4, Nz, Na)
    # 4  : ((mags, cols, mags_cols, extra_pars))
    # Nz : metallicity files (values)
    # Na : log(age) values per metallicity file

    return np.array(all_interp_data)
