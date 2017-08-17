
import numpy as np
import read_isochs


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
                a.append(age[i])
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


def interp_isoch_data(data, N=2000):
    '''
    Interpolate extra values for all the parameters in the theoretic
    isochrones.
    '''
    interp_data = []
    # For each metallicity value list.
    for met in data:
        m = []
        # For each age value list.
        for age in met:
            a = []
            # For each filter/color/extra parameter in list.
            for fce in age:
                t, xp = np.linspace(0, 1, N), np.linspace(0, 1, len(fce))
                a.append(np.interp(t, xp, fce))
            m.append(a)
        interp_data.append(m)

    return interp_data


def main(met_f_filter, age_values, cmd_evol_tracks, evol_track,
         all_syst_filters, cmd_systs, filters, colors, fundam_params,
         **kwargs):
    '''
    Read isochrones and parameters if best fit function is set to run.
    '''
    # Print info about tracks.
    print("Processing {} theoretical isochrones".format(
        cmd_evol_tracks[evol_track][1]))

    for syst in all_syst_filters:
        print("in the '{}' photometric system.\n".format(
            cmd_systs[syst[0]][0]))
    fs = ', '.join(_[1] for _ in filters)
    cs = ', '.join('(' + _[1].replace(',', '-') + ')' for _ in colors)
    print("Filter: {}".format(fs))
    print("Color:  {}\n".format(cs))

    # Get isochrones and their extra parameters (mass, etc.).
    isoch_list, extra_pars = read_isochs.main(met_f_filter, age_values,
                                              evol_track, all_syst_filters)

    # Take the synthetic data from the unique filters read, create the
    # necessary colors, and position the magnitudes and colors in the
    # same sense they are read from the cluster's data file.
    # The mags_cols_theor list contains the magnitudes used to create the
    # defined colors. This is necessary to properly add binarity to the
    # synthetic clusters later on.
    mags_theor, cols_theor, mags_cols_theor = arrange_filters(
        isoch_list, all_syst_filters, filters, colors)

    # Interpolate extra points into all the filters, colors, filters of colors,
    # and extra parameters (masses, etc)
    a = interp_isoch_data(mags_theor)
    b = interp_isoch_data(cols_theor)
    c = interp_isoch_data(mags_cols_theor)
    d = interp_isoch_data(extra_pars)
    # Create list structured as:
    # theor_tracks = [m1, m2, .., mN]
    # mX = [age1, age2, ..., ageM]
    # ageX = [f1, f2, ..., c1, c2, ..., fc1, fc2, ..., m_ini, .., m_bol]
    # where fX are the individual filters (mags), cX are the colors, fcX are
    # the filters that make up the colors (where c1=(fc1-fc2), c2=(fc3-fc4)),
    # and the final lists are the six extra parameters.
    # Create empty lists for each metallicity, and empty sublists for each age.
    theor_tracks = [[[] for _ in a[0]] for _ in a]
    for l in [a, b, c, d]:
        for i, mx in enumerate(l):
            for j, ax in enumerate(mx):
                theor_tracks[i][j] = theor_tracks[i][j] + ax

    # Obtain number of models in the solutions space.
    lens = [len(_) for _ in fundam_params]
    total = reduce(lambda x, y: x * y, lens, 1)
    print(
        "Number of values per parameter:\n"
        "  {} metallicity values (z),\n"
        "  {} age values (per z),\n"
        "  {} reddening values,\n"
        "  {} distance values,\n"
        "  {} mass values,\n"
        "  {} binary fraction values.".format(*lens))
    print("  = {:.1e} approx total models.\n".format(total))

    return theor_tracks
