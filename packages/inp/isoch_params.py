
import numpy as np
import read_isochs


def arrange_filters(isoch_list, all_syst_filters, filters, colors, **kwargs):
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
    mags_the = []
    for met in isoch_list:
        m = []
        for age in met:
            a = []
            for i in fi:
                a.append(age[i])
            m.append(a)
        mags_the.append(m)

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
    cols_the = []
    for met in isoch_list:
        m = []
        for age in met:
            a = []
            for ic in fci:
                # Generate color in the sense it was given in
                # 'params_input.dat'.
                a.append(np.array(age[ic[0]]) - np.array(age[ic[1]]))
            m.append(a)
        cols_the.append(m)

    return mags_the, cols_the


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
            for f in age:
                t, xp = np.linspace(0, 1, N), np.linspace(0, 1, len(f))
                a.append(np.interp(t, xp, f))
            m.append(a)
        interp_data.append(m)

    return interp_data


def main(pd, met_f_filter, age_values):
    '''
    Read isochrones and parameters if best fit function is set to run.
    '''
    # Print info about tracks.
    print("Processing {} theoretical isochrones".format(
        pd['cmd_evol_tracks'][pd['evol_track']][1]))

    for syst in pd['all_syst_filters']:
        print("in the '{}' photometric system.\n".format(
            pd['cmd_systs'][syst[0]][0]))

    # Get isochrones and their extra parameters (mass, etc.).
    isoch_list, extra_pars = read_isochs.main(met_f_filter, age_values,
                                              **pd)

    # Take the synthetic data from the unique filters read, create the
    # necessary colors, and position the magnitudes and colors in the
    # same sense they are read from the cluster's data file.
    mags_the, cols_the = arrange_filters(isoch_list, **pd)

    # Interpolate extra points into all the isochrones.
    pd['mags_interp'] = interp_isoch_data(mags_the)
    pd['cols_interp'] = interp_isoch_data(cols_the)
    pd['extra_pars_interp'] = interp_isoch_data(extra_pars)

    # Obtain number of models in the solutions space.
    lens = [len(_) for _ in pd['param_values']]
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

    return pd
