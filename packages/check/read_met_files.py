
import numpy as np
from packages.inp import readZA
from packages.inp import read_isochs
from packages.inp import interp_isochs


def check_get(pd):
    """
    Process all the metallicity files and the ages stored in them. To save
    time, we process and store all the theoretical isochrones data here.
    """

    # Only read files if best fit method is set to run. Else pass empty list.
    pd['theor_tracks'] = []

    if pd['bf_flag']:
        # Print info about tracks.
        nt = '' if len(pd['all_syst_filters']) == 1 else 's'
        print("Processing {} theoretical isochrones\n"
              "in the photometric system{}:".format(
                  pd['all_evol_tracks'][pd['evol_track']][1], nt))
        for syst in pd['all_syst_filters']:
            print(" * {}".format(pd['cmd_systs'][syst[0]][0]))

        # Get all metallicity files and their values, and the log(age) values.
        met_files, met_vals_all, age_vals_all, ages_strs = readZA.main(**pd)

        # Store the common grid values for the metallicity and age.
        pd['fundam_params'][:2] = met_vals_all, age_vals_all

        # Get isochrones and their extra parameters (mass, etc.).
        isoch_list, extra_pars = read_isochs.main(
            met_files, ages_strs, pd['evol_track'], pd['CMD_extra_pars'],
            pd['all_syst_filters'])

        # Take the synthetic data from the unique filters read, create the
        # necessary colors, and position the magnitudes and colors in the
        # same order as they are read from the cluster's data file.
        # The mags_cols_theor list contains the magnitudes used to create the
        # defined colors. This is necessary to properly add binarity to the
        # synthetic clusters below.
        mags_theor, cols_theor, mags_cols_theor = arrange_filters(
            isoch_list, pd['all_syst_filters'], pd['filters'], pd['colors'])

        # Interpolate all the data in the isochrones (including the binarity
        # data)
        pd['theor_tracks'] = interp_isochs.main(
            mags_theor, cols_theor, mags_cols_theor, extra_pars,
            pd['fundam_params'][5], pd['bin_mr'])

        print("\nGrid values: {} [z], {} [log(age)]".format(
            len(met_vals_all), len(age_vals_all)))

    return pd


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
