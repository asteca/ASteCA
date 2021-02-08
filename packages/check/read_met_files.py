
import numpy as np
from packages.inp import readZA
from packages.inp import read_isochs


def check_get(pd):
    """
    Process all the metallicity files and the ages stored in them. To save
    time, we process and store all the theoretical isochrones data here.
    """

    # Only read files if best fit method is set to run. Else pass empty list.
    pd['m_ini_idx'] = -1

    if pd['best_fit_algor'] != 'n':
        # Print info about tracks.
        nt = '' if len(pd['all_syst_filters']) == 1 else 's'
        print("Processing {} theoretical isochrones\n"
              "in the photometric system{}:".format(pd['evol_track'], nt))
        for syst in pd['all_syst_filters']:
            print(" * {}".format(syst[0]))

        # Get all metallicity files and their values, and the log(age) values.
        met_files, met_vals_all, age_vals_all, ages_strs = readZA.main(**pd)

        # Store the common grid values for the metallicity and age.
        pd['fundam_params'][:2] = met_vals_all, age_vals_all

        # Get isochrones and their extra parameters (mass, etc.).
        isoch_list, extra_pars = read_isochs.main(
            met_files, ages_strs, pd['evol_track'], pd['CMD_extra_pars'],
            pd['all_syst_filters'])

        # Check equality of the initial mass across photometric systems.
        miniCheck(extra_pars, met_vals_all, age_vals_all)

        # Once the above check has passed, remove the extra 'M_ini' array
        # from 'extra_pars'.
        # TODO this will need the change when/if more extra parameters are
        # stored beyond 'M_ini'
        N_extra_pars = 1
        Nz, Na = len(met_vals_all), len(age_vals_all)
        extra_pars2 = [[[[] for _ in range(Na)] for _ in range(Nz)]
                       for _ in range(N_extra_pars)]
        for i, z in enumerate(extra_pars):
            for j, a in enumerate(z):
                extra_pars2[0][i][j] = np.array(a[0])
        extra_pars = extra_pars2

        # extra_pars2 = [[] for _ in met_vals_all]
        # for i, z in enumerate(extra_pars):
        #     masses = [[] for _ in age_vals_all]
        #     for j, a in enumerate(z):
        #         masses[j].append(a[0])
        #     extra_pars2[i] = masses
        # extra_pars = extra_pars2

        # Take the synthetic data from the unique filters read, create the
        # necessary colors, and position the magnitudes and colors in the
        # same order as they are read from the cluster's data file.
        # The mags_cols_theor list contains the magnitudes used to create the
        # defined colors. This is used to properly add binarity to the
        # synthetic clusters.
        mags_theor, cols_theor, mags_cols_theor = arrange_filters(
            isoch_list, pd['all_syst_filters'], pd['filters'], pd['colors'])

        pd['data_tracks'] = [
            mags_theor, cols_theor, extra_pars, mags_cols_theor]

        print("\nGrid values")
        print("z        : {:<5} [{}, {}]".format(
            len(met_vals_all), pd['fundam_params'][0][0],
            pd['fundam_params'][0][-1]))
        print("log(age) : {:<5} [{}, {}]".format(
            len(age_vals_all), pd['fundam_params'][1][0],
            pd['fundam_params'][1][-1]))

    return pd


def miniCheck(extra_pars, met_vals_all, age_vals_all):
    """
    The extra isochrone parameter 'M_ini' is assumed to be equal across
    photometric systems, for a given metallicity and age. We check here that
    this is the case.

    extra_pars.shape = (#met_vals, #log(ages), #phot_systs)
    """
    extra_pars = np.array(extra_pars, dtype=object)
    # If a single z and log(age) are defined, this array will have a shape
    # (#met_vals, #log(ages), #phot_systs, #stars). Hence the '[:3]'.
    Nz, Na, Ndim = extra_pars.shape[:3]
    if Ndim == 1:
        # Single photometric system defined. Nothing to check.
        return
    else:
        txt = "initial mass values are not equal across the\n" +\
            "photometric systems for the isochrone: z={}, log(age)={}"
        for d in range(1, Ndim):
            for z in range(Nz):
                for a in range(Na):
                    arr0, arrd = extra_pars[z, a, 0], extra_pars[z, a, d]
                    if not np.array_equal(arr0, arrd):
                        # Check across (z, log(age)) values for each
                        # photometric system.
                        raise ValueError(
                            txt.format(met_vals_all[z], age_vals_all[a]))


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
    # mags_theor = []
    # for met in isoch_list:
    #     m = []
    #     for age in met:
    #         a = []
    #         for i in fi:
    #             age_arr = np.array(age[i])
    #             a.append(age_arr)
    #         m.append(a)
    #     mags_theor.append(m)
    mags_theor = []
    for i in fi:
        m = []
        for met in isoch_list:
            a = []
            for age in met:
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
    # cols_theor = []
    # for met in isoch_list:
    #     m = []
    #     for age in met:
    #         a = []
    #         for ic in fci:
    #             # Generate color in the sense it was given in
    #             # 'params_input.dat'.
    #             col_arr = np.array(age[ic[0]]) - np.array(age[ic[1]])
    #             a.append(col_arr)
    #         m.append(a)
    #     cols_theor.append(m)
    cols_theor = []
    for ic in fci:
        m = []
        for met in isoch_list:
            a = []
            for age in met:
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
    # such that: color1 = filter1 - filter2, color2 = filter3 - filter4, ...
    # mags_cols_theor0 = []
    # for met in isoch_list:
    #     m = []
    #     for age in met:
    #         a = []
    #         # For each color defined.
    #         for ic in fci:
    #             # For each filter of this color.
    #             a.append(age[ic[0]])
    #             a.append(age[ic[1]])
    #         m.append(a)
    #     mags_cols_theor0.append(m)

    mags_cols_theor = []
    for ic in fci:
        for c in ic:
            m = []
            for met in isoch_list:
                a = []
                for age in met:
                    a.append(np.array(age[c]))
                m.append(a)
            mags_cols_theor.append(m)

    return mags_theor, cols_theor, mags_cols_theor
