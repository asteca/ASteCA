
import numpy as np
from . import extin_coefs
from . import imf
from . import binarity


def main(
    td, synth_rand_seed, cmd_systs, filters, colors, IMF_name,
        all_syst_filters, N_interp, min_bmass_ratio, lkl_method, **kwargs):
    """
    Return list structured as:

    theor_tracks = [m1, m2, .., mN]
    mX = [age1, age2, ..., ageM]
    ageX = [f1,.., c1, c2,.., Mini, f1b,.., c1b, c2b,.., bp, Mb]
    where:
    fX  : individual filters (mags). Currently a single main mag is allowed.
    cX  : colors
    Mini: initial mass
    fXb : filters with binary data  ------
    cXb : colors with the binary data    | Only stored if binarity is set as a
    Mb  : binary masses -----------------| free param, or fixed to a value > 0.

    theor_tracks.shape = (Nz, Na, Nd, Ni)
    Nz: number of metallicities
    Na: number of log(age)s
    Nd: number of data columns
    Ni: number of interpolated values
    """
    # Set the random seed in numpy
    from . import set_rand_seed
    set_rand_seed.main(synth_rand_seed)

    # Extract some data
    all_met_vals, all_age_vals = td['fundam_params'][:2]
    all_masses, binar_fracs = td['fundam_params'][4], td['fundam_params'][5]
    Nmets, Max_mass = len(all_met_vals), all_masses[-1]

    # Obtain extinction coefficients.
    td['ext_coefs'] = extin_coefs.main(cmd_systs, filters, colors)

    # Set the binary flag
    binar_flag = False
    if len(binar_fracs) > 1 or binar_fracs[0] > 0.:
        binar_flag = True
    td['binar_flag'] = binar_flag

    # Store the number of defined filters and colors.
    td['N_fc'] = [len(filters), len(colors)]
    # Index of 'M_ini' (theoretical initial mass), stored in the
    # interpolated isochrones: right after the magnitude and color(s)
    td['m_ini_idx'] = td['N_fc'][0] + td['N_fc'][1]

    # Obtain mass distribution using the selected IMF.
    td['st_dist_mass'] = imf.main(IMF_name, Nmets, Max_mass)
    # tot_mem = 0.
    # for stm in st_dist_mass:
    #     tot_mem += stm[0].nbytes / 1024.**2 + stm[1].nbytes / 1024.**2

    # Combine all data into a single array of shape:
    # (N_z, N_age, N_data, N_interp), where 'N_data' depends on the number
    # of colors defines, and whether the binarity is on/off

    # Create the necessary colors, and position the magnitudes and colors
    # in the same order as they are read from the cluster's data file.
    interp_tracks, mags_cols_intp = interpIsochs(
        td['isoch_list'], td['extra_pars'], all_syst_filters,
        filters, colors, N_interp)

    if binar_flag:
        td['theor_tracks'], td['mean_bin_mr'] = binarity.binarGen(
            min_bmass_ratio, td['m_ini_idx'], td['N_fc'], interp_tracks,
            mags_cols_intp, all_met_vals, all_age_vals)
    else:
        td['theor_tracks'], td['mean_bin_mr'] = interp_tracks, 0.
    # Size of array
    # print("{:.0f} Mbs".format(theor_tracks.nbytes / 1024.**2))

    td['err_norm_rand'], td['binar_probs'], td['ext_unif_rand'] = randVals(
        lkl_method, td['st_dist_mass'], td['theor_tracks'])

    return td


def interpIsochs(
        isoch_list, extra_pars, all_syst_filters, filters, colors, N_interp):
    """
    Take the list of filters stored, create the necessary colors, arrange and
    interpolate all magnitudes and colors according to the order given to the
    photometric data read from file.

    The mags_cols_theor list contains the magnitudes used to create the
    defined colors. This is used to properly add binarity to the
    synthetic clusters.

    If more than one photometric system was used, then the 'extra_pars' array
    has extra 'M_ini' arrays (used to check) that need to be removed.
    TODO this will need the change when/if more extra parameters are
    stored beyond 'M_ini'
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

    # Store the index of each filter for each color read from data, as they
    # are stored in 'isoch_list'.
    fci = []
    for c in colors:
        ci = []
        for f in c[1].split(','):
            ci.append(all_filts.index(f))
        fci.append(ci)

    # Calculate the maximum number of stars in all isochrones, and add '10'
    # before interpolating. This is so that all tracks have the same number
    # of points.
    N_pts_max = 0
    for z in isoch_list:
        for a in z:
            N_pts_max = max(N_pts_max, len(a[0]))
    N_pts_max = N_pts_max + N_interp
    xx = np.linspace(0., 1., N_pts_max)

    interp_tracks = []
    for i, met in enumerate(isoch_list):
        m = []
        for j, age in enumerate(met):
            a = []
            for im in fi:
                xp = np.linspace(0, 1, len(age[im]))
                a.append(np.interp(xx, xp, age[im]))
            for ic in fci:
                col = np.array(age[ic[0]]) - np.array(age[ic[1]])
                xp = np.linspace(0, 1, len(col))
                a.append(np.interp(xx, xp, col))

            # Use the first 'M_ini' array, hence the '[0]'
            p = extra_pars[i][j][0]
            xp = np.linspace(0, 1, len(p))
            a.append(np.interp(xx, xp, p))
            m.append(a)
        interp_tracks.append(m)

    interp_tracks = np.array(interp_tracks)

    # Create list of theoretical colors, in the same orders as they are
    # read from the cluster's data file.
    # mags_cols_intp = [met1, met2, ..., metN]
    # metX = [age1, age2, ..., age_M]
    # ageX = [filter1, filter2, filter3, filter4, ..., filterQ]
    # such that: color1 = filter1 - filter2, color2 = filter3 - filter4, ...
    mags_cols_intp = []
    for met in isoch_list:
        m = []
        for age in met:
            a = []
            # For each color defined.
            for ic in fci:
                # For each filter of this color.
                xp = np.linspace(0, 1, len(age[ic[0]]))
                a.append(np.interp(xx, xp, age[ic[0]]))
                xp = np.linspace(0, 1, len(age[ic[1]]))
                a.append(np.interp(xx, xp, age[ic[1]]))
            m.append(a)
        mags_cols_intp.append(m)
    mags_cols_intp = np.array(mags_cols_intp)

    return interp_tracks, mags_cols_intp


def randVals(lkl_method, st_dist_mass, theor_tracks):
    """
    Generate lists of random values used by the synthetic cluster generating
    function.
    """
    # This is the maximum number of stars that will ever be interpolated into
    # an isochrone
    N_mass = 0
    for sdm in st_dist_mass:
        N_mass = max(len(sdm[0]), N_mass)

    # Number of metallicity values defined
    N_mets = theor_tracks.shape[0]

    err_norm_rand, binar_probs, ext_unif_rand = [], [], []
    for _ in range(N_mets):

        if lkl_method == 'tolstoy':
            # Tolstoy likelihood considers uncertainties, there's no need to
            # add it to the synthetic clusters.
            err_norm_rand.append(np.zeros(N_mass))
        else:
            err_norm_rand.append(np.random.normal(0., 1., N_mass))

        # Binary probabilities, one per metallicity
        binar_probs.append(np.random.uniform(0., 1., N_mass))

        # For the move_isoch() function. In place for #174
        ext_unif_rand.append(np.random.uniform(0., 1., N_mass))

    return err_norm_rand, binar_probs, ext_unif_rand
