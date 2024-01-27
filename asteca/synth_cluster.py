import numpy as np
from scipy.interpolate import make_interp_spline
from scipy import stats

#
from .synth_clust import imfs
from .synth_clust import zaWAverage
from .synth_clust import move_isochrone
from .synth_clust import cut_max_mag
# from .synth_clust import mass_distribution
from .synth_clust import mass_interp
from .synth_clust import extinction
from .synth_clust import binarity
# from .synth_clust import completeness_rm
from .synth_clust import add_errors


def synthcl_load(self):
    """ """
    isochs_dict = self.isochs.isochs_dict

    synthcl_dict = {"met_age_dict": isochs_dict["met_age_dict"]}

    N_colors = len(isochs_dict["color_filter_name"])
    if self.binarity is True:
        synthcl_dict["theor_tracks"] = add_binarity(
            self.gamma,
            isochs_dict["theor_tracks"],
            isochs_dict["mags_cols_intp"],
            N_colors,
            isochs_dict["color_filter_name"],
        )
    else:
        synthcl_dict["theor_tracks"] = isochs_dict["theor_tracks"]

    synthcl_dict["ext_coefs"] = extinction_coeffs(
        isochs_dict["mag_filter_name"],
        isochs_dict["color_filter_name"],
        isochs_dict["filter_lambdas"],
    )

    synthcl_dict["st_dist_mass"] = imf(
        synthcl_dict["theor_tracks"], self.IMF_name, self.max_mass
    )

    synthcl_dict["rand_norm_vals"], synthcl_dict["rand_unif_vals"] = randVals(
        synthcl_dict["theor_tracks"], synthcl_dict["st_dist_mass"]
    )

    return synthcl_dict


def add_binarity(gamma, theor_tracks, mags_cols_intp, N_colors, color_filter_name):
    """
    For each theoretical isochrone defined.
        1. Draw random secondary masses for *all* stars.
        2. Interpolate magnitude values for these masses
        3. Calculate the magnitudes for the binary system.
        4. Calculate the colors for the binary system.
    """
    #
    print("Generating binary data")

    mag_idx = 0
    m_ini_idx = mag_idx + N_colors + 1

    # Extend to accommodate binary data
    Nz, Na, Nd, Ni = theor_tracks.shape
    theor_tracks = np.concatenate((theor_tracks, np.zeros([Nz, Na, Nd, Ni])), axis=2)

    # For each metallicity defined.
    for mx, met in enumerate(theor_tracks):
        # For each age defined.
        for ax, isoch in enumerate(met):
            # Extract initial masses for this isochrone.
            mass_ini = isoch[m_ini_idx]

            # Mass-ratio distribution
            mass_ratios = qDistribution(mass_ini, gamma)

            # Calculate random secondary masses of these binary stars
            # between bin_mass_ratio*m1 and m1, where m1 is the primary
            # mass.
            m2 = mass_ratios * mass_ini

            # Calculate unresolved binary magnitude
            mag = isoch[mag_idx]
            mag_m2 = np.interp(m2, mass_ini, mag)
            theor_tracks[mx][ax][m_ini_idx + 1] = mag_combine(mag, mag_m2)

            # Calculate unresolved color for each color defined.
            for ic, color in enumerate(color_filter_name):
                mc1 = mags_cols_intp[mx][ax][color[0]]
                mc2 = mags_cols_intp[mx][ax][color[1]]
                f_m1 = np.interp(m2, mass_ini, mc1)
                f1 = mag_combine(mc1, f_m1)
                f_m2 = np.interp(m2, mass_ini, mc2)
                f2 = mag_combine(mc2, f_m2)
                theor_tracks[mx][ax][m_ini_idx + 1 + 1 + ic] = f1 - f2

            # Secondary masses
            theor_tracks[mx][ax][-1] = m2

    return theor_tracks


def qDistribution(M1, gamma):
    """
    D&K      : Distribution of mass-ratios versus primary masses
    (Duchene & Kraus 2013). Mass dependent.
    powerlaw : Power-law distribution with shape parameter 'gamma'. Not mass
    dependent.


    Use 'gamma + 1' in the power-law distribution below because in D&K this
    distribution is defined as f(q)~q^gamma, while numpy's distribution is
    defined as a*x^(a-1).
    """
    N = M1.size

    try:
        gamma = float(gamma)
        mass_ratios = np.random.power(gamma + 1, N)
    except ValueError:

        if gamma == 'D&K':
            msk1, gamma1 = M1 <= 0.1, 4.2
            msk2, gamma2 = (M1 > 0.1) & (M1 <= 0.6), 0.4
            msk3, gamma3 = (M1 > 0.6) & (M1 <= 1.4), 0.3
            msk4, gamma4 = (M1 > 1.4) & (M1 <= 6.5), -0.5
            msk5, gamma5 = (M1 > 6.5) & (M1 <= 16), 0.0  # <- Not sure. Use uniform
            msk6, gamma6 = M1 > 16, 0.0  # <- Not sure. Use uniform

            mass_ratios = np.zeros(N)
            for msk, gammaX in (
                (msk1, gamma1),
                (msk2, gamma2),
                (msk3, gamma3),
                (msk4, gamma4),
                (msk5, gamma5),
                (msk6, gamma6),
            ):
                q = np.random.power(gammaX + 1, msk.sum())
                mass_ratios[msk] = q

        # Fisher's distribution
        elif gamma == 'fisher_stepped':
            # Fisher, Schröder & Smith (2005), 10.1111/j.1365-2966.2005.09193.x;
            # Table 3, stepped
            xk = np.linspace(0., 1., 10)
            pk = np.array([29., 29., 30., 32., 31., 32., 36., 45., 27., 76.])

        elif gamma == 'fisher_peaked':
            # Fisher, Schröder & Smith (2005), 10.1111/j.1365-2966.2005.09193.x;
            # Table 3, peaked
            xk = np.linspace(0., 1., 10)
            pk = np.array([27., 30., 34., 33., 29., 26., 27., 33., 41., 89.])

        elif gamma == 'raghavan':
            # Raghavan et al. (2010), 10.1088/0067-0049/190/1/1; Fig 16 (left)
            xk = np.linspace(0., 1., 20)
            pk = np.array([
                0.53,  2.61,  0.53,  4.67,  7.81,  3.64,  9.89,  5.71,  4.69,
                5.73,  4.67,  6.76,  5.75,  5.73,  2.61,  5.71,  4.72,  5.71,
                3.64, 12.99])

        if gamma != 'D&K':
            def fQ(xk, pk):
                """
                Discrete function
                """
                pk /= pk.sum()
                fq = stats.rv_discrete(a=0., b=1., values=(xk, pk))
                return fq
            fq = fQ(xk, pk)
            # 'ppf' is the inverse CDF
            mass_ratios = fq.ppf(np.random.uniform(0., 1., N))

    return mass_ratios


def mag_combine(m1, m2):
    """
    Combine two magnitudes. This is a faster re-ordering of the standard
    formula:

    -2.5 * np.log10(10 ** (-0.4 * m1) + 10 ** (-0.4 * m2))

    """
    c = 10**-0.4
    mbin = -2.5 * (-0.4 * m1 + np.log10(1.0 + c ** (m2 - m1)))
    return mbin


def extinction_coeffs(mag_filter_name, color_filter_name, filter_lambdas):
    """
    Obtain extinction coefficients for all the observed filters and colors,
    in the order in which they are stored in theor_tracks.

    ext_coefs = [ec_mag, ec_col1, ...]
    """
    print("Obtaining extinction coefficients")

    # For the magnitude.
    # Effective wavelength in Armstrong.
    eff_wave = filter_lambdas[mag_filter_name]
    # CCM coefficient. Effective wavelength in inverse microns.
    ext_coefs = [ccm_model(10000.0 / eff_wave)]

    # For colors.
    for color in color_filter_name:
        c_filt1, c_filt2 = color
        # Effective wavelength in Armstrong.
        eff_wave1, eff_wave2 = filter_lambdas[c_filt1], filter_lambdas[c_filt2]
        ext_coefs.append(
            [ccm_model(10000.0 / eff_wave1), ccm_model(10000.0 / eff_wave2)]
        )

    return ext_coefs


def ccm_model(mw):
    """
    Cardelli, Clayton, and Mathis (1989 ApJ. 345, 245) model for extinction
    coefficients with updated coefficients for near-UV from O'Donnell
    (1994, ApJ, 422, 158).

    ccm_coef = a + b / Rv

    Implementation taken from:

    http://idlastro.gsfc.nasa.gov/ftp/pro/astro/ccm_unred.pro

    There appears to be an error in the Far-UV range in the original IDL
    routine where the maximum inverse wavelength is 11 and it should be 10
    according to Cardelli et al. 1989 (pag 251, Eq (5,a,b)).
    """

    if 0.3 <= mw < 1.1:
        # Infrared.
        a, b = 0.574 * (mw**1.61), -0.527 * (mw**1.61)

    elif 1.1 <= mw < 3.3:
        # Optical/NIR.
        # Original coefficients from CCM89
        # c1 = [1., 0.17699, -0.50447, -0.02427, 0.72085, 0.01979, -0.77530,
        #       0.32999]
        # c2 = [0., 1.41338, 2.28305, 1.07233, -5.38434, -0.62251, 5.30260,
        #       -2.09002]
        # New coefficients from O'Donnell (1994)
        c1 = [1.0, 0.104, -0.609, 0.701, 1.137, -1.718, -0.827, 1.647, -0.505]
        c2 = [0.0, 1.952, 2.908, -3.989, -7.985, 11.102, 5.491, -10.805, 3.347]
        y = mw - 1.82
        # Reverse because polyval starts from the highest degree.
        c1.reverse(), c2.reverse()
        a, b = np.polyval(c1, y), np.polyval(c2, y)

    elif 3.3 <= mw < 8.0:
        # Mid-UV
        F_a, F_b = 0.0, 0.0
        if mw >= 5.9:
            y = mw - 5.9
            F_a = -0.04473 * y**2 - 0.009779 * y**3
            F_b = 0.2130 * y**2 + 0.1207 * y**3
        a = 1.752 - 0.316 * mw - (0.104 / ((mw - 4.67) ** 2 + 0.341)) + F_a
        b = -3.090 + 1.825 * mw + (1.206 / ((mw - 4.62) ** 2 + 0.263)) + F_b

    elif 8.0 <= mw <= 10.0:
        # Far-UV
        c1 = [-1.073, -0.628, 0.137, -0.070]
        c2 = [13.670, 4.257, -0.420, 0.374]
        y = mw - 8.0
        c1.reverse(), c2.reverse()
        a, b = np.polyval(c1, y), np.polyval(c2, y)
    else:
        raise ValueError(
            "The effective wavelength is {} [1/micron], beyond "
            "the CCM model limit (10 [1/micron]).".format(mw)
        )

    return a, b


def imf(theor_tracks, IMF_name, max_mass):
    """
    Returns the number of stars per interval of mass for the selected IMF.

    Parameters
    ----------
    IMF_name : str
      Name of the IMF to be used.
    max_mass: float
      Maximum mass defined.

    Returns
    -------
    st_dist_mass : list
      Tuple that contains: a given number of stars sampled from the selected
      IMF, that (approximately) sum to the associated total mass; the
      cumulative sum of those masses.
      One list per metallicity value defined is returned. This is to add some
      variety to the sampled IMF. Ideally this should be done every time a
      new isochrone is populated, but this is not possible due to performance
      costs and reproducibility.

    """
    print(f"Sampling selected IMF ({IMF_name})")
    inv_cdf = invTrnsfSmpl(IMF_name)

    # Number of metallicity values
    Nmets, Nages = theor_tracks.shape[:2]

    # m_ini_idx = N_colors + 1
    # m_lowest = get_lowest_mass(theor_tracks, max_obs_mag, max_dm, N_colors)

    st_dist_mass = []
    for _ in range(Nmets):
        met_lst = []
        for _ in range(Nages):
            sampled_IMF = sampleInv(max_mass, inv_cdf)
            # st_dist_mass += [[sampled_IMF, np.cumsum(sampled_IMF)]]
            met_lst.append(sampled_IMF)
        st_dist_mass.append(met_lst)

        # sort_m_idxs = []
        # for j in range(Nages):

        #     mass_ini = theor_tracks[i][j][m_ini_idx]
        #     mmin, mmax = mass_ini.min(), mass_ini.max()
        #     msk_m = (sampled_IMF >= mmin) & (sampled_IMF <= mmax)
        #     sampled_IMF = sampled_IMF[msk_m]

        #     sort_m_idxs.append([sampled_IMF, np.cumsum(sampled_IMF), np.searchsorted(mass_ini, sampled_IMF)])
        # st_dist_mass.append(sort_m_idxs)

        # # In place for #545
        # # print("\n\n----> BE CAREFUL HERE IN THE IMF FUNCTION <----\n\n")
        # sampled_IMF, mass_removed = remove_low_mass(sampled_IMF, m_lowest)
        # st_dist_mass += [[sampled_IMF, mass_removed + np.cumsum(sampled_IMF)]]

    return st_dist_mass


# def get_lowest_mass(theor_tracks, max_obs_mag, max_dm, N_colors):
#     """ """
#     m_ini_idx = N_colors + 1
#     for met in theor_tracks:
#         for age in met:
#             mag = age[0]
#             i = np.argmin(abs(max_obs_mag - (mag + max_dm + 1).max()))
#             breakpoint()
#             print((mag+max_dm).max())
#             # mag_max = mag + max_dm
#     breakpoint()
#     return m_lowest


# def remove_low_mass(sampled_IMF, m_lowest):
#     """
#     Remove the lowest mass stars from the array, to save memory
#     """
#     msk = sampled_IMF > m_lowest
#     mass_removed = sampled_IMF[~msk].sum()
#     return sampled_IMF[msk], mass_removed


def invTrnsfSmpl(IMF_name, m_low=0.08, m_high=150, mass_step=0.05):
    """
    IMF inverse transform sampling.

    Asked here: https://stackoverflow.com/q/21100716/1391441
    """
    # IMF mass interpolation step and grid values.
    mass_values = np.arange(m_low, m_high, mass_step)

    def IMF_func(m_star, IMF_name):
        return imfs.get_imf(IMF_name, m_star)

    # CDF_samples = []
    # for m in mass_values:
    #     x = quad(IMF_func, m_low, m, args=(IMF_name))[0]
    #     CDF_samples.append(x)
    # # Normalize values
    # CDF_samples = np.array(CDF_samples) / max(CDF_samples)
    # # These are (0, 1)
    # # CDF_min, CDF_max = CDF_samples.min(), CDF_samples.max()
    # # Inverse CDF
    # from scipy.interpolate import interp1d
    # inv_cdf = interp1d(CDF_samples, mass_values)  # Almost deprecated

    # x_old, CDF_samples = -np.inf, []
    # for m in mass_values:
    #     x = quad(IMF_func, m_low, m, args=(IMF_name))[0]
    #     # This ensures that CDF_samples is monotonically increasing, bypassing
    #     # rounding errors
    #     if x <= x_old:
    #         x = x_old + 1e-15
    #     CDF_samples.append(x)
    #     x_old = x

    # # Normalize values
    # CDF_samples = np.array(CDF_samples) / max(CDF_samples)
    # # These are (0, 1)
    # # CDF_min, CDF_max = CDF_samples.min(), CDF_samples.max()

    # # k=1 is important otherwise the interpolator can become unstable for values
    # # close to 1 and return huge masses and even negative ones
    # inv_cdf = make_interp_spline(CDF_samples, mass_values, k=1)

    # The lower mass region needs to be sampled more accurately
    mass_values = list(np.linspace(m_low, .5, 50))
    mass_values += list(np.linspace(.501, 10, 250))
    mass_values += list(np.linspace(10.01, m_high, 250))

    CDF_samples = []
    IMF_old, m_old, area_CDF = IMF_func(m_low, IMF_name), m_low, 0.
    for m in mass_values[1:]:
        # Approximate integral with rectangular area, and add to previous total area
        IMF_new = IMF_func(m, IMF_name)
        area_CDF += .5*(IMF_new + IMF_old) * (m - m_old)
        CDF_samples.append(area_CDF)
        IMF_old, m_old = IMF_new, m
    CDF_samples = np.array(CDF_samples) / max(CDF_samples)

    # k=1 is important otherwise the interpolator can become unstable for values
    # close to 1 and return huge masses and even negative ones
    inv_cdf = make_interp_spline(CDF_samples, mass_values[1:], k=1)

    return inv_cdf


def sampleInv(Max_mass, inv_cdf):
    """
    Sample the inverse CDF up to `Max_mass`
    """
    def sampled_inv_cdf(N):
        mr = np.random.rand(N)
        return inv_cdf(mr)

    # Sample in chunks until the maximum defined mass is reached.
    N_chunk = max(100, int(2.5 * Max_mass / 100.0))
    Mass_tot, mass_samples = 0, []
    while Mass_tot < Max_mass:
        masses = sampled_inv_cdf(N_chunk).tolist()
        mass_samples += masses
        Mass_tot += sum(masses)
    sampled_IMF = np.array(mass_samples)

    return sampled_IMF


def randVals(theor_tracks, st_dist_mass):
    """
    Generate lists of random values used by the synthetic cluster generating
    function.
    """
    # This is the maximum number of stars that will ever be interpolated into
    # an isochrone
    N_isoch, N_mass = theor_tracks.shape[-1], 0
    for sdm in st_dist_mass:
        N_mass = max(len(sdm[0]), N_mass, N_isoch)

    # Used by `move_isochrone()` and `add_errors`
    rand_norm_vals = np.random.normal(0.0, 1.0, (2, N_mass))

    # Used by `move_isochrone()`, `binarity()`, `completeness_rm()`
    rand_unif_vals = np.random.uniform(0.0, 1.0, (3, N_mass))

    return rand_norm_vals, rand_unif_vals


def synthcl_generate(
    self, model_fit, model_fixed, cluster_dict, full_synth_arr=False
):
    """
    Takes an isochrone and returns a synthetic cluster created according to
    a certain mass distribution.

    The synthetic cluster returned has the shape:

    * full_synth_arr = False

    synth_clust = [mag, c1, (c2)]

    where c1 and c2 colors defined.

    * full_synth_arr = True

    synth_clust = [mag, c1, (c2), m_ini_1, mag_b, c1_b, (c2_b), m_ini_2]

    where 'm_ini_1, m_ini_2' are the primary and secondary masses of the
    binary systems. The single systems only have a '0' stored in 'm_ini_2'.
    """
    synthcl_dict = self.synthcl_dict
    m_ini_idx = cluster_dict["N_colors"] + 1

    # Return proper values for fixed parameters and parameters required
    # for the (z, log(age)) isochrone averaging.
    met, loga, beta, av, dr, rv, dm, ml, mh, al, ah = properModel(
        synthcl_dict["met_age_dict"], model_fixed, model_fit
    )

    # Generate a weighted average isochrone from the (z, log(age)) values in
    # the 'model'.
    isochrone = zaWAverage.main(
        synthcl_dict["theor_tracks"],
        synthcl_dict["met_age_dict"],
        m_ini_idx,
        met,
        loga,
        ml,
        mh,
        al,
        ah,
    )

    # Move theoretical isochrone using the distance modulus
    isoch_moved = move_isochrone.main(isochrone, m_ini_idx, dm)

    # Empty list to pass if at some point no stars are left.
    synth_clust = np.array([])

    # Remove isochrone stars beyond the maximum magnitude
    isoch_cut, _ = cut_max_mag.main(isoch_moved, cluster_dict["max_mag_syn"])
    if not isoch_cut.any():
        return synth_clust

    # Interpolate IMF's sampled masses into the isochrone.
    isoch_mass = mass_interp.main(
        isoch_cut, m_ini_idx, synthcl_dict["st_dist_mass"][ml][al])
    if not isoch_mass.any():
        return synth_clust

    # Assignment of binarity.
    isoch_binar = binarity.main(
        self.binarity, self.alpha, beta, m_ini_idx,
        synthcl_dict["rand_unif_vals"][1], isoch_mass
    )

    # Apply extinction correction
    isoch_extin = extinction.main(
        isoch_binar,
        av,
        dr,
        rv,
        synthcl_dict["ext_coefs"],
        self.DR_dist,
        self.DR_percentage,
        synthcl_dict["rand_norm_vals"][0],
        synthcl_dict["rand_unif_vals"][0],
        m_ini_idx,
    )

    # Remove stars moved beyond the maximum magnitude
    isoch_extin, N_phot2 = cut_max_mag.main(isoch_extin, cluster_dict["max_mag_syn"])
    if not isoch_extin.any():
        return synth_clust

    # # Completeness removal of stars.
    # isoch_compl, N_phot3 = completeness_rm.main(
    #     isoch_extin, self.completeness_f, synthcl_dict["rand_unif_vals"][2]
    # )
    # if not isoch_compl.any():
    #     return synth_clust

    # Keep the same number of synthetic stars as observed stars
    isoch_extin = isoch_extin[:, : cluster_dict["N_obs_stars"]]

    # Assign errors according to errors distribution.
    synth_clust = add_errors.main(
        isoch_extin, cluster_dict["err_lst"], synthcl_dict["rand_norm_vals"][1]
    )

    if full_synth_arr:
        return synth_clust
    return synth_clust[:m_ini_idx]


def properModel(met_age_dict, model_fixed, model_fit):
    """
    Define the 'proper' model with values for (z, a) taken from its grid,
    and filled values for those parameters that are fixed.

    Parameters
    ----------
    model : array
      Array of *free* fundamental parameters only (ie: in varIdxs).

    Returns
    -------
    model_fixed  : array
      All fundamental parameters, including the fixed parameters that are
      missing from 'model'.
    ml, mh, al, ah : ints
      Indexes of the (z, a) values in the grid that define the box that enclose
      the proper (z, a) values.

    """
    model_comb = model_fixed | model_fit
    model_full = np.array(
        [model_comb[k] for k in ["z", "loga", "beta", "Av", "DR", "Rv", "dm"]]
    )

    if len(met_age_dict["z"]) == 1:
        ml = mh = 0
    else:
        par = met_age_dict["z"]
        mh = min(len(par) - 1, np.searchsorted(par, model_full[0]))
        ml = mh - 1

    if len(met_age_dict["a"]) == 1:
        al = ah = 0
    else:
        par = met_age_dict["a"]
        ah = min(len(par) - 1, np.searchsorted(par, model_full[1]))
        al = ah - 1

    return *model_full, ml, mh, al, ah
