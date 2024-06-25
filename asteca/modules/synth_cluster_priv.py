import numpy as np
from scipy import stats
from .imfs import invTrnsfSmpl, sampleInv


def error_distribution(mag, e_mag, e_colors, rand_norm_vals):
    """
    Extract the magnitude and color(s) uncertainties to use as error model for the
    synthetic clusters.
    """

    def filnans(data):
        msk = np.isnan(data)
        data[msk] = np.interp(np.flatnonzero(msk), np.flatnonzero(~msk), data[~msk])
        return data

    # Replace nan values with interpolated values
    mag = filnans(mag)
    e_mag = filnans(e_mag)
    e_colors = [filnans(_) for _ in e_colors]

    # The minus generates reversed sorting in magnitude. This is important so that
    # synthetic clusters with a smaller magnitude range are assigned uncertainties
    # from the bottom up (i.e.: from the largest to the smallest magnitudes)
    idx = np.argsort(-mag)

    # Multiplying by a random normal float centered at 0 with STDDEV=1 generates the
    # uncertainty values ready to be added to the synthetic photometry.
    N = len(mag)
    err_dist = [rand_norm_vals[:N] * e_mag[idx]]
    for e_col in e_colors:
        err_dist.append(rand_norm_vals[:N] * e_col[idx])

    return err_dist


def add_binarity(self) -> np.ndarray:
    """For each theoretical isochrone defined.
        1. Draw random secondary masses for *all* stars.
        2. Interpolate magnitude values for these masses
        3. Calculate the magnitudes for the binary system.
        4. Calculate the colors for the binary system.

    theor_tracks.shape = (Nz, Na, N_cols, N_interp)

    If binarity is processed:
        N_cols: magnitude + color(s) + initial masses + unresolved mag +
            unresolved color(s) + secondary masses
    else:
        N_cols: magnitude + colors + initial mass

    """
    all_colors = [self.isochs.color]
    if self.isochs.color2 is not None:
        all_colors.append(self.isochs.color2)
    N_colors = len(all_colors)

    mag_idx = 0
    m_ini_idx = mag_idx + N_colors + 1

    # Extend to accommodate binary data
    Nz, Na, Nd, Ni = self.isochs.theor_tracks.shape
    theor_tracks = np.concatenate(
        (self.isochs.theor_tracks, np.zeros([Nz, Na, Nd, Ni])), axis=2
    )

    # For each metallicity defined.
    for mx, met in enumerate(theor_tracks):
        # For each age defined.
        for ax, isoch in enumerate(met):
            # Extract initial masses for this isochrone.
            mass_ini = isoch[m_ini_idx]

            # Mass-ratio distribution
            mass_ratios = qDistribution(mass_ini, self.gamma, self.seed + mx + ax)
            # Secondary masses
            m2 = mass_ratios * mass_ini

            # Calculate unresolved binary magnitude
            mag = isoch[mag_idx]
            mag_m2 = np.interp(m2, mass_ini, mag)
            theor_tracks[mx][ax][m_ini_idx + 1] = mag_combine(mag, mag_m2)

            # Calculate unresolved color for each color defined.
            for ic, color in enumerate(all_colors):
                mc1 = self.isochs.color_filters[mx][ax][color[0]]
                mc2 = self.isochs.color_filters[mx][ax][color[1]]
                f_m1 = np.interp(m2, mass_ini, mc1)
                f1 = mag_combine(mc1, f_m1)
                f_m2 = np.interp(m2, mass_ini, mc2)
                f2 = mag_combine(mc2, f_m2)
                theor_tracks[mx][ax][m_ini_idx + 1 + 1 + ic] = f1 - f2

            # Secondary masses
            theor_tracks[mx][ax][-1] = m2

    return theor_tracks


def ccmo_ext_coeffs(self) -> list:
    """Obtain extinction coefficients for all the observed filters and colors,
    in the order in which they are stored in theor_tracks.

    ext_coefs = [ec_mag, ec_col1, ...]

    """
    ext_coefs = [[], []]
    if self.ext_law != "GAIADR3":
        # For the magnitude.
        # Effective wavelength in Armstrong.
        eff_wave = self.isochs.magnitude_effl
        # Effective wavelength in inverse microns.
        ext_coefs[0] = ccmo_model(10000.0 / eff_wave)

        # For colors.
        eff_wave1, eff_wave2 = self.isochs.color_effl
        ext_coefs[1] = [
            ccmo_model(10000.0 / eff_wave1),
            ccmo_model(10000.0 / eff_wave2),
        ]

    if self.isochs.color2_effl is not None:
        eff_wave1, eff_wave2 = self.isochs.color2_effl
        ext_coefs += [
            [ccmo_model(10000.0 / eff_wave1), ccmo_model(10000.0 / eff_wave2)]
        ]

    return ext_coefs


def ccmo_model(mw: float) -> list:
    """Cardelli, Clayton, and Mathis (1989 ApJ. 345, 245) model for extinction
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


def sample_imf(self, Nmets: int, Nages: int) -> list:
    """Returns the number of stars per interval of mass for the selected IMF.

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
      One list per metallicity and age values defined is returned. This is to add some
      variety to the sampled IMF.

    """
    inv_cdf = invTrnsfSmpl(self.IMF_name)

    st_dist_mass, st_dist_mass_ordered = [], []
    for i in range(Nmets):
        met_lst, met_lst_ord = [], []
        for j in range(Nages):
            sampled_IMF = sampleInv(i + j + self.seed, self.max_mass, inv_cdf)
            met_lst.append(sampled_IMF)
            met_lst_ord.append(np.sort(sampled_IMF))
        st_dist_mass.append(met_lst)
        st_dist_mass_ordered.append(met_lst_ord)

    return st_dist_mass, st_dist_mass_ordered


def randVals(self) -> dict:
    """Generate lists of random values used by the synthetic cluster generating
    function.

    """
    # This is the maximum number of stars that will ever be interpolated into
    # an isochrone
    N_isoch, N_mass = self.theor_tracks.shape[-1], 0
    for sdm in self.st_dist_mass:
        N_mass = max(len(sdm[0]), N_mass, N_isoch)

    # Used by `move_isochrone()` and `add_errors`
    # rand_norm_vals = np.random.normal(0.0, 1.0, (2, N_mass))
    rand_norm_vals = np.random.default_rng(self.seed).normal(0.0, 1.0, (2, N_mass))

    # Used by `move_isochrone()`, `binarity()`
    # rand_unif_vals = np.random.uniform(0.0, 1.0, (2, N_mass))
    rand_unif_vals = np.random.default_rng(self.seed).uniform(0.0, 1.0, (2, N_mass))

    rand_floats = {"norm": rand_norm_vals, "unif": rand_unif_vals}
    return rand_floats


def qDistribution(M1: np.ndarray, gamma: [float, str], seed: int) -> np.ndarray:
    """Distribution of q=m2/m1 for binary systems
    float : Power-law distribution with shape parameter 'gamma'. Not mass
    dependent.
    D&K : Distribution of mass-ratios versus primary masses
    (Duchene & Kraus 2013). Mass dependent.

    Use 'gamma + 1' in the power-law distribution below because in D&K this
    distribution is defined as f(q)~q^gamma, while numpy's distribution is
    defined as a*x^(a-1).

    """
    N = M1.size

    try:
        gamma = float(gamma)
        # mass_ratios = np.random.power(gamma + 1, N)
        mass_ratios = np.random.default_rng(seed).power(gamma + 1, N)
    except ValueError:
        if gamma == "D&K":
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
                # q = np.random.power(gammaX + 1, msk.sum())
                q = np.random.default_rng(seed).power(gammaX + 1, msk.sum())
                mass_ratios[msk] = q

        # Fisher's distribution
        elif gamma == "fisher_stepped":
            # Fisher, Schröder & Smith (2005); Table 3, stepped
            # https://doi.org/10.1111/j.1365-2966.2005.09193.x
            xk = np.linspace(0.0, 1.0, 10)
            pk = np.array([29.0, 29.0, 30.0, 32.0, 31.0, 32.0, 36.0, 45.0, 27.0, 76.0])

        elif gamma == "fisher_peaked":
            # Fisher, Schröder & Smith (2005); Table 3, peaked
            # https://doi.org/10.1111/j.1365-2966.2005.09193.x
            xk = np.linspace(0.0, 1.0, 10)
            pk = np.array([27.0, 30.0, 34.0, 33.0, 29.0, 26.0, 27.0, 33.0, 41.0, 89.0])

        elif gamma == "raghavan":
            # Raghavan et al. (2010); Fig 16 (left)
            # https://iopscience.iop.org/article/10.1088/0067-0049/190/1/1
            xk = np.linspace(0.0, 1.0, 20)
            pk = np.array(
                [
                    0.53,
                    2.61,
                    0.53,
                    4.67,
                    7.81,
                    3.64,
                    9.89,
                    5.71,
                    4.69,
                    5.73,
                    4.67,
                    6.76,
                    5.75,
                    5.73,
                    2.61,
                    5.71,
                    4.72,
                    5.71,
                    3.64,
                    12.99,
                ]
            )

        if gamma != "D&K":

            def fQ(xk, pk):
                """
                Discrete function
                """
                pk /= pk.sum()
                fq = stats.rv_discrete(a=0.0, b=1.0, values=(xk, pk))
                return fq

            fq = fQ(xk, pk)
            # 'ppf' is the inverse CDF
            mass_ratios = fq.ppf(
                # np.random.uniform(0.0, 1.0, N)
                np.random.default_rng(seed).uniform(0.0, 1.0, N)
            )

    return mass_ratios


def mag_combine(m1: np.ndarray, m2: np.ndarray) -> np.ndarray:
    """Combine two magnitudes.
    This is a faster re-ordering of the standard formula:

    -2.5 * np.log10(10 ** (-0.4 * m1) + 10 ** (-0.4 * m2))

    """
    c = 10**-0.4
    mbin = -2.5 * (-0.4 * m1 + np.log10(1.0 + c ** (m2 - m1)))
    return mbin


def properModel(
    met_age_dict: dict, fix_params: dict, fit_params: dict
) -> tuple[float, float, float, float, float, float, float, int, int, int, int]:
    """Define the 'proper' model with values for (z, a) taken from its grid,
    and filled values for those parameters that are fixed.

    Parameters
    ----------
    model : array
      Array of *free* fundamental parameters only (ie: in varIdxs).

    Returns
    -------
    fix_params  : array
      All fundamental parameters, including the fixed parameters that are
      missing from 'model'.
    ml, mh, al, ah : ints
      Indexes of the (z, a) values in the grid that define the box that enclose
      the proper (z, a) values.

    """
    # If any parameter is repeated in both dictionaries, this order gives priority
    # to the 'fix_params' dict, in the sense that its values will be the ones
    # written to 'model_comb'
    model_comb = fit_params | fix_params
    model_full = np.array(
        [
            model_comb[k]
            for k in ["met", "loga", "alpha", "beta", "Av", "DR", "Rv", "dm"]
        ]
    )

    if len(met_age_dict["met"]) == 1:
        ml = mh = 0
    else:
        par = met_age_dict["met"]
        mh = min(len(par) - 1, np.searchsorted(par, model_full[0]))
        ml = mh - 1

    if len(met_age_dict["loga"]) == 1:
        al = ah = 0
    else:
        par = met_age_dict["loga"]
        ah = min(len(par) - 1, np.searchsorted(par, model_full[1]))
        al = ah - 1

    return *model_full, ml, mh, al, ah


def zaWAverage(theor_tracks, met_age_dict, m_ini_idx, z_model, a_model, ml, mh, al, ah):
    """Generate a new "weighted" isochrone from the four closest points in the
    (z, a) grid.

    theor_tracks = [m1, m2, .., mN]
    mX = [age1, age2, ..., ageM]
    ageX = [f1,.., c1, c2,.., M_ini, f1b,.., c1b, c2b,.., M_b]
    where:
    fX:  individual filters (mags)
    cX:  colors
    M_ini: initial mass
    fXb: filters with binary data
    cXb: colors with the binary data
    M_b:  binary masses

    It is important that the returned isochrone is a *copy* of the
    theor_tracks[mx][ax] array if no averaging is done. Otherwise the
    'theor_tracks' are modified.

    """

    # If (z, a) are both fixed, just return the single processed isochrone
    if ml == al == mh == ah == 0:
        # The np.array() is important to avoid overwriting 'theor_tracks'
        return np.array(theor_tracks[ml][al])

    # The four points in the (z, age) grid that define the box that contains
    # the model value (z_model, a_model)
    z1, z2 = met_age_dict["met"][ml], met_age_dict["met"][mh]
    a1, a2 = met_age_dict["loga"][al], met_age_dict["loga"][ah]
    pts = np.array([(z1, a1), (z1, a2), (z2, a1), (z2, a2)])

    # Order: (z1, a1), (z1, a2), (z2, a1), (z2, a2)
    isochs = np.array(
        [
            theor_tracks[ml][al],
            theor_tracks[ml][ah],
            theor_tracks[mh][al],
            theor_tracks[mh][ah],
        ]
    )

    # Distances between the (z, a) points in the 'model', and the four
    # points in the (z, a) grid that contain the model point.
    # Fast euclidean distance: https://stackoverflow.com/a/47775357/1391441
    a_min_b = np.array([(z_model, a_model)]) - pts
    # Don't take the square root, it's not necessary
    dist = np.einsum("ij,ij->i", a_min_b, a_min_b)

    # If the model has a 0. distance in (z, a) to the closest isochrone,
    # then just return that isochrone.
    try:
        idx = np.where(dist == 0.0)[0][0]
        return isochs[idx]
    except IndexError:
        pass

    # Weighted average by the (inverse) distance to the four (z, a) grid
    # points. This way is faster than using 'np.average()'.
    # Inverse of the distance.
    inv_d = 1.0 / dist
    weights = inv_d / sum(inv_d)
    isochrone = (
        isochs[0] * weights[0]
        + isochs[1] * weights[1]
        + isochs[2] * weights[2]
        + isochs[3] * weights[3]
    )

    # DO NOT average the masses or their distribution will be lost. We use the
    # values of the closest isochrone.
    idx = np.argmin(dist)
    isochrone[m_ini_idx] = isochs[idx][m_ini_idx]
    # Now for the secondary masses
    isochrone[-1] = isochs[idx][-1]

    # met_idx = np.array([ml, mh, ml, mh])[idx]
    # age_idx = np.array([al, ah, al, ah])[idx]

    # isochrone = theor_tracks[ml][al] * weights[0] +\
    #     theor_tracks[ml][ah] * weights[1] +\
    #     theor_tracks[mh][al] * weights[2] +\
    #     theor_tracks[mh][ah] * weights[3]

    # #
    # isochrone = np.average(isochs, weights=weights, axis=0)

    # idxs = (weights * mass_idx).astype(int)
    # # Order: (z1, a1), (z1, a2), (z2, a1), (z2, a2)
    # isochrone = np.concatenate([
    #     theor_tracks[ml][al][:, :idxs[0]],
    #     theor_tracks[ml][ah][:, idxs[0]:idxs[0] + idxs[1]],
    #     theor_tracks[mh][al][:, idxs[0] + idxs[1]:idxs[0] + idxs[1] + idxs[2]],
    #     theor_tracks[mh][ah][:, idxs[0] + idxs[1] + idxs[2]:mass_idx]],
    #     axis=1)

    # import matplotlib.pyplot as plt
    # print(z_model, a_model, ml, mh, al, ah)
    # plt.subplot(121)
    # plt.scatter(*pts.T, c='r')
    # plt.scatter(z_model, a_model, marker='x', c='g')
    # # First color
    # plt.subplot(122)
    # plt.scatter(theor_tracks[ml][al][1], theor_tracks[ml][al][0], c='b', alpha=.25)
    # plt.scatter(theor_tracks[ml][ah][1], theor_tracks[ml][ah][0], c='r', alpha=.25)
    # plt.scatter(theor_tracks[mh][al][1], theor_tracks[mh][al][0], c='cyan', alpha=.25)
    # plt.scatter(theor_tracks[mh][ah][1], theor_tracks[mh][ah][0], c='orange', alpha=.25)
    # plt.scatter(isochrone[1], isochrone[0], c='k', marker='x', alpha=.5, label='avrg')
    # plt.gca().invert_yaxis()
    # plt.legend()
    # # Second color
    # # plt.subplot(133)
    # # plt.scatter(theor_tracks[ml][al][2], theor_tracks[ml][al][0], c='b')
    # # plt.scatter(theor_tracks[ml][ah][2], theor_tracks[ml][ah][0], c='r')
    # # plt.scatter(theor_tracks[mh][al][2], theor_tracks[mh][al][0], c='cyan')
    # # plt.scatter(theor_tracks[mh][ah][2], theor_tracks[mh][ah][0], c='orange')
    # # plt.scatter(isochrone[2], isochrone[0], c='g', ls='--')
    # # plt.gca().invert_yaxis()
    # plt.show()
    # # breakpoint()

    return isochrone


def move_isochrone(isochrone, m_ini_idx, dm):
    """
    Receives an isochrone of a given age and metallicity and modifies
    its magnitude values according to a given distance modulus.
    """
    # Move magnitude
    isochrone[0] += dm
    # Move binary magnitude if it exists
    if isochrone.shape[0] > m_ini_idx + 1:
        isochrone[m_ini_idx + 1] += dm

    return isochrone


def extinction(
    ext_law,
    ext_coefs,
    rand_norm,
    rand_unif,
    DR_distribution,
    m_ini_idx,
    binar_flag,
    Av,
    dr,
    Rv,
    isochrone,
):
    """
    Modifies magnitude and color(s) according to given values for the
    total absorption Av. Using this parameter instead of the E(B-V) extinction
    reduces the correlation with Rv.

    The distance modulus was applied before this function.

    isochrone = [mag, c1, (c2), .., Mini, mag_b, c1b, .., Mini_b]

    For the CCMO model:

    ext_coefs = [mag_ec, c1_ec, ...]

    where:
    mag_ec = [a, b]  ; cX_ec = [[a1, b1], [a2, b2]]

    and:
    ccm_coef = a + b / Rv = ext_coefs[i][0] + ext_coefs[i][1] / Rv

    Ax = ef * Av
    m_obs = M_int + Ax + dist_mod
          = M_int + ef * R_V * E(B-V) + dist_mod
          = M_int + (a + b / Rv) * R_V * E(B-V) + dist_mod

    E(m1 - m2) = A_m1 - A_m2
               = (ef_m1 - ef_m2) * Av
               = [(a1 + b1/Rv) - (a2 + b2/Rv)] * Av
               = [(a1 - a2) + (b1 - b2)/Rv] * Av
               = (a12 + b12/Rv) * Av
               = (a12 + b12/Rv) * R_V * E(B-V)
    (m1 - m2)_obs = (m1 - m2)_int + E(m1 - m2)
    (m1 - m2)_obs = (m1 - m2)_int + (a12 + b12/Rv) * R_V * E(B-V)

    """

    if dr > 0.0:
        Ns = isochrone.shape[-1]

        if DR_distribution == "uniform":
            # Av_dr = rand_unif[:Ns] * dr
            Av_dr = (2 * rand_unif[:Ns] - 1) * dr
        elif DR_distribution == "normal":
            # Av_dr = abs(rand_norm[:Ns]) * dr
            Av_dr = rand_norm[:Ns] * dr

        # In place in case I ever want to implement the percentage of stars affected.
        # Without this, all stars are affected by the DR.
        # Av_dr[rand_unif[:Ns] > DR_percentage] = 0.0

        # Clip at 0
        Av = np.clip(Av + Av_dr, a_min=0, a_max=np.inf)

    if ext_law == "CCMO":
        # Magnitude
        ec_mag = ext_coefs[0][0] + ext_coefs[0][1] / Rv
        # First color
        ec_col1 = (ext_coefs[1][0][0] + ext_coefs[1][0][1] / Rv) - (
            ext_coefs[1][1][0] + ext_coefs[1][1][1] / Rv
        )
    elif ext_law == "GAIADR3":
        # If this model is used the first color is always expected to be BP-RP
        # BP_RP = isochrone[1]
        ec_mag, ec_col1 = dustapprox(isochrone[1], Av)

    Ax = ec_mag * Av
    isochrone[0] += Ax
    Ex1 = ec_col1 * Av
    isochrone[1] += Ex1

    # Move binary data.
    if binar_flag:
        isochrone[m_ini_idx + 1] += Ax  # Magnitude
        isochrone[m_ini_idx + 2] += Ex1  # First color

    # Second color
    if len(ext_coefs) > 2:
        ec_col2 = (ext_coefs[2][0][0] + ext_coefs[2][0][1] / Rv) - (
            ext_coefs[2][1][0] + ext_coefs[2][1][1] / Rv
        )
        Ex2 = ec_col2 * Av
        isochrone[2] += Ex2
        # Move color with binary data.
        if binar_flag:
            isochrone[m_ini_idx + 3] += Ex2

    return isochrone


def dustapprox(X_, Av):
    """
    The 'coeffs' values are the main sequence values taken from:

    https://www.cosmos.esa.int/web/gaia/edr3-extinction-law

    The order of the coefficients is:

    Intercept   X   X2  X3  A   A2  A3  XA  AX2 XA2

    X_ == BP-RP
    """
    coeffs = {
        "G": (
            0.995969721536602,
            -0.159726460302015,
            0.0122380738156057,
            0.00090726555099859,
            -0.0377160263914123,
            0.00151347495244888,
            -2.52364537395142e-05,
            0.0114522658102451,
            -0.000936914989014318,
            -0.000260296774134201,
        ),
        "BP": (
            1.15363197483424,
            -0.0814012991657388,
            -0.036013023976704,
            0.0192143585568966,
            -0.022397548243016,
            0.000840562680547171,
            -1.31018008013549e-05,
            0.00660124080271006,
            -0.000882247501989453,
            -0.000111215755291684,
        ),
        "RP": (
            0.66320787941067,
            -0.0179847164933981,
            0.000493769449961458,
            -0.00267994405695751,
            -0.00651422146709376,
            3.30179903473159e-05,
            1.57894227641527e-06,
            -7.9800898337247e-05,
            0.000255679812110045,
            1.10476584967393e-05,
        ),
    }

    X_2 = X_**2
    X_3 = X_**3
    Av_2 = Av**2
    Av_3 = Av**3

    def ext_coeff(k):
        """
        https://www.cosmos.esa.int/web/gaia/edr3-extinction-law
        """
        # X   X2  X3  A   A2  A3  XA  AX2 XA2
        ay = coeffs[k][0]
        for i, Xk in enumerate([X_, X_2, X_3]):
            ay += coeffs[k][1 + i] * Xk
        for i, Ak in enumerate([Av, Av_2, Av_3]):
            ay += coeffs[k][4 + i] * Ak

        ay += (
            coeffs[k][7] * X_ * Av
            + coeffs[k][9] * X_ * Av_2  # This index not a mistake
            + coeffs[k][8] * X_2 * Av
        )
        return ay

    ec_G = ext_coeff("G")
    ec_BPRP = ext_coeff("BP") - ext_coeff("RP")

    return ec_G, ec_BPRP


def cut_max_mag(isoch_moved, max_mag_syn):
    """
    Remove stars from isochrone with magnitude values larger that the maximum
    observed value.
    """
    # Discard stars in isochrone beyond max_mag_syn limit.
    msk = isoch_moved[0] < max_mag_syn
    return isoch_moved[:, msk]


def mass_interp(isoch_cut, m_ini_idx, st_dist_mass, N_obs_stars):
    """
    For each mass in the sampled IMF mass distribution, interpolate its value
    (and those of all the sub-arrays in 'isoch_cut') into the isochrone.

    Masses that fall outside of the isochrone's mass range have been previously
    rejected.
    """
    # Assumes `mass_ini=isoch_cut[m_ini_idx]` is sorted min to max <-- IMPORTANT
    mass_ini = isoch_cut[m_ini_idx]

    # Filter masses in the IMF sampling that are outside of the mass
    # range given by 'isoch_cut' (st_dist_mass[0]: sampled masses from IMF)
    # msk_min = (st_dist_mass[0] >= mass_ini.min())
    # msk_max = (st_dist_mass[0] <= mass_ini.max())
    # (~msk_min).sum(): stars lost below the minimum mass (photometric)
    # (~msk_max).sum(): stars lost above the maximum mass (evolutionary)
    msk_m = (st_dist_mass >= mass_ini.min()) & (st_dist_mass <= mass_ini.max())

    # Interpolate the same number of observed stars into the isochrone
    mass_dist = st_dist_mass[msk_m][:N_obs_stars]
    if not mass_dist.any():
        return np.array([])

    # # This is about 2x faster than interp1d() but it comes at the cost of a more
    # # coarse distribution of stars throughout the isochrone. To fix this, we have to
    # # apply a random noise to the magnitude(s) proportional to the percentage that
    # # the masses sampled differ from those in the isochrone. This lowers the
    # # performance to an increase of ~22%
    # idx = np.searchsorted(mass_ini, mass_dist)
    # isoch_mass = isoch_cut[:, idx]
    # m_perc = (isoch_mass[m_ini_idx]-mass_dist)/mass_dist
    # isoch_mass[0] += isoch_mass[0]*m_perc
    # if isoch_mass.shape[0] > m_ini_idx:
    #     isoch_mass[m_ini_idx+1] += isoch_mass[m_ini_idx+1]*m_perc
    # return isoch_mass

    # # This is equivalent to the the block below but slower for more than 5
    # # dimensions, and *very* slightly faster for lower dimensions
    # isoch_mass = np.empty([isoch_cut.shape[0], mass_dist.size])
    # for i, arr in enumerate(isoch_cut):
    #     isoch_mass[i] = np.interp(mass_dist, mass_ini, arr)

    # Interpolate the sampled stars (masses) into the isochrone
    isoch_mass = interp_mass_isoch(isoch_cut, mass_ini, mass_dist)

    return isoch_mass


def interp_mass_isoch(isoch_cut, mass_ini, mass_dist):
    """
    Find where in the original data, the values to interpolate would be inserted.

    NOTE: I already tried to speed this block up using numba (@jit(nopython=True))
    but it does not help. The code runs even slower.
    """
    # Note: If x_new[n] == x[m], then m is returned by searchsorted.
    x_new_indices = np.searchsorted(mass_ini, mass_dist)

    # Calculate the slope of regions that each x_new value falls in.
    lo = x_new_indices - 1
    x_lo = mass_ini[lo]
    x_hi = mass_ini[x_new_indices]
    y_lo = isoch_cut[:, lo]
    y_hi = isoch_cut[:, x_new_indices]
    slope = (y_hi - y_lo) / (x_hi - x_lo)
    #
    x_diff = mass_dist - x_lo
    y_diff = slope * x_diff
    # Calculate the actual value for each entry in x_new.
    isoch_mass = y_diff + y_lo

    return isoch_mass


def binarity(alpha, beta, binar_flag, m_ini_idx, rand_unif_vals, isoch_mass):
    """
    Select a fraction of stars to be binaries, given a chosen method.
    """
    # No binarity process defined
    if binar_flag is False:
        # Update the binary systems' masses so that the secondary masses for
        # SINGLE systems are identified with a 'nan' value.
        isoch_mass[-1] = np.nan
        return isoch_mass

    Ns = isoch_mass.shape[-1]
    mass = isoch_mass[m_ini_idx]

    # Distribution of probability of binarity (multiplicity fraction) versus
    # primary masses.

    # Offner et al. (2022); Fig 1 (left), Table 1
    # b_p = np.clip(alpha + beta * np.log(mass), a_min=0, a_max=1)
    # b_p = np.clip(alpha + beta * np.arctan(mass), a_min=0, a_max=1)
    b_p = np.clip(alpha + beta * (1 / (1 + 1.4 / mass)), a_min=0, a_max=1)

    # Stars (masses) with the largest binary probabilities are selected
    # proportional to their probability
    bin_indxs = b_p > rand_unif_vals[:Ns]

    # Index of the binary magnitude: mag_binar = m_ini_idx + 1
    # Update array with new values of magnitudes and colors for the binary
    # systems. This does not change the primary mass values, just the magnitude and
    # color(s).
    if bin_indxs.any():
        for i in range(m_ini_idx):
            isoch_mass[i][bin_indxs] = isoch_mass[m_ini_idx + 1 + i][bin_indxs]

    # Update the binary systems' masses so that the secondary masses for
    # SINGLE systems are identified with a 'nan' value.
    isoch_mass[-1][~bin_indxs] = np.nan

    return isoch_mass


def add_errors(isoch_binar, err_dist):
    """
    Add random synthetic uncertainties to the magnitude and color(s)
    """

    N = len(isoch_binar[0])
    mag_sort = np.argsort(-isoch_binar[0])
    for i, sigma in enumerate(err_dist):
        isoch_binar[i][mag_sort] += sigma[:N]

    return isoch_binar


# def _rm_low_masses(self, dm_min):
#     """
#     dm_min: float | None = None

#     dm_min : float, optional, default=None
#         Value for the minimum distance modulus. Used to constrain the lower masses
#         in the theoretical isochrones to make the generating process more
#         efficient.
#     """
#     min_masses = []
#     for met_arr in self.theor_tracks:
#         met_lst = []
#         for age_arr in met_arr:
#             mag, mass = age_arr[0], age_arr[self.m_ini_idx]
#             i = np.argmin(abs(self.max_mag_syn - (mag + dm_min)))
#             met_lst.append(mass[i])
#         min_masses.append(met_lst)

#     st_dist_mass_lmass = []
#     for i, met_arr in enumerate(self.st_dist_mass):
#         met_lst = []
#         for j, mass_sample in enumerate(met_arr):
#             min_mass = min_masses[i][j]
#             msk = mass_sample > min_mass
#             sampled_IMF = mass_sample[msk]
#             met_lst.append(sampled_IMF)
#         st_dist_mass_lmass.append(met_lst)

#     # Make copy of original array, used for mass estimation in cluster() class
#     self.st_dist_mass_full = self.st_dist_mass.copy()
#     # Update this parameter with the new array
#     self.st_dist_mass = st_dist_mass_lmass
