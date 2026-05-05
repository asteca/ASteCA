import numpy as np
import numpy.typing as npt
from scipy import stats

from .imfs import invTrnsfSmpl, sampleInv


def sample_imf(
    rng: np.random.Generator, IMF_name: str, max_mass: float, Nmets: int, Nages: int
) -> tuple[list, list]:
    """Returns arrays of sampled stars for the selected IMF.

    :param rng: Random number generator.
    :type rng: np.random.Generator
    :param IMF_name: Name of the IMF to be used.
    :type IMF_name: str
    :param max_mass: Maximum mass defined.
    :type max_mass: float
    :param Nmets: Number of metallicity values.
    :type Nmets: int
    :param Nages: Number of age values.
    :type Nages: int

    :returns: A tuple containing two lists. The first list contains the sampled masses,
     and the second list contains the ordered sampled masses.
    :rtype: tuple[list, list]
    """
    inv_cdf = invTrnsfSmpl(IMF_name)

    # Sample in chunks until the maximum defined mass is reached. A simple
    # analysis points to this being the optimal N_chunk
    N_chunk = max(100, int(max_mass / 40))

    st_dist_mass, st_dist_mass_ordered = [], []
    for i in range(Nmets):
        met_lst, met_lst_ord = [], []
        for j in range(Nages):
            sampled_IMF = sampleInv(rng, max_mass, inv_cdf, N_chunk)
            met_lst.append(sampled_IMF)
            met_lst_ord.append(np.sort(sampled_IMF))
        st_dist_mass.append(met_lst)
        st_dist_mass_ordered.append(met_lst_ord)

    return st_dist_mass, st_dist_mass_ordered


def add_binarity(
    rng: np.random.Generator,
    gamma: float | str,
    color: tuple,
    color2: tuple | None,
    theor_tracks: np.ndarray,
    color_filters: list,
) -> np.ndarray:
    """For each theoretical isochrone defined.
        1. Draw random secondary masses for *all* stars.
        2. Interpolate magnitude values for these masses
        3. Calculate the magnitudes for the binary system.
        4. Calculate the colors for the binary system.

    theor_tracks.shape = (Nz, Na, N_cols, N_interp)

    If binarity is processed:
        N_cols: magnitude + color(s) + initial masses + secondary masses +
                unresolved mag + unresolved color(s)
    else:
        N_cols: magnitude + colors + initial mass

    :param rng: Random number generator.
    :type rng: np.random.Generator
    :param gamma: Mass-ratio distribution.
    :type gamma: float | str
    :param color: First color defined
    :type color: tuple
    :param color2: Second color defined (optional)
    :type color2: tuple | None
    :param theor_tracks: Array with the processed theoretical isochrones
    :type theor_tracks: np.ndarray
    :param color_filters: Magnitude values required to generate the colors separately
    :type color_filters: list

    :returns: Array of isochrones with binary data.
    :rtype: np.ndarray
    """
    all_colors = [color]
    if color2 is not None:
        all_colors.append(color2)
    N_colors = len(all_colors)

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
            mass_ratios = qDistribution(mass_ini, gamma, rng)
            # Secondary masses
            m2 = mass_ratios * mass_ini

            # Calculate unresolved binary magnitude
            mag = isoch[mag_idx]
            mag_m2 = np.interp(m2, mass_ini, mag)
            theor_tracks[mx][ax][m_ini_idx + 2] = mag_combine(mag, mag_m2)

            # Calculate unresolved color for each color defined.
            for ic, color in enumerate(all_colors):
                mc1 = color_filters[mx][ax][color[0]]
                mc2 = color_filters[mx][ax][color[1]]
                f_m1 = np.interp(m2, mass_ini, mc1)
                f1 = mag_combine(mc1, f_m1)
                f_m2 = np.interp(m2, mass_ini, mc2)
                f2 = mag_combine(mc2, f_m2)
                theor_tracks[mx][ax][m_ini_idx + 2 + 1 + ic] = f1 - f2

            # Secondary masses
            theor_tracks[mx][ax][m_ini_idx + 1] = m2

    return theor_tracks


def qDistribution(
    M1: np.ndarray, gamma: float | str, rng: np.random.Generator
) -> np.ndarray:
    """Distribution of q=m2/m1 for binary systems

    :param M1: Array of primary masses.
    :type M1: np.ndarray
    :param gamma: Mass-ratio distribution. Can either be a float value (in which case
     a non mass-dependent power-law distribution with shape parameter 'gamma' will be
     used), 'D&K' (mass dependent distribution of mass-ratios versus primary masses
      from Duchene & Kraus 2013), or one of Fisher or Raghavan distributions (non
      mass-dependent).
    :type gamma: float | str
    :param rng: Random number generator.
    :type rng: np.random.Generator

    :returns: Array of mass ratios.
    :rtype: np.ndarray
    """
    N = M1.size

    try:
        gamma = float(gamma)
        # mass_ratios = np.random.power(gamma + 1, N)  # DEPRECATED 13/01/25

        # Use 'gamma + 1' in the power-law distribution because in D&K this
        # distribution is defined as f(q)~q^gamma, while numpy's distribution is
        # defined as a*x^(a-1).
        mass_ratios = rng.power(gamma + 1, N)

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
                q = rng.power(gammaX + 1, msk.sum())
                mass_ratios[msk] = q
        else:

            def fQ(xk: np.ndarray, pk: np.ndarray):
                """
                Discrete function
                """
                pk /= pk.sum()
                # pyright error: https://github.com/scipy/scipy/issues/22327
                fq = stats.rv_discrete(a=0.0, b=1.0, values=(xk, pk))
                return fq

            # Fisher's distribution
            # "fisher_stepped":
            # Fisher, Schröder & Smith (2005); Table 3, stepped
            # https://doi.org/10.1111/j.1365-2966.2005.09193.x
            #
            # "fisher_peaked":
            # Fisher, Schröder & Smith (2005); Table 3, peaked
            # https://doi.org/10.1111/j.1365-2966.2005.09193.x
            #
            # "raghavan":
            # Raghavan et al. (2010); Fig 16 (left)
            # https://iopscience.iop.org/article/10.1088/0067-0049/190/1/1

            xy_dict = {
                "fisher_stepped": [
                    np.linspace(0.0, 1.0, 10),
                    np.array(
                        [29.0, 29.0, 30.0, 32.0, 31.0, 32.0, 36.0, 45.0, 27.0, 76.0]
                    ),
                ],
                "fisher_peaked": [
                    np.linspace(0.0, 1.0, 10),
                    np.array(
                        [27.0, 30.0, 34.0, 33.0, 29.0, 26.0, 27.0, 33.0, 41.0, 89.0]
                    ),
                ],
                "raghavan": [
                    np.linspace(0.0, 1.0, 20),
                    np.array(
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
                    ),
                ],
            }

            xk, pk = xy_dict[str(gamma)]

            # # Fisher's distribution
            # if gamma == "fisher_stepped":
            #     # Fisher, Schröder & Smith (2005); Table 3, stepped
            #     # https://doi.org/10.1111/j.1365-2966.2005.09193.x
            #     xk = np.linspace(0.0, 1.0, 10)
            #     pk = np.array([29.0, 29.0, 30.0, 32.0, 31.0, 32.0, 36.0, 45.0, 27.0, 76.0])

            # elif gamma == "fisher_peaked":
            #     # Fisher, Schröder & Smith (2005); Table 3, peaked
            #     # https://doi.org/10.1111/j.1365-2966.2005.09193.x
            #     xk = np.linspace(0.0, 1.0, 10)
            #     pk = np.array([27.0, 30.0, 34.0, 33.0, 29.0, 26.0, 27.0, 33.0, 41.0, 89.0])

            # elif gamma == "raghavan":
            #     # Raghavan et al. (2010); Fig 16 (left)
            #     # https://iopscience.iop.org/article/10.1088/0067-0049/190/1/1
            #     xk = np.linspace(0.0, 1.0, 20)
            #     pk = np.array(
            #         [
            #             0.53,
            #             2.61,
            #             0.53,
            #             4.67,
            #             7.81,
            #             3.64,
            #             9.89,
            #             5.71,
            #             4.69,
            #             5.73,
            #             4.67,
            #             6.76,
            #             5.75,
            #             5.73,
            #             2.61,
            #             5.71,
            #             4.72,
            #             5.71,
            #             3.64,
            #             12.99,
            #         ]
            #     )

            fq = fQ(xk, pk)
            # 'ppf' is the inverse CDF
            mass_ratios = fq.ppf(
                # np.random.uniform(0.0, 1.0, N)
                rng.uniform(0.0, 1.0, N)
            )

    return mass_ratios


def mag_combine(
    m1: npt.NDArray[np.floating], m2: npt.NDArray[np.floating]
) -> np.ndarray:
    """Combine two magnitudes.

    This is a faster re-ordering of the standard formula:
    -2.5 * np.log10(10 ** (-0.4 * m1) + 10 ** (-0.4 * m2))

    :param m1: Array of magnitudes.
    :type m1: npt.NDArray[np.floating]
    :param m2: Array of magnitudes.
    :type m2: npt.NDArray[np.floating]

    :returns: Array of combined magnitudes.
    :rtype: np.ndarray
    """
    c = 10**-0.4
    mbin = -2.5 * (-0.4 * m1 + np.log10(1.0 + c ** (m2 - m1)))
    return mbin


def ccmo_ext_coeffs(
    magnitude_effl: float,
    color_effl: tuple,
    color2_effl: tuple | None,
) -> list:
    """Obtain extinction coefficients for all the observed filters and colors,
    in the order in which they are stored in theor_tracks.

    ext_coefs = [ec_mag, ec_col1, ...]

    :param magnitude_effl: Effective lambda (in Angstrom) for the magnitude filter
    :type magnitude_effl: float
    :param color_effl: Effective lambdas for the filters that make up the first color
    :type color_effl: tuple
    :param color2_effl: Effective lambdas for the filters that make up the second color
    :type color2_effl: tuple | None

    :returns: List of extinction coefficients.
    :rtype: list

    """
    # Effective wavelength in Armstrong.
    eff_wave = magnitude_effl
    eff_wave1, eff_wave2 = color_effl
    # Effective wavelength in inverse microns.
    ext_coefs = [
        ccmo_model(10000.0 / eff_wave),
        [
            ccmo_model(10000.0 / eff_wave1),
            ccmo_model(10000.0 / eff_wave2),
        ],
    ]

    if color2_effl is not None:
        eff_wave1, eff_wave2 = color2_effl
        ext_coefs += [
            [ccmo_model(10000.0 / eff_wave1), ccmo_model(10000.0 / eff_wave2)]
        ]

    return ext_coefs


def ccmo_model(mw: float) -> tuple[float, float]:
    """Cardelli, Clayton, and Mathis (1989 ApJ. 345, 245) model for extinction
    coefficients with updated coefficients for near-UV from O'Donnell
    (1994, ApJ, 422, 158).

    ccm_coef = a + b / Rv

    Implementation taken from:

    http://idlastro.gsfc.nasa.gov/ftp/pro/astro/ccm_unred.pro

    There appears to be an error in the Far-UV range in the original IDL
    routine where the maximum inverse wavelength is 11 and it should be 10
    according to Cardelli et al. 1989 (pag 251, Eq (5,a,b)).

    :param mw: Wavelength in inverse microns.
    :type mw: float

    :raises ValueError: If the effective wavelength is beyond the CCM model limit

    :returns: Extinction coefficients a and b.
    :rtype: tuple[float, float]
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
        c1.reverse()
        c2.reverse()
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
        c1.reverse()
        c2.reverse()
        a, b = np.polyval(c1, y), np.polyval(c2, y)
    else:
        raise ValueError(
            "The effective wavelength is {} [1/micron], beyond "
            "the CCM model limit (10 [1/micron]).".format(mw)
        )

    return float(a), float(b)


def randVals(
    rng: np.random.Generator, theor_tracks: np.ndarray, st_dist_mass: list
) -> dict:
    """Generate lists of random values used by the synthetic cluster generating
    function.

    :param rng: Random number generator.
    :type rng: np.random.Generator
    :param theor_tracks: Array of theoretical isochrones.
    :type theor_tracks: np.ndarray
    :param st_dist_mass: List of sampled masses.
    :type st_dist_mass: list

    :returns: Dictionary of random values.
    :rtype: dict
    """
    # This is the maximum number of stars that will ever be interpolated into
    # an isochrone
    N_isoch, N_mass = theor_tracks.shape[-1], 0
    for sdm in st_dist_mass:
        N_mass = max(len(sdm[0]), N_mass, N_isoch)

    # Used by `move_isochrone()` and `add_errors`
    # rand_norm_vals = np.random.normal(0.0, 1.0, (2, N_mass))
    rand_norm_vals = rng.normal(0.0, 1.0, (2, N_mass))

    # Used by `move_isochrone()`, `binarity()`
    # rand_unif_vals = np.random.uniform(0.0, 1.0, (2, N_mass))
    rand_unif_vals = rng.uniform(0.0, 1.0, (2, N_mass))

    rand_floats = {"norm": rand_norm_vals, "unif": rand_unif_vals}
    return rand_floats


def error_distribution(
    mag: np.ndarray,
    e_mag: np.ndarray,
    e_color: np.ndarray,
    e_color2: np.ndarray | None,
    rand_norm_vals: np.ndarray,
) -> list[np.ndarray]:
    """Extract the magnitude and color(s) uncertainties to use as error model for the
    synthetic clusters.

    :param mag: Array of magnitudes.
    :type mag: np.ndarray
    :param e_mag: Array of magnitude uncertainties.
    :type e_mag: np.ndarray
    :param e_color: Array of color uncertainties.
    :type e_color: np.ndarray
    :param e_color2: Array of color uncertainties.
    :type e_color2: np.ndarray | None
    :param rand_norm_vals: Array of random normal values.
    :type rand_norm_vals: np.ndarray

    :returns: List of arrays of error distributions.
    :rtype: list[np.ndarray]
    """

    def filnans(data):
        msk = np.isnan(data)
        if msk.any() and (~msk).any():
            data = np.array(data, copy=True)
            data[msk] = np.interp(np.flatnonzero(msk), np.flatnonzero(~msk), data[~msk])
        return data

    # Replace nan values with interpolated values
    mag = filnans(mag)
    e_mag = filnans(e_mag)
    e_colors = [filnans(e_color)]
    if e_color2 is not None:
        e_colors += [filnans(e_color2)]

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
