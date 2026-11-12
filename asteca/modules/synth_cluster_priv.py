import numpy as np
import numpy.typing as npt
from scipy import stats

from . import HAS_NUMBA as _HAS_NUMBA
from . import njit
from .imfs import invCDF_lut, invTrnsfSmpl, sampleInv


def sample_imf(
    rng: np.random.Generator, IMF_name: str, max_mass: float, Nmets: int, Nages: int
) -> list:
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

    :returns: Sampled masses.
    :rtype: list
    """
    inv_cdf = invTrnsfSmpl(IMF_name)

    # Sample in chunks until the maximum defined mass is reached. A simple
    # analysis points to this being the optimal N_chunk
    N_chunk = max(100, int(max_mass / 40))

    # Number of values drawn at once: the expected number of chunks required to
    # reach 'max_mass' (the mean mass is the integral of the inverse CDF in [0, 1])
    # plus a 20% margin
    mean_mass = float(inv_cdf.integrate(0.0, 1.0))
    N_draw = (int(1.2 * max_mass / (mean_mass * N_chunk)) + 2) * N_chunk

    # Knot lookup table and work array, shared by all the samples
    lut = invCDF_lut(inv_cdf.t[1:-1])
    buf = np.empty(N_draw)

    st_dist_mass = []
    for i in range(Nmets):
        met_lst = []
        for j in range(Nages):
            sampled_IMF = sampleInv(
                rng, max_mass, inv_cdf, N_chunk, N_draw, lut=lut, buf=buf
            )
            met_lst.append(sampled_IMF)
        st_dist_mass.append(met_lst)

    return st_dist_mass


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

            # Calculate unresolved binary magnitude. Secondary masses below the
            # isochrone's minimum mass are assigned an infinite magnitude (i.e.: no
            # flux), so that the binary system's photometry equals that of the primary
            mag = isoch[mag_idx]
            mag_m2 = np.interp(m2, mass_ini, mag, left=np.inf)
            theor_tracks[mx][ax][m_ini_idx + 2] = mag_combine(mag, mag_m2)

            # Calculate unresolved color for each color defined.
            for ic, _color in enumerate(all_colors):
                mc1 = color_filters[mx][ax][_color[0]]
                mc2 = color_filters[mx][ax][_color[1]]
                f_m1 = np.interp(m2, mass_ini, mc1, left=np.inf)
                f1 = mag_combine(mc1, f_m1)
                f_m2 = np.interp(m2, mass_ini, mc2, left=np.inf)
                f2 = mag_combine(mc2, f_m2)
                theor_tracks[mx][ax][m_ini_idx + 2 + 1 + ic] = f1 - f2

            # Secondary masses
            theor_tracks[mx][ax][m_ini_idx + 1] = m2

    return theor_tracks


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
    return -2.5 * (-0.4 * m1 + np.log10(1.0 + c ** (m2 - m1)))


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
            # The Fisher and Raghavan distributions are histograms of the number of
            # systems in equal width 'q' bins covering (0, 1]. Values of 'q' are
            # sampled uniformly within the bins, using the inverse of the
            # (piecewise linear) cumulative distribution.
            #
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
                "fisher_stepped": np.array(
                    [29.0, 29.0, 30.0, 32.0, 31.0, 32.0, 36.0, 45.0, 27.0, 76.0]
                ),
                "fisher_peaked": np.array(
                    [27.0, 30.0, 34.0, 33.0, 29.0, 26.0, 27.0, 33.0, 41.0, 89.0]
                ),
                "raghavan": np.array(
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
            }

            pk = xy_dict[str(gamma)]
            # Bin edges and cumulative distribution evaluated at the edges
            q_edges = np.linspace(0.0, 1.0, pk.size + 1)
            cdf = np.concatenate(([0.0], np.cumsum(pk) / pk.sum()))
            # Use (0, 1] uniform values so that q=0 is never drawn (the first bin
            # has a non-zero weight in all the distributions)
            u = 1.0 - rng.uniform(0.0, 1.0, N)
            mass_ratios = np.interp(u, cdf, q_edges)

    return mass_ratios


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
            f"The effective wavelength is {mw} [1/micron], beyond "
            "the CCM model limit (10 [1/micron])."
        )

    return float(a), float(b)


def randVals(
    rng: np.random.Generator,
    theor_tracks: np.ndarray,
    st_dist_mass: list,
    DR_distribution: str,
) -> dict:
    """Generate lists of random values used by the synthetic cluster generating
    function.

    :param rng: Random number generator.
    :type rng: np.random.Generator
    :param theor_tracks: Array of theoretical isochrones.
    :type theor_tracks: np.ndarray
    :param st_dist_mass: List of sampled masses.
    :type st_dist_mass: list
    :param DR_distribution: Differential reddening model, see
        :py:func:`dr_fractions`.
    :type DR_distribution: str

    :returns: Dictionary of random values.
    :rtype: dict
    """
    # This is the maximum number of stars that will ever be interpolated into
    # an isochrone: the largest of all the IMF samples (one per metallicity and
    # age), or the number of isochrone points. All the samples must be checked, not
    # just the first age of each metallicity: their lengths differ, and the numba
    # functions that use these arrays do not check bounds
    N_isoch = theor_tracks.shape[-1]
    N_mass = max(N_isoch, max(arr.size for met in st_dist_mass for arr in met))

    # Used by `move_isochrone()` and `add_errors`
    # rand_norm_vals = np.random.normal(0.0, 1.0, (2, N_mass))
    rand_norm_vals = rng.normal(0.0, 1.0, (2, N_mass))

    # Used by `move_isochrone()`, `binarity()`
    # rand_unif_vals = np.random.uniform(0.0, 1.0, (2, N_mass))
    rand_unif_vals = rng.uniform(0.0, 1.0, (2, N_mass))

    # Used by `extinction()`
    dr_frac = dr_fractions(rand_unif_vals[0], DR_distribution)

    rand_floats = {"norm": rand_norm_vals, "unif": rand_unif_vals, "dr_frac": dr_frac}
    return rand_floats


# Lower limits of the parameters that can not be negative (Rv must also be > 0,
# see synth_model.properModel())
PARAMS_MIN = {"alpha": 0.0, "beta": 0.0, "Av": 0.0, "DR": 0.0, "dm": 0.0, "Rv": 0.0}


def truncnorm_sample(
    rng: np.random.Generator,
    mean: float,
    std: float,
    low: float,
    high: float,
    N: int,
) -> np.ndarray:
    """Sample a normal distribution truncated to the range [low, high].

    Unlike clipping the samples of a normal distribution, this does not pile up
    values at the limits of the range.

    :param rng: Random number generator.
    :type rng: np.random.Generator
    :param mean: Mean of the (untruncated) normal distribution.
    :type mean: float
    :param std: Standard deviation of the (untruncated) normal distribution. If
        ``0``, all the values are equal to ``mean``.
    :type std: float
    :param low: Lower limit of the range.
    :type low: float
    :param high: Upper limit of the range.
    :type high: float
    :param N: Number of samples.
    :type N: int

    :returns: Array of samples.
    :rtype: np.ndarray
    """
    if std <= 0.0:
        return np.full(N, float(mean))
    a, b = (low - mean) / std, (high - mean) / std
    return stats.truncnorm.rvs(a, b, loc=mean, scale=std, size=N, random_state=rng)


def dr_fractions(rand_unif: np.ndarray, DR_distribution: str) -> np.ndarray:
    """Fractions 'f' (in [0, 1]) of the differential reddening applied to each
    point, A_i = Av + f_i * DR.

    internal: stars and dust uniformly mixed in a sphere, behind a uniform
    foreground screen (Av). DR is the absorption through the centre of the
    cluster, front to back, and 'f' the fraction of that column in front of each
    star. Its distribution is p(f) = 3/2 (1 - f^2), with inverse CDF
    f = 2 sin(arcsin(u) / 3).

    foreground: dust-free cluster behind a patchy foreground screen. Av is its
    minimum absorption and DR the full range of its variation, i.e.: 'f' follows
    a normal distribution N(0.5, 1/6) truncated at +/-3 sigma.

    :param rand_unif: Array of random uniform values in [0, 1).
    :type rand_unif: np.ndarray
    :param DR_distribution: Differential reddening model, one of
        ``internal, foreground``.
    :type DR_distribution: str

    :raises ValueError: If the differential reddening model is not recognized.

    :returns: Array of fractions.
    :rtype: np.ndarray
    """
    if DR_distribution == "internal":
        return 2.0 * np.sin(np.arcsin(rand_unif) / 3.0)
    if DR_distribution == "foreground":
        # Inverse CDF of the normal, restricted to +/-3 sigma
        p_lo, p_hi = stats.norm.cdf(-3.0), stats.norm.cdf(3.0)
        return 0.5 + stats.norm.ppf(p_lo + rand_unif * (p_hi - p_lo)) / 6.0
    raise ValueError(f"Differential reddening model '{DR_distribution}' not recognized")


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


def mass_idx_lookup(
    theor_tracks: np.ndarray, st_dist_mass: list, m_ini_idx: int
) -> dict:
    """Pre-compute the search indexes of each sampled IMF mass distribution into
    the mass grid of the isochrones, used by :py:func:`synth_model.mass_interp`.

    In :py:func:`synth_model.get_model` the masses of the isochrone are always taken
    from a single node of the (z, log(age)) grid: the closest node (m0, a0) to the
    model, with m0 in (ml, mh) and a0 in (al, ah). The IMF sample used is
    ``st_dist_mass[ml][al]``. Hence, for each IMF sample (ml, al) only the (up to)
    four nodes (ml or ml+1, al or al+1) are required.

    The indexes are ``np.searchsorted(mass_grid, st_dist_mass[ml][al])`` where
    ``mass_grid = theor_tracks[m0][a0][m_ini_idx]``, i.e.: the full mass grid
    before the maximum magnitude cut. These are stored with the smallest integer
    type that can hold them, to reduce memory usage. With numba the (identical)
    indexes are obtained by :py:func:`_mass_idx_nb`, using numba's threads.

    :param theor_tracks: Array of theoretical isochrones.
    :type theor_tracks: np.ndarray
    :param st_dist_mass: Sampled IMF mass distributions.
    :type st_dist_mass: list
    :param m_ini_idx: Index of the initial mass.
    :type m_ini_idx: int

    :returns: Dictionary with keys (ml, al, m0, a0) and integer arrays of indexes
     as values.
    :rtype: dict
    """
    N_z, N_a = theor_tracks.shape[:2]
    N_interp = theor_tracks.shape[-1]
    dtype = np.int16 if N_interp <= np.iinfo(np.int16).max else np.int32

    if not _HAS_NUMBA:
        lookup = {}
        for ml in range(N_z):
            for al in range(N_a):
                for m0 in sorted({ml, min(ml + 1, N_z - 1)}):
                    for a0 in sorted({al, min(al + 1, N_a - 1)}):
                        lookup[(ml, al, m0, a0)] = np.searchsorted(
                            theor_tracks[m0][a0][m_ini_idx], st_dist_mass[ml][al]
                        ).astype(dtype)
        return lookup

    from numba.typed.typedlist import List

    # Output arrays are allocated here and filled by the compiled function. They
    # are passed in a typed list where the indexes for the (ml, al) sample and the
    # node (m0, a0) are at position ((ml * N_a + al) * 2 + (m0 - ml)) * 2 + (a0 - al),
    # the same layout used by synth_batch.batch_data(). Nodes beyond the grid's
    # edges repeat the key (and the array) of the node inside the grid
    lookup, st_list, midx_list = {}, List(), List()
    for ml in range(N_z):
        for al in range(N_a):
            st_list.append(st_dist_mass[ml][al])
            for dz in (0, 1):
                for da in (0, 1):
                    key = (ml, al, min(ml + dz, N_z - 1), min(al + da, N_a - 1))
                    if key not in lookup:
                        lookup[key] = np.empty(st_dist_mass[ml][al].size, dtype)
                    midx_list.append(lookup[key])

    _mass_idx_nb(theor_tracks, m_ini_idx, st_list, midx_list, 4 * N_interp)

    return lookup


if _HAS_NUMBA:
    from numba import prange

    @njit(inline="always")
    def _sorted_idx_nb(g, lut, g0, inv_w, N_b, v):
        """Left insertion index of 'v' in the sorted array 'g', i.e.:
        ``np.searchsorted(g, v)``. The bucket table 'lut' only gives the
        starting point, the walks that follow make the result exact."""
        N = g.size
        f = (v - g0) * inv_w
        if f < 0.0:
            b = 0
        elif f >= N_b:
            b = N_b
        else:
            b = int(f)
        i = lut[b]
        while i < N and g[i] < v:
            i += 1
        while i > 0 and g[i - 1] >= v:
            i -= 1
        return i

    @njit(parallel=True, cache=True)
    def _mass_idx_nb(theor_tracks, m_ini_idx, st_list, midx_list, N_b):
        """Compiled (numba) version of :py:func:`mass_idx_lookup`, writing the
        indexes into the arrays of 'midx_list'.

        The loop runs over the nodes (m0, a0), distributed among numba's threads by
        metallicity. Each node's mass grid is searched by the (up to) four IMF
        samples (ml, al) with ml in (m0 - 1, m0) and al in (a0 - 1, a0). The search
        starts from a table of ``N_b + 1`` equal width buckets over the grid's mass
        range, built once per node, instead of a binary search over all the grid.
        """
        N_z, N_a = theor_tracks.shape[:2]
        N = theor_tracks.shape[-1]
        for m0 in prange(N_z):
            lut = np.empty(N_b + 1, np.int64)
            for a0 in range(N_a):
                g = theor_tracks[m0, a0, m_ini_idx]
                g0 = g[0]
                w = (g[N - 1] - g0) / N_b
                inv_w = 1.0 / w if w > 0.0 else 0.0
                # lut[b] = np.searchsorted(g, g0 + b * w)
                j = 0
                for b in range(N_b + 1):
                    edge = g0 + b * w
                    while j < N and g[j] < edge:
                        j += 1
                    lut[b] = j
                for ml in range(max(m0 - 1, 0), m0 + 1):
                    for al in range(max(a0 - 1, 0), a0 + 1):
                        s = ml * N_a + al
                        masses = st_list[s]
                        out = midx_list[(s * 2 + m0 - ml) * 2 + a0 - al]
                        for i in range(masses.size):
                            out[i] = _sorted_idx_nb(g, lut, g0, inv_w, N_b, masses[i])
