import numpy as np
import numpy.typing as npt
from scipy import stats

from .imfs import imf_inverse_cdf

# Number of candidate stars available to each synthetic cluster, per observed star,
# plus a fixed number (see n_candidates()). A synthetic cluster can only have fewer
# stars than the observed one if fewer than 1 in IMF_N_PER_STAR candidates remain
# brighter than the maximum magnitude after the extinction
IMF_N_PER_STAR = 20
IMF_N_EXTRA = 1000
# Number of buckets of the lookup table of the IMF CDF of each grid node (see
# imf_cdf_tables())
IMF_LUT_BUCKETS = 4096
# Differential reddening models, see dr_fractions()
DR_MODELS = ("foreground", "internal")
# Mean of the differential reddening fractions of each model (see dr_fractions()):
# 1/2 for the normal truncated at +/-3 sigma around 0.5, and 3/8 for the density
# 3/2 (1 - f^2). The mean absorption of the stars is Av + DR_MEAN_FRAC * DR
DR_MEAN_FRAC = {"foreground": 0.5, "internal": 0.375}
# Number of quantiles of the differential reddening fractions stored by
# draw_rand_floats() (see generate.lost_fraction())
DR_FRAC_N_QUANTILES = 65


def n_candidates(N_stars: int, N_isoch: int) -> int:
    """Number of uniform random values drawn for the candidate stars of each
    synthetic cluster (and for their binarity and differential reddening), for
    synthetic clusters of ``N_stars`` stars: ``IMF_N_PER_STAR * N_stars +
    IMF_N_EXTRA``, and at least the number of isochrone points

    :param N_stars: Number of synthetic stars
    :type N_stars: int
    :param N_isoch: Number of points of the isochrones
    :type N_isoch: int

    :returns: Number of candidates
    :rtype: int
    """
    return max(IMF_N_PER_STAR * int(N_stars) + IMF_N_EXTRA, int(N_isoch))


def _stream(seed: int, *key: int) -> np.random.Generator:
    """Independent random generator for the random values identified by ``key``,
    obtained from ``seed``. The values drawn from it do not depend on how many are
    drawn, nor on any other random values, which makes the random values of the
    candidate stars nested: those obtained for a larger number of candidates start
    with those obtained for a smaller one (see :py:func:`draw_rand_floats`)

    :param seed: Random seed
    :type seed: int
    :param *key: Integers that identify the stream
    :type *key: int

    :returns: Random number generator
    :rtype: np.random.Generator
    """
    return np.random.default_rng(np.random.SeedSequence(seed, spawn_key=key))


def add_binary_tracks(
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
            q = mass_ratios(mass_ini, gamma, rng)
            # Secondary masses
            m2 = q * mass_ini

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


def mass_ratios(
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

        # Use 'gamma + 1' in the power-law distribution because in D&K this
        # distribution is defined as f(q)~q^gamma, while numpy's distribution is
        # defined as a*x^(a-1).
        ratios = rng.power(gamma + 1, N)

    except ValueError:
        if gamma == "D&K":
            msk1, gamma1 = M1 <= 0.1, 4.2
            msk2, gamma2 = (M1 > 0.1) & (M1 <= 0.6), 0.4
            msk3, gamma3 = (M1 > 0.6) & (M1 <= 1.4), 0.3
            msk4, gamma4 = (M1 > 1.4) & (M1 <= 6.5), -0.5
            msk5, gamma5 = (M1 > 6.5) & (M1 <= 16), 0.0  # <- Not sure. Use uniform
            msk6, gamma6 = M1 > 16, 0.0  # <- Not sure. Use uniform

            ratios = np.zeros(N)
            for msk, gammaX in (
                (msk1, gamma1),
                (msk2, gamma2),
                (msk3, gamma3),
                (msk4, gamma4),
                (msk5, gamma5),
                (msk6, gamma6),
            ):
                q = rng.power(gammaX + 1, msk.sum())
                ratios[msk] = q
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
            ratios = np.interp(u, cdf, q_edges)

    return ratios


def ccmo_filter_coefs(
    magnitude_effl: float,
    color_effl: tuple,
    color2_effl: tuple | None,
) -> np.ndarray:
    """Obtain extinction coefficients for all the observed filters and colors,
    in the order in which they are stored in theor_tracks.

    The coefficient of a filter is ``a + b / Rv`` (see :py:func:`ccmo_model`). The
    array has a row ``[a, b, 0, 0]`` for the magnitude, and a row
    ``[a1, b1, a2, b2]`` for each color (one pair for each of its two filters), see
    :py:func:`generate.ccmo_ratios`.

    :param magnitude_effl: Effective lambda (in Angstrom) for the magnitude filter
    :type magnitude_effl: float
    :param color_effl: Effective lambdas for the filters that make up the first color
    :type color_effl: tuple
    :param color2_effl: Effective lambdas for the filters that make up the second color
    :type color2_effl: tuple | None

    :returns: Array of extinction coefficients, shape ``(N_colors + 1, 4)``
    :rtype: np.ndarray

    """
    # Effective wavelengths in Angstrom, converted to inverse microns
    rows = [[*ccmo_model(10000.0 / magnitude_effl), 0.0, 0.0]]
    for effl in (color_effl, color2_effl):
        if effl is not None:
            rows.append(
                [*ccmo_model(10000.0 / effl[0]), *ccmo_model(10000.0 / effl[1])]
            )
    return np.array(rows, dtype=float)


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


def imf_cdf_tables(
    IMF_name: str, theor_tracks: np.ndarray, m_ini_idx: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """IMF data used to draw the masses of the candidate stars of the synthetic
    clusters (see :py:func:`generate.draw_candidates`).

    The inverse CDF of the IMF (:py:func:`imfs.imf_inverse_cdf`) is linear between its
    knots ``(x, y)`` (CDF values, masses). The CDF ``F`` is evaluated at the initial
    masses of every grid node of the isochrones (``np.interp(mass, y, x)``), and a
    lookup table of ``IMF_LUT_BUCKETS + 1`` equal buckets over ``[0, 1]`` stores
    ``np.searchsorted(F, b / IMF_LUT_BUCKETS)`` for each node, the starting point of
    the search of a CDF value in the compiled functions.

    :param IMF_name: Name of the IMF
    :type IMF_name: str
    :param theor_tracks: Theoretical isochrones, shape ``(Nz, Na, N_cols, N_pts)``
    :type theor_tracks: np.ndarray
    :param m_ini_idx: Index of the initial mass row
    :type m_ini_idx: int

    :returns: Knots of the inverse CDF ``(x, y)``, the CDF at the masses of each
        grid node (shape ``(Nz, Na, N_pts)``), and the lookup tables (shape
        ``(Nz, Na, IMF_LUT_BUCKETS + 1)``)
    :rtype: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
    """
    inv_cdf = imf_inverse_cdf(IMF_name)
    x = np.ascontiguousarray(inv_cdf.t[1:-1], dtype=float)
    y = np.ascontiguousarray(inv_cdf.c, dtype=float)

    Nz, Na = theor_tracks.shape[:2]
    N_pts = theor_tracks.shape[-1]
    dtype = np.int16 if N_pts <= np.iinfo(np.int16).max else np.int32
    edges = np.arange(IMF_LUT_BUCKETS + 1) / IMF_LUT_BUCKETS
    imf_cdf = np.empty((Nz, Na, N_pts))
    imf_lut = np.empty((Nz, Na, IMF_LUT_BUCKETS + 1), dtype=dtype)
    for i in range(Nz):
        for j in range(Na):
            F = np.interp(theor_tracks[i, j, m_ini_idx], y, x)
            imf_cdf[i, j] = F
            imf_lut[i, j] = np.searchsorted(F, edges)
    return x, y, imf_cdf, imf_lut


def draw_rand_floats(seed: int, N_cand: int, N_nodes: int) -> dict:
    """Random values used to generate the synthetic clusters, for ``N_cand``
    candidate stars (see :py:func:`n_candidates`).

    Each row of values is drawn from its own random stream (see :py:func:`_stream`),
    so that the values obtained for a larger ``N_cand`` start with those obtained
    for a smaller one (except ``imf_offsets``, drawn in ``[0, N_cand)``).

    :param seed: Random seed.
    :type seed: int
    :param N_cand: Number of candidate stars
    :type N_cand: int
    :param N_nodes: Number of grid nodes of the isochrones
    :type N_nodes: int

    :returns: Dictionary of random values: ``norm`` (error distribution, see
        :py:func:`error_distribution`), ``unif`` (binarity and differential
        reddening), ``dr_frac`` (a dictionary with the fractions of the differential
        reddening for each model in ``DR_MODELS``, all obtained from the same uniform
        values), ``dr_frac_q`` (``DR_FRAC_N_QUANTILES`` quantiles of each ``dr_frac``
        array, from ``0`` to ``1``), ``imf_u`` (``2 * N_cand`` uniform values for the
        masses of the candidates) and ``imf_offsets`` (the start of the window of
        ``N_cand`` values of ``imf_u`` used by each grid node)
    :rtype: dict
    """
    # Row 1: used by error_distribution() (Synthetic.calibrate()). Row 0 is not
    # used (it keeps the layout of the array)
    rand_norm_vals = np.empty((2, N_cand))
    # Row 0: used by dr_fractions() (below). Row 1: used by assign_binaries()
    rand_unif_vals = np.empty((2, N_cand))
    for row in (0, 1):
        _stream(seed, 1, row).standard_normal(out=rand_norm_vals[row])
        _stream(seed, 2, row).random(out=rand_unif_vals[row])

    # Used by `apply_extinction()`. The fractions of every differential reddening model
    # are obtained (from the same uniform values) so that the model can be selected
    # when a synthetic cluster is generated, at no cost
    dr_frac = {m: dr_fractions(rand_unif_vals[0], m) for m in DR_MODELS}
    # Quantiles of the fractions, used to estimate how many candidate stars are lost
    # to the differential reddening (see generate.lost_fraction())
    q_levels = np.linspace(0.0, 1.0, DR_FRAC_N_QUANTILES)
    dr_frac_q = {m: np.quantile(v, q_levels) for m, v in dr_frac.items()}

    # Uniform values for the masses of the candidate stars. Each grid node (the
    # (z, a) box of the model) uses its own window of N_cand of them, at a random
    # start, so that the masses of the synthetic clusters of different nodes are
    # (nearly) independent, as with an independent IMF sample per node
    imf_u = _stream(seed, 3).random(2 * N_cand)
    imf_offsets = _stream(seed, 4).integers(0, N_cand, N_nodes)

    rand_floats = {
        "norm": rand_norm_vals,
        "unif": rand_unif_vals,
        "dr_frac": dr_frac,
        "dr_frac_q": dr_frac_q,
        "imf_u": imf_u,
        "imf_offsets": imf_offsets,
    }
    return rand_floats


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
    values at the limits of the range. If ``low = high`` all the values are equal to
    it.

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
    if low == high:
        return np.full(N, float(low))
    if std <= 0.0:
        return np.full(N, float(mean))
    a, b = (low - mean) / std, (high - mean) / std
    return stats.truncnorm.rvs(a, b, loc=mean, scale=std, size=N, random_state=rng)


def sample_params(
    rng: np.random.Generator,
    model: dict[str, float],
    model_std: dict[str, float],
    limits: dict[str, tuple[float, float]],
    N_models: int,
) -> list[dict[str, float]]:
    """Sample sets of parameters from normal distributions centered on ``model``
    with STDDEVs ``model_std``, truncated to the valid range of each parameter (see
    :py:func:`truncnorm_sample`).

    :param rng: Random number generator
    :type rng: np.random.Generator
    :param model: Dictionary with the values of all the parameters
    :type model: dict[str, float]
    :param model_std: Dictionary with the STDDEV of each parameter in ``model``
    :type model_std: dict[str, float]
    :param limits: Valid range of each parameter, see
        :py:attr:`SynthEngine.sampling_ranges <engine.SynthEngine.sampling_ranges>`
    :type limits: dict[str, tuple[float, float]]
    :param N_models: Number of sets of parameters
    :type N_models: int

    :return: List of ``N_models`` dictionaries, one per set of parameters
    :rtype: list[dict[str, float]]
    """
    sampled = {}
    for k, f_val in model.items():
        low, high = limits[k]
        sampled[k] = truncnorm_sample(rng, f_val, model_std[k], low, high, N_models)
    return [dict(zip(sampled, t)) for t in zip(*sampled.values())]


def dr_fractions(rand_unif: np.ndarray, DR_distribution: str) -> np.ndarray:
    """Fractions 'f' (in [0, 1]) of the differential reddening applied to each
    synthetic star, A_i = Av + f_i * DR.

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
) -> np.ndarray:
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

    :returns: Error distributions, one row per magnitude and color
    :rtype: np.ndarray
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

    return np.array(err_dist)
