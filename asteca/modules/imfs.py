import numpy as np
from scipy.interpolate import BSpline, make_interp_spline

from . import HAS_NUMBA as _HAS_NUMBA
from . import njit


def invTrnsfSmpl(IMF_name: str, m_low: float = 0.08, m_high: float = 100) -> BSpline:
    """IMF inverse transform sampling.

    :param IMF_name: Name of the IMF to use.
    :type IMF_name: str
    :param m_low: Lower mass limit, defaults to 0.08
    :type m_low: float
    :param m_high: Higher mass limit, defaults to 100
    :type m_high: float

    :return: Interpolated function of the inverse cumulative distribution function.
    :rtype: BSpline
    """

    # The lower mass region needs to be sampled more accurately
    mass_values = list(np.linspace(m_low, 0.5, 100))
    mass_values += list(np.linspace(0.501, 2, 75))
    mass_values += list(np.linspace(2.01, 10, 50))
    mass_values += list(np.linspace(10.01, m_high, 25))
    mass_values = np.array(mass_values)

    # DEPRECATED 12/01/25
    #
    # IMF_old = get_imf(IMF_name, m_low)
    # CDF_samples, m_old, area_CDF = [], m_low, 0.0
    # for m in mass_values[1:]:
    #     # Approximate integral with rectangular area, and add to previous total area
    #     IMF_new = get_imf(IMF_name, m)
    #     area_CDF += 0.5 * (IMF_new + IMF_old) * (m - m_old)
    #     CDF_samples.append(area_CDF)
    #     IMF_old, m_old = IMF_new, m
    # CDF_samples = np.array(CDF_samples) / max(CDF_samples)

    # Approximate integral with rectangular area, and add to previous total area
    # This is the length of the base of the rectangle
    length = mass_values[1:] - mass_values[:-1]
    # This is the average height of the rectangle
    imf_val = get_imf(IMF_name, mass_values)
    height = 0.5 * (imf_val[1:] + imf_val[:-1])
    # Rectangle area
    r_area = length * height
    # Cumulative sum
    CDF_samples = np.cumsum(r_area)
    # Normalize
    CDF_samples = CDF_samples / max(CDF_samples)

    # k=1 is important otherwise the interpolator can become unstable for values
    # close to 1 and return huge masses and even negative ones
    inv_cdf = make_interp_spline(CDF_samples, mass_values[1:], k=1)

    return inv_cdf


def get_imf(IMF_name: str, m_star_array: np.ndarray) -> np.ndarray:
    """Define any number of IMFs.
    The package https://github.com/keflavich/imf has some more (I think).

    :param IMF_name: Name of the IMF to use.
    :type IMF_name: str
    :param m_star_array: Array of stellar masses.
    :type m_star_array: np.ndarray

    :raises ValueError: If the IMF is not implemented.

    :return: Array of IMF values.
    :rtype: np.ndarray
    """
    if IMF_name == "salpeter_1955":
        # Salpeter (1955)  IMF.
        # https://ui.adsabs.harvard.edu/abs/1955ApJ...121..161S/
        imf_vals = m_star_array**-2.35

    elif IMF_name == "kroupa_2001":
        # Kroupa (2001), Mon. Not. R. Astron. Soc. 322, 231-246 (2001); Eq. 2
        # https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..231K/
        alpha = [-0.3, -1.3, -2.3]
        m0, m1, m2 = [0.01, 0.08, 0.5]
        # Continuity factor taken from Popescu & Hanson (2009; MASSCLEAN),
        # Eq. (2) & (3), p. 1725
        factor = [
            (1.0 / m1) ** alpha[0],
            (1.0 / m1) ** alpha[1],
            ((m2 / m1) ** alpha[1]) * ((1.0 / m2) ** alpha[2]),
        ]

        # Conditions for m_star_array <= m0 and m_star_array > m0
        lower_msk = (m0 <= m_star_array) & (m_star_array <= m1)
        middle_msk = (m1 < m_star_array) & (m_star_array <= m2)
        upper_msk = m_star_array > m2

        # Initialize output array
        imf_vals = np.zeros_like(m_star_array)
        imf_vals[lower_msk] = factor[0] * (m_star_array[lower_msk] ** alpha[0])
        imf_vals[middle_msk] = factor[1] * (m_star_array[middle_msk] ** alpha[1])
        imf_vals[upper_msk] = factor[2] * (m_star_array[upper_msk] ** alpha[2])

    elif IMF_name == "chabrier_2014":
        # Chabrier et al. (2014)
        # https://ui.adsabs.harvard.edu/abs/2014ApJ...796...75C/ ; Eq (34)
        nc, mc = 11, 0.18
        m0 = nc * mc
        Ah, x = 0.649, 1.35
        Al = Ah * nc ** (x / 2)
        sigma_2 = np.log10(nc) / (x * np.log(10))
        c_array = 0.434294 / m_star_array  # np.log10(e)/m ; This is the transformation
        # from dN/log(m) --> dN/dm

        # Conditions for m_star_array <= m0 and m_star_array > m0
        lower_msk = m_star_array <= m0
        upper_msk = m_star_array > m0

        # Initialize output array
        imf_vals = np.zeros_like(m_star_array)
        # Compute imf_val for m_star <= m0
        imf_vals[lower_msk] = (
            c_array[lower_msk]
            * Al
            * m0 ** (-x)
            * np.exp(
                -((np.log10(m_star_array[lower_msk]) - np.log10(mc)) ** 2)
                / (2 * sigma_2)
            )
        )
        # Compute imf_val for m_star > m0
        imf_vals[upper_msk] = c_array[upper_msk] * Ah * m_star_array[upper_msk] ** (-x)
    else:
        raise ValueError(f"IMF {IMF_name} not implemented.")

    return imf_vals


def sampleInv(
    rng: np.random.Generator,
    Max_mass: float,
    inv_cdf: BSpline,
    N_chunk: int,
    N_draw: int,
    lut: np.ndarray | None = None,
    buf: np.ndarray | None = None,
) -> np.ndarray:
    """Sample the inverse CDF up to `Max_mass`

    Masses are sampled in chunks of `N_chunk` stars until the total mass reaches
    `Max_mass`. To avoid sampling one chunk at a time, `N_draw` uniform values are
    drawn at once, and the masses are kept up to the end of the first chunk that
    reaches `Max_mass`. The random generator is then rewound and advanced by the
    number of values kept, so that the results (and the state of `rng`) are the
    same as those obtained sampling one chunk at a time.

    :param rng: Random number generator.
    :type rng: np.random.Generator
    :param Max_mass: Maximum mass to sample.
    :type Max_mass: float
    :param inv_cdf: Interpolated function of the inverse cumulative distribution function.
    :type inv_cdf: BSpline
    :param N_chunk: Number of stars sampled in each chunk
    :type N_chunk: int
    :param N_draw: Number of uniform values drawn at once (a multiple of `N_chunk`).
        It is doubled if it is not enough to reach `Max_mass`
    :type N_draw: int
    :param lut: Bucket lookup table of the inverse CDF knots, as returned by
        :py:func:`invCDF_lut` (numba only). Computed if not given; pass it to
        avoid re-computing it when this function is called repeatedly.
    :type lut: np.ndarray | None
    :param buf: Work array of at least `N_draw` elements, re-used to hold the
        uniform values (avoids allocating a new array in each call).
    :type buf: np.ndarray | None

    :return: Array of sampled masses.
    :rtype: np.ndarray
    """
    if _HAS_NUMBA:
        # The spline is of order k=1 (linear) so its knots without the
        # boundary repetitions and its coefficients are the (x, y) points
        x, y = inv_cdf.t[1:-1], inv_cdf.c
        if lut is None:
            lut = invCDF_lut(x)

    state = rng.bit_generator.state
    while True:
        if buf is None or buf.size < N_draw:
            buf = np.empty(N_draw)
        # Same values as rng.uniform(0.0, 1.0, N_draw)
        u_vals = rng.random(out=buf[:N_draw])
        if _HAS_NUMBA:
            N_kept = _sample_chunks_nb(u_vals, x, y, lut, Max_mass, N_chunk)  # pyright: ignore[reportPossiblyUnboundVariable]
            masses = u_vals
        else:
            masses = inv_cdf(u_vals)
            # Total mass at the end of each chunk
            chunk_mass = np.cumsum(masses.reshape(-1, N_chunk).sum(axis=1))
            N_kept = (int(np.searchsorted(chunk_mass, Max_mass)) + 1) * N_chunk
            if N_kept > N_draw:
                N_kept = -1
        if N_kept > 0:
            break
        # Not enough values drawn to reach 'Max_mass': rewind and draw more
        rng.bit_generator.state = state
        N_draw *= 2

    # Rewind and advance the generator by the number of values kept. Each
    # uniform float consumes a single 64 bit value of the generator
    rng.bit_generator.state = state
    # Not all bit generators implement 'advance()' (e.g.: PCG64 and Philox do)
    advance = getattr(rng.bit_generator, "advance", None)
    if advance is not None:
        advance(N_kept)
    else:
        rng.uniform(0.0, 1.0, N_kept)

    # Copy to release the values not kept (and the work array)
    return masses[:N_kept].copy()


def invCDF_lut(x: np.ndarray, N_buckets: int = 4096) -> np.ndarray:
    """Bucket lookup table for the knots of the inverse CDF.

    The [0, 1] range is split into `N_buckets` equal buckets, and for each bucket
    edge ``b / N_buckets`` the table stores the largest index ``k`` such that
    ``x[k] < b / N_buckets`` (i.e.: ``np.searchsorted(x, b / N_buckets) - 1``).
    This gives a starting point for the knot search of a uniform value that is
    (almost always) at most a few knots away from the result, replacing the
    binary search over all the knots.

    :param x: Sorted knots of the inverse CDF (CDF values).
    :type x: np.ndarray
    :param N_buckets: Number of buckets, defaults to 4096
    :type N_buckets: int

    :return: Lookup table of `N_buckets + 1` indexes.
    :rtype: np.ndarray
    """
    edges = np.arange(N_buckets + 1) / N_buckets
    return (np.searchsorted(x, edges) - 1).astype(np.int64)


if _HAS_NUMBA:

    @njit(cache=True)
    def _sample_chunks_nb(
        u_vals: np.ndarray,
        x: np.ndarray,
        y: np.ndarray,
        lut: np.ndarray,
        Max_mass: float,
        N_chunk: int,
    ) -> int:
        """Compiled (numba) version of the chunk sampling in :py:func:`sampleInv`.

        The uniform values in 'u_vals' are replaced in place by their masses,
        linearly interpolated (and extrapolated, as the k=1 BSpline does) in the
        (x, y) points of the inverse CDF. The masses are summed per chunk, as
        the original one chunk at a time loop did, and the conversion stops at
        the end of the first chunk that reaches 'Max_mass'.

        The knot of each value is located starting from the bucket lookup table
        'lut' (see :py:func:`invCDF_lut`) and walking to the largest 'k' with
        x[k] < u, the same index given by ``np.searchsorted(x, u) - 1``. The
        backward walk only guards against the rounding of ``u * N_buckets``.

        Returns the number of masses converted, or -1 if 'Max_mass' was not
        reached.
        """
        N_x = x.size
        N_b = lut.size - 1
        Mass_tot = 0.0
        i = 0
        while i + N_chunk <= u_vals.size:
            chunk_mass = 0.0
            for _ in range(N_chunk):
                u = u_vals[i]
                k = lut[min(int(u * N_b), N_b)]
                while k + 1 < N_x and x[k + 1] < u:
                    k += 1
                while k >= 0 and x[k] >= u:
                    k -= 1
                k = min(max(k, 0), N_x - 2)
                m = y[k] + (u - x[k]) * (y[k + 1] - y[k]) / (x[k + 1] - x[k])
                u_vals[i] = m
                chunk_mass += m
                i += 1
            Mass_tot += chunk_mass
            if Mass_tot >= Max_mass:
                return i
        return -1
