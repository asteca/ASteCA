import numpy as np
from scipy.interpolate import make_interp_spline, BSpline


def invTrnsfSmpl(
    IMF_name: str, m_low: float = 0.08, m_high: float = 100
) -> BSpline:
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
) -> np.ndarray:
    """Sample the inverse CDF up to `Max_mass`

    :param rng: Random number generator.
    :type rng: np.random.Generator
    :param Max_mass: Maximum mass to sample.
    :type Max_mass: float
    :param inv_cdf: Interpolated function of the inverse cumulative distribution function.
    :type inv_cdf: BSpline
    :param N_chunk: Integer that controls the while loop that samples masses in chunks
    :type N_chunk: int

    :return: Array of sampled masses.
    :rtype: np.ndarray
    """

    # DEPRECATED 13/01/25
    # def sampled_inv_cdf(N):
    #     # The internal seed can not be 'self.seed' because otherwise all the
    #     # floats returned are equal
    #     # mr = np.random.default_rng(ijkseed).uniform(0.0, 1.0, N)
    #     return inv_cdf(mr)

    # TODO: make this faster?
    k, Mass_tot, mass_samples = 0, 0, []
    while Mass_tot < Max_mass:
        masses = list(inv_cdf(rng.uniform(0.0, 1.0, N_chunk)))
        mass_samples += masses
        Mass_tot += sum(masses)
        k += 1

    return np.array(mass_samples)
