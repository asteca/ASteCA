import astropy.coordinates as coord
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.spatial import KDTree


def ranModels(
    fit_params: dict,
    model_std: dict,
    N_models: int,
    rng: np.random.Generator,
) -> list[dict]:
    """Generate the 'N_models' models via sampling a Gaussian centered on
    'fit_params', with standard deviation given by 'model_std'.

    :param fit_params: Dictionary of parameters to fit.
    :type fit_params: dict
    :param model_std: Dictionary of standard deviations for each parameter.
    :type model_std: dict
    :param N_models: Number of models to generate.
    :type N_models: int
    :param rng: Random number generator.
    :type rng: np.random.Generator

    :return: List of dictionaries, each containing a set of sampled parameters.
    :rtype: list[dict]
    """
    models_ran = {}
    for k, f_val in fit_params.items():
        std = model_std[k]
        models_ran[k] = rng.normal(f_val, std, N_models)
    ran_models = [dict(zip(models_ran, t)) for t in zip(*models_ran.values())]

    return ran_models


def get_close_idxs(
    m_ini_idx: int, obs_phot: np.ndarray, isoch: np.ndarray
) -> np.ndarray:
    """Indexes of the closest synthetic stars to observed stars

    :param m_ini_idx: Index of the initial mass column
    :type m_ini_idx: int
    :param obs_phot: Observed photometry.
    :type obs_phot: np.ndarray
    :param isoch: Isochrone data.
    :type isoch: np.ndarray

    :return: Indexes of the closest synthetic stars.
    :rtype: np.ndarray
    """
    synth_photom = isoch[:m_ini_idx].T
    tree = KDTree(synth_photom)
    _, idxs = tree.query(obs_phot, k=1)
    return np.array(idxs)


def get_m1m2(
    m_ini_idx: int, isoch: np.ndarray, idxs: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Assign primary and secondary (synthetic) masses to each observed star.

    :param m_ini_idx: Index of the initial mass column
    :type m_ini_idx: int
    :param isoch: Isochrone data.
    :type isoch: np.ndarray
    :param idxs: Indexes of the closest synthetic stars.
    :type idxs: np.ndarray

    :return: Primary and secondary masses for each observed star.
    :rtype: tuple[np.ndarray, np.ndarray]
    """
    mass_1, mass_2 = isoch[m_ini_idx], isoch[-1]
    m1_obs, m2_obs = mass_1[idxs], mass_2[idxs]

    return m1_obs, m2_obs


def galactic_coords(
    sampled_models: list[dict],
    radec_c: tuple[float, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert equatorial coordinates to cylindrical, and obtain the vertical
    distance Z and the galactocentric distance R_GC.

    :param sampled_models: List of dictionaries, each containing a set of sampled
     parameters associated to non empty isochrones.
    :type sampled_models: list[dict]
    :param radec_c: Right ascension and declination of the cluster center.
    :type radec_c: tuple[float, float]

    :return: Vertical distance Z, galactocentric distance R_GC, and projected
        galactocentric distance R_xy.
    :rtype: tuple[np.ndarray, np.ndarray, np.ndarray]
    """
    c = SkyCoord(ra=radec_c[0] * u.degree, dec=radec_c[1] * u.degree)  # pyright: ignore
    lon, lat = c.galactic.l, c.galactic.b  # pyright: ignore
    dist_pc = []
    for model in sampled_models:
        dist_pc.append(10 ** (0.2 * (model["dm"] + 5)))
    cgal = SkyCoord(l=lon, b=lat, distance=dist_pc * u.pc, frame="galactic")  # pyright: ignore
    c_GC = cgal.transform_to(coord.Galactocentric())

    X, Y, Z = np.array(c_GC.x), np.array(c_GC.y), np.array(c_GC.z)
    R_GC = np.sqrt(X**2 + Y**2 + Z**2)
    R_xy = np.sqrt(X**2 + Y**2)

    return Z, R_GC, R_xy


def get_M_actual(
    rng: np.random.Generator,
    m_ini_idx: int,
    st_dist_mass: list[list],
    st_dist_mass_ordered: list[list],
    sampled_synthcl: np.ndarray,
) -> tuple[float, float, float]:
    """Estimate the actual mass using the observed mass and the fraction of
    mass estimated to be beyond the maximum observed magnitude.

    The final actual mass is just:

    M_actual = M_obs + M_phot

    :param rng: Random number generator
    :type rng: np.random.Generator
    :param m_ini_idx: Index of the initial mass column
    :type m_ini_idx: int
    :param st_dist_mass: List of sampled masses
    :type st_dist_mass: list[list]
    :param st_dist_mass_ordered: List of ordered sampled masses
    :type st_dist_mass_ordered: list[list]
    :param sampled_synthcl: Sampled synthetic cluster data.
    :type sampled_synthcl: np.ndarray

    :return: Observed, photometric, and actual mass.
    :rtype: tuple[float, float, float]
    """
    mass_ini = sampled_synthcl[m_ini_idx]
    mass_2nd = sampled_synthcl[-1]
    M_obs = np.nansum([mass_ini, mass_2nd])

    Nmets, Nages = len(st_dist_mass), len(st_dist_mass[0])
    i = rng.integers(Nmets)
    j = rng.integers(Nages)
    sorted_masses = st_dist_mass_ordered[i][j]

    idx_min = np.argmin(abs(sorted_masses - mass_ini.min()))
    idx_max = np.argmin(abs(sorted_masses - mass_ini.max()))
    M_phot_sample = max(1, sorted_masses[:idx_min].sum())
    M_obs_sample = max(1, sorted_masses[idx_min:idx_max].sum())
    factor = M_phot_sample / M_obs_sample
    M_phot = factor * M_obs

    M_a = M_obs + M_phot

    return M_obs, M_phot, M_a


def stellar_evol_mass_loss(z_met: float, loga: float) -> float:
    """Fraction of the initial cluster mass (M_ini) lost by stellar evolution.

    Source: Lamers, Baumgardt & Gieles (2010); Table B2
    (http://adsabs.harvard.edu/abs/2010MNRAS.409..305L)

    :param z_met: Metallicity.
    :type z_met: float
    :param loga: Logarithm of the age.
    :type loga: float

    :return: Fraction of mass lost by stellar evolution.
    :rtype: float
    """
    mu_coeffs = {
        "Z": np.array([0.0004, 0.0010, 0.0040, 0.0080, 0.0200]),
        "a0": np.array([1.0541, 1.0469, 1.0247, 1.0078, 0.9770]),
        "a1": np.array([-0.10912, -0.10122, -0.08307, -0.07456, -0.05709]),
        "a2": np.array([-0.01082, -0.01349, -0.01845, -0.02002, -0.02338]),
        "a3": np.array([0.00285, 0.00306, 0.00336, 0.00340, 0.00348]),
    }
    i = np.argmin(abs(mu_coeffs["Z"] - z_met))
    a0 = mu_coeffs["a0"][i]
    a1 = mu_coeffs["a1"][i]
    a2 = mu_coeffs["a2"][i]
    a3 = mu_coeffs["a3"][i]

    t_Myr = 10**loga
    x = np.log10(t_Myr / 1e6)
    mu_ev = a0 + a1 * x + a2 * x**2 + a3 * x**3

    return mu_ev


def ambient_density(
    M_B: float,
    r_B: float,
    M_D: float,
    a: float,
    b: float,
    r_s: float,
    M_s: float,
    Z: np.ndarray,
    R_GC: np.ndarray,
    R_xy: np.ndarray,
) -> np.ndarray:
    """Calculate the ambient density.

    Source: Angelo et al. (2023); 10.1093/mnras/stad1038

    The 4*pi constant is left for the final evaluation

    :param M_B: Bulge mass.
    :type M_B: float
    :param r_B: Bulge radius.
    :type r_B: float
    :param M_D: Disk mass.
    :type M_D: float
    :param a: Disk parameter a.
    :type a: float
    :param b: Disk parameter b.
    :type b: float
    :param r_s: Dark matter halo radius.
    :type r_s: float
    :param M_s: Dark matter halo mass.
    :type M_s: float
    :param Z: Vertical distance.
    :type Z: np.ndarray
    :param R_GC: Galactocentric distance.
    :type R_GC: np.ndarray
    :param R_xy: Projected galactocentric distance.
    :type R_xy: np.ndarray

    :return: Ambient density.
    :rtype: np.ndarray
    """
    Phi_B_Laplacian = 2 * M_B * r_B / (R_GC * (R_GC + r_B) ** 3)
    numerator = (
        M_D
        * b**2
        * (
            a * R_xy**2
            + (a + 3 * np.sqrt(Z**2 + b**2)) * (a + np.sqrt(Z**2 + b**2)) ** 2
        )
    )
    denominator = (b**2 + Z**2) ** (3 / 2) * (
        R_xy**2 + (a + np.sqrt(b**2 + Z**2)) ** 2
    ) ** (5 / 2)
    Phi_D_laplacian = numerator / denominator
    A = -M_s / (np.log(2) - 0.5)
    Phi_H_Laplacian = -A / (R_GC * (R_GC + r_s) ** 2)
    rho_amb = (1 / (4 * np.pi)) * (Phi_B_Laplacian + Phi_D_laplacian + Phi_H_Laplacian)

    return rho_amb


def dissolution_param(
    C_env: float, epsilon: float, gamma: float, rho_amb: float
) -> float:
    """Calculate the dissolution parameter.

    Lamers, Gieles & Zwart (2005), "Disruption time scales of star clusters in
    different galaxies" introduces the "disruption time" 't_dis' in Eq 8.

    The 'C_env' constant is originally set to 810 Myr and described as: "indicates
    the time when 95 per cent of the initial cluster mass is lost". In Eq 9 it is
    generalized to the range 300-800 Myr.

    In Eq 10 't_dis' is written as t_4*(Mi/10^4)^(gamma) where 't_4' is "the disruption
    time (in yrs) of a cluster with an initial mass of 10^4 Mo". This is the equivalent
    of 't0' in later articles.

    The relation with 'rho_amb' is described as "t_4 is expected to scale
    with the inverse square-root of the mean density in the host galaxy".
    For our Galaxy the authors estimate (see Table 1):
    t_4 ~ 10^8.75 ; rho_amb ~ 10^-1

    Lamers, Gieles, Bastian, Baumgardt, Kharchen & Zwart (2005), "An analytical
    description of the disruption of star clusters in tidal fields with an application
    to Galactic open clusters"
    The 't0' constant is "a constant that depends on the tidal field of the particular
    galaxy in which the cluster moves and on the ellipticity of its orbit". The
    relation with 'rho_amb' is described as: "'t0' is expected to depend on the
    ambient density at the location of the clusters in that galaxy as
    t0~rho_amb^(-1/2)". 'C_env' is said to be in the range 300-800 Myr.

    'mu_ev' is "fraction of the initial mass of the cluster that would have
    remained at age t, if stellar evolution would have been the only mass loss
    mechanism"

    Lamers, Baumgardt & Gieles (2010), "Mass-loss rates and the mass evolution of
    star clusters"

    "Theoretical considerations suggest that gamma ~ 0.65 to 0.85."

    "dissolution parameter (which is the hypothetical dissolution time-scale of a
    cluster of 1 M)"

    :param C_env: Constant related to the disruption time.
    :type C_env: float
    :param epsilon: Parameter related to the tidal field.
    :type epsilon: float
    :param gamma: Parameter related to the mass-loss rate.
    :type gamma: float
    :param rho_amb: Ambient density.
    :type rho_amb: float

    :return: Dissolution parameter.
    :rtype: float
    """
    t0 = C_env * (1 - epsilon) * 10 ** (-4 * gamma) * rho_amb ** (-0.5)

    return t0


def minit_LGB05(
    loga: float, M_actual: float, gamma: float, t0: float, mu_ev: float
) -> float:
    """Estimate the initial mass from Lamers et al. 2005.

    For loga>10 the M_init value grows *very* fast

    :param loga: Logarithm of the age.
    :type loga: float
    :param M_actual: Actual mass.
    :type M_actual: float
    :param gamma: Parameter related to the mass-loss rate.
    :type gamma: float
    :param t0: Dissolution parameter.
    :type t0: float
    :param mu_ev: Fraction of mass lost by stellar evolution.
    :type mu_ev: float

    :return: Initial mass.
    :rtype: float
    """
    t = 10**loga
    M_init = ((M_actual**gamma + gamma * (t / t0)) ** (1 / gamma)) / mu_ev

    return M_init
