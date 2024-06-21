import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
from scipy.spatial import KDTree


def ranModels(fit_params: dict, model_std: dict, N_models: int, seed: int) -> list:
    """
    Generate the 'N_models' models via sampling a Gaussian centered on 'fit_params',
    with standard deviation given by 'model_std'.
    """
    models_ran, int_seed = {}, np.random.default_rng(seed).integers(1)
    for k, f_val in fit_params.items():
        std = model_std[k]
        # models_ran[k] = np.random.normal(f_val, std, N_models)
        models_ran[k] = np.random.default_rng(seed + int_seed).normal(
            f_val, std, N_models
        )
        int_seed += 1
    # Transpose dict of arrays into list of dicts
    ran_models = [dict(zip(models_ran, t)) for t in zip(*models_ran.values())]

    return ran_models


def get_close_idxs(self, obs_phot, isoch):
    """Indexes of the closest synthetic stars to observed stars"""
    synth_photom = isoch[: self.m_ini_idx].T
    tree = KDTree(synth_photom)
    _, idxs = tree.query(obs_phot, k=1)
    return idxs


def get_m1m2(self, isoch: np.array, idxs: np.array):
    """ """
    # Masses
    mass_1, mass_2 = isoch[self.m_ini_idx], isoch[-1]

    # Assign primary and secondary (synthetic) masses to each observed star
    m1_obs, m2_obs = mass_1[idxs], mass_2[idxs]

    return m1_obs, m2_obs


def get_bpr(self, isoch: np.array, idxs: np.array):
    """ """
    # Secondary masses
    mass_2 = isoch[-1]

    # Assign secondary (synthetic) masses to each observed star
    m2_obs = mass_2[idxs]

    # Single systems are identified with m2=np.nan. This mask identifies observed
    # stars identified as binary systems
    m2_msk = ~np.isnan(m2_obs)

    # Total binary fraction for the observed cluster
    b_fr = m2_msk.sum() / len(m2_obs)

    return b_fr


def galactic_coords(
    synthcl,
    radec_c: float,
) -> tuple[np.array, np.array]:
    """Convert equatorial coordinates to cylindrical, and obtain the vertical distance
    Z and the galactocentric distance R_GC
    """
    c = SkyCoord(ra=radec_c[0] * u.degree, dec=radec_c[1] * u.degree)
    lon, lat = c.galactic.l, c.galactic.b
    dist_pc = []
    for model in synthcl.sampled_models:
        # Extract dm
        model_comb = synthcl.fix_params | model
        dist_pc.append(10 ** (0.2 * (model_comb["dm"] + 5)))
    cgal = SkyCoord(l=lon, b=lat, distance=dist_pc * u.pc, frame="galactic")
    c_GC = cgal.transform_to(coord.Galactocentric())

    X, Y, Z = c_GC.x.value, c_GC.y.value, c_GC.z.value
    R_GC = np.sqrt(X**2 + Y**2 + Z**2)
    R_xy = np.sqrt(X**2 + Y**2)

    return Z, R_GC, R_xy


def get_M_actual(synthcl, isoch, int_seed) -> tuple[float, float]:
    """Estimate the actual mass using the observed mass and the fraction of
    mass estimated to be beyond the maximum observed magnitude.
    """

    mass_ini = isoch[synthcl.m_ini_idx]
    mass_2nd = isoch[-1]
    # Add secondary masses
    M_obs = np.nansum([mass_ini, mass_2nd])

    # Select a random IMF sampling array (faster than sampling the IMF in place)
    Nmets, Nages = len(synthcl.st_dist_mass), len(synthcl.st_dist_mass[0])
    i = np.random.default_rng(synthcl.seed + int_seed).integers(Nmets)
    j = np.random.default_rng(synthcl.seed + int_seed).integers(Nages)
    sorted_masses = synthcl.st_dist_mass_ordered[i][j]

    idx_min = np.argmin(abs(sorted_masses - mass_ini.min()))
    idx_max = np.argmin(abs(sorted_masses - mass_ini.max()))
    M_phot_sample = sorted_masses[:idx_min].sum()
    M_obs_sample = sorted_masses[idx_min:idx_max].sum()

    # This is the ratio of the sampled mass below the minimum mass value,
    # over the sampled mass within the observed mass range
    factor = M_phot_sample / M_obs_sample

    # This is the 'photometric mass', or the mass that is lost beyond the minimum
    # observed mass
    M_phot = factor * M_obs

    return M_obs, M_phot


def stellar_evol_mass_loss(z_met, loga) -> float:
    """Fraction of the initial cluster mass (M_ini) lost by stellar evolution

    Source: Lamers, Baumgardt & Gieles (2010); Table B2
    (http://adsabs.harvard.edu/abs/2010MNRAS.409..305L)
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


def ambient_density(M_B, r_B, M_D, a, b, r_s, M_s, Z, R_GC, R_xy):
    """
    Source: Angelo et al. (2023); 10.1093/mnras/stad1038

    The 4*pi constant is left for the final evaluation
    """

    # Hernquist potential (bulge)
    # https://galaxiesbook.org/chapters/I-01.-Potential-Theory-and-Spherical-Mass-Distributions.html; Eq 3.58
    # Phi_B_Laplacian = rho_0 * r_B/(r*(1+r/r_B)**3)
    # https://docs.galpy.org/en/latest/reference/potentialhernquist.html
    # "note that amp is 2 x [total mass] for the chosen definition of the Two Power Spherical potential"
    # Phi_B_Laplacian = 2*M_B/(4*pi*r_B**3) * 1/((r/r_B)*(1+r/r_B)**3)
    # The two above are equivalent if rho_0 = 2*M_B/(4*pi*r_B**3)
    Phi_B_Laplacian = 2 * M_B * r_B / (R_GC * (R_GC + r_B) ** 3)

    # Miyamoto & Nagai potential (disk)
    # https://galaxiesbook.org/chapters/II-01.-Flattened-Mass-Distributions.html#Thickened-disk:-the-Miyamoto-Nagai-model
    # https://articles.adsabs.harvard.edu/pdf/1975PASJ...27..533M
    # https://www.astro.utu.fi/~cflynn/galdyn/lecture4.html
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

    # Sanderson, Hartke & Helmi (2017) potential (dark matter halo)
    A = -M_s / (np.log(2) - 0.5)
    Phi_H_Laplacian = -A / (R_GC * (R_GC + r_s) ** 2)

    # Ambient density
    rho_amb = (1 / (4 * np.pi)) * (Phi_B_Laplacian + Phi_D_laplacian + Phi_H_Laplacian)

    return rho_amb


def dissolution_param(C_env, epsilon, gamma, rho_amb):
    """
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


    """
    # Dissolution parameter
    t0 = C_env * (1 - epsilon) * 10 ** (-4 * gamma) * rho_amb ** (-0.5)

    return t0


def minit_LGB05(loga, M_actual, gamma, t0, mu_ev):
    """
    Initial mass estimation from Lamers et al. 2005
    """
    t = 10**loga
    M_init = (M_actual**gamma + gamma * (t / t0)) ** (1 / gamma) / mu_ev

    return M_init
