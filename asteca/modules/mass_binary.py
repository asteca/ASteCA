import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
from scipy.spatial import KDTree


def ranModels(fit_params, model_std, seed, N_models=200):
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


def get_m1m2_bpr(self, isoch, obs_phot):
    """ """
    # Masses
    mass_1, mass_2 = isoch[self.m_ini_idx], isoch[-1]

    # Indexes of closets synthetic stars to observed stars
    synth_photom = isoch[: self.m_ini_idx].T
    tree = KDTree(synth_photom)
    _, idxs = tree.query(obs_phot, k=1)
    # Assign primary and secondary (synthetic) masses to each observed star
    m1_obs, m2_obs = mass_1[idxs], mass_2[idxs]

    # Single systems are identified with m2=np.nan. This mask identifies observed
    # stars identified as binary systems
    m2_msk = ~np.isnan(m2_obs)
    # Binary fraction for the observed cluster
    b_fr = m2_msk.sum() / len(m2_obs)

    return m1_obs, m2_obs, b_fr


def get_masses(self, model, ra_c, dec_c, isoch, int_seed):
    """ """
    # Extract dm and loga
    model_comb = self.fix_params | model
    loga, dm = model_comb["loga"], model_comb["dm"]

    # N_obs = cl_dict["N_obs_stars"]
    mass_ini = isoch[self.m_ini_idx]
    M_obs = mass_ini.sum()
    mass_min, mass_max = mass_ini.min(), mass_ini.max()

    # Select a random IMF sampling array
    Nmets, Nages = len(self.st_dist_mass), len(self.st_dist_mass[0])
    # i = np.random.randint(Nmets)
    i = np.random.default_rng(self.seed + int_seed).integers(Nmets)
    # j = np.random.randint(Nages)
    j = np.random.default_rng(self.seed + int_seed).integers(Nages)
    mass_samples = self.st_dist_mass[i][j]

    mass_tot = np.cumsum(mass_samples)

    gamma = 0.62
    t = 10**loga
    qev_t = stellar_evol_mass_loss(loga)
    term1 = (1 - qev_t) ** gamma
    t0 = minit_LGB05(ra_c, dec_c, dm)

    def func_optm(N_max, flag_mass=False):
        M_init_sample = mass_samples[: int(N_max)]
        M_init = mass_tot[int(N_max)]

        # Masks for the min-max mass range
        msk_max = M_init_sample > mass_max
        M_ev = M_init_sample[msk_max].sum()
        M_actual_dyn = M_init - M_ev
        msk_min = M_init_sample < mass_min
        # This is the percentage of mass below the photometric minimum mass limit
        M_ratio = M_init_sample[msk_min].sum() / M_actual_dyn

        term2 = (gamma / (M_init**gamma)) * (t / t0)
        M_actual, M_actual_range = 0, 0
        if term1 > term2:
            M_actual = M_init * (term1 - term2) ** (1 / gamma)
            M_actual_range = M_actual - M_actual * M_ratio

        if flag_mass:
            return M_actual, M_init, abs(M_actual_range - M_obs)
        return abs(M_actual_range - M_obs)

    # Perform grid search
    optimal_param = grid_search_optimal_parameter(func_optm, 100, len(mass_samples))
    M_actual, M_init, mass_diff = func_optm(optimal_param, True)

    # print(round(loga, 2), int(M_actual), int(M_init), int(mass_diff))

    return M_init, M_actual


def grid_search_optimal_parameter(
    func, lower_bound, upper_bound, tolerance=500, max_iterations=5
) -> float:
    """
    Perform a grid search to find the optimal parameter within a given range.

    Parameters:
    - func: The objective function to optimize.
    - lower_bound: The lower bound of the parameter range.
    - upper_bound: The upper bound of the parameter range.
    - tolerance: The tolerance level to determine convergence.
    - max_iterations: Maximum number of iterations.

    Returns:
    - optimal_parameter: The optimal parameter within the specified range.
    """

    iteration = 0
    while iteration < max_iterations and (upper_bound - lower_bound) > tolerance:
        mid_point = (lower_bound + upper_bound) / 2
        par_range = upper_bound - lower_bound
        left_point = mid_point - par_range / 4
        right_point = mid_point + par_range / 4

        func_mid = func(mid_point)
        if func(left_point) < func_mid:
            upper_bound = mid_point
        elif func(right_point) < func_mid:
            lower_bound = mid_point
        else:
            lower_bound = left_point
            upper_bound = right_point

        iteration += 1

    optimal_parameter = (lower_bound + upper_bound) / 2

    return optimal_parameter


def stellar_evol_mass_loss(loga) -> float:
    """Fraction of the initial cluster mass (Mini) lost by stellar evolution"""
    a, b, c = 7, 0.26, -1.8
    q_ev = 10 ** ((max(7.1, loga) - a) ** b + c)
    return q_ev


def minit_LGB05(ra, dec, dm, epsilon=0.08):
    """
    Laplacian in spherical coordinates:

    https://planetmath.org/derivationofthelaplacianfromrectangulartosphericalcoordinates
    https://www.math.cmu.edu/~rcristof/pdf/Teaching/Spring2017/The%20Laplacian%20in%20polar%20and%20spherical%20coordinates(1).pdf

    I need polar coordinates in r so I just disregard the derivatives in the two
    angles.
    """
    # Constants for all clusters
    gamma = 0.62
    C_env0 = 810e6

    # Constants for MW potentials
    M_B = 2.5e10
    r_B = 0.5e3
    M_D = 7.5e10
    a = 5.4e3
    b = 0.3e3
    r_s = 15.19e3
    M_s = 1.87e11

    # Extract z value
    c = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)
    lon, lat = c.galactic.l, c.galactic.b
    dist_pc = 10 ** (0.2 * (dm + 5))
    c = SkyCoord(l=lon, b=lat, distance=dist_pc * u.pc, frame="galactic")
    c.representation_type = "cylindrical"
    z = c.z.value
    # Estimate R_GC
    gc = c.transform_to(coord.Galactocentric(galcen_distance=8 * u.kpc, z_sun=0 * u.pc))
    r = np.sqrt(gc.x.value**2 + gc.y.value**2)

    # Hernquist potential (bulge)
    # https://galaxiesbook.org/chapters/I-01.-Potential-Theory-and-Spherical-Mass-Distributions.html; Eq 3.58
    # This source above gives a density that does not match the two below formulas:
    # r_B/(r*(1+r/R_B)**3) <-- ?
    # https://docs.galpy.org/en/latest/reference/potentialhernquist.html
    # "note that amp is 2 x [total mass] for the chosen definition of the Two Power Spherical potential"
    # Phi_B_Laplacian = 2*M_B/r_B**3 * 1/((r/r_B)*(1+r/r_B)**3)
    # The above is equivalent to this one by
    # https://academic.oup.com/mnras/article/428/4/2805/992063; Eq 1
    Phi_B_Laplacian = 2 * M_B * r_B / (r * (r + r_B) ** 3)

    # Miyamoto & Nagai potential (disk)
    # https://galaxiesbook.org/chapters/II-01.-Flattened-Mass-Distributions.html#Thickened-disk:-the-Miyamoto-Nagai-model
    # https://articles.adsabs.harvard.edu/pdf/1975PASJ...27..533M
    # https://www.astro.utu.fi/~cflynn/galdyn/lecture4.html
    numerator = (
        M_D
        * b**2
        * (a * r**2 + (a + 3 * np.sqrt(z**2 + b**2)) * (a + np.sqrt(z**2 + b**2)) ** 2)
    )
    denominator = (b**2 + z**2) ** (3 / 2) * (
        r**2 + (a + np.sqrt(b**2 + z**2)) ** 2
    ) ** (5 / 2)
    Phi_D_laplacian = numerator / denominator

    # Sanderson potential (dark matter halo)
    A = -M_s / (np.log(2) - 0.5)
    Phi_H_Laplacian = -A / (r * (r + r_s) ** 2)

    # Ambient density
    rho_amb = (1 / (4 * np.pi)) * (Phi_B_Laplacian + Phi_D_laplacian + Phi_H_Laplacian)
    t0 = C_env0 * (1 - epsilon) * 10 ** (-4 * gamma) * rho_amb ** (-0.5)

    return t0
    # Minit = (M_actual**gamma + gamma*(t/t0))**(1/gamma)/(1-q_ev)
