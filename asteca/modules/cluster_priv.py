import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
from scipy.spatial import KDTree


def load(self):
    """
    The photometry is store with a '_p' to differentiate from the self.magnitude,
    self.color, etc that are defined with the class is called.
    """
    print("Reading and processing cluster data")

    self.mag_p = np.array(self.cluster_df[self.magnitude])
    self.e_mag_p = np.array(self.cluster_df[self.e_mag])
    
    self.colors_p = [np.array(self.cluster_df[self.color])]
    if self.color2 is not None:
        self.colors_p.append(np.array(self.cluster_df[self.color2]))
    self.e_colors_p = [np.array(self.cluster_df[self.e_color])]
    if self.e_color2 is not None:
        self.e_colors_p.append(np.array(self.cluster_df[self.e_color2]))


def ranModels(N_models, fit_params, model_std):
    """
    Generate the requested models via sampling a Gaussian centered on the
    selected solution, with standard deviation given by the attached
    uncertainty.

    N_models: number of models to generate (HARDCODED)
    """
    models_ran = {}
    for k, f_val in fit_params.items():
        std = model_std[k]
        models_ran[k] = np.random.normal(f_val, std, N_models)
    # Transpose dict of arrays into list of dicts
    ran_models = [dict(zip(models_ran, t)) for t in zip(*models_ran.values())]

    return ran_models


def xxx(
    Nm, st_mass_mean, st_mass_var, Nm_binar, obs_phot, m_ini_idx,
    st_mass_mean_binar, st_mass_var_binar, prob_binar, binar_vals, alpha, isoch
):
    """
    Estimate the mean and variance for each star via recurrence.
    """
    # Masses, binary mask
    mass_primary = isoch[m_ini_idx]
    mass_secondary = isoch[-1]
    # shape: (N_stars, Ndim)
    photom = isoch[:m_ini_idx].T

    if alpha is not None:
        # Binaries have M2>0
        binar_idxs = isoch[-1] > 0.0
        binar_frac = binar_idxs.sum() / isoch.shape[-1]
    else:
        # No binaries were defined
        binar_idxs = np.full(isoch.shape[1], False)
        binar_frac = 0.
    binar_vals.append(binar_frac)

    # For non-binary systems
    photom_single = photom[~binar_idxs]
    if photom_single.any():
        Nm += 1
        obs_mass, lkl_p = photomMatch(
            obs_phot, photom_single, mass_primary[~binar_idxs]
        )
        # Estimate mean and variance
        st_mass_mean, st_mass_var = recurrentStats(
            Nm, st_mass_mean, st_mass_var, obs_mass)

        # For binary systems
        if binar_idxs.sum() > 0:
            photom_binar = photom[binar_idxs]
            # If there are no binary systems, skip
            if photom_binar.any():
                Nm_binar += 1
                obs_mass, lkl_b = photomMatch(
                    obs_phot, photom_binar, mass_secondary[binar_idxs]
                )
                st_mass_mean_binar, st_mass_var_binar = recurrentStats(
                    Nm, st_mass_mean_binar, st_mass_var_binar, obs_mass
                )

                # Bayesian probability
                new_prob_binar = 1.0 / (1.0 + (lkl_p / lkl_b))
                prob_binar = recurrentStats(Nm, prob_binar, None, new_prob_binar)

    return Nm, Nm_binar, st_mass_mean, st_mass_var, st_mass_mean_binar, \
           st_mass_var_binar, binar_vals, prob_binar


def photomMatch(obs_phot, photom, mass_ini):
    """
    For each observed star in 'obs_phot', find the closest synthetic star in
    the (synthetic) photometric space 'photom', and assign the mass of that
    synthetic star to the observed star
    """
    tree = KDTree(photom)
    dd, ii = tree.query(obs_phot, k=1)

    # Assign masses to each observed star
    obs_mass = mass_ini[ii]

    # Likelihood is defined as the inverse of the distance
    lkl = 1.0 / dd

    return obs_mass, lkl


def recurrentStats(Nm, mean, var, newValue):
    """
    Source: en.wikipedia.org/wiki/
            Algorithms_for_calculating_variance#Welford's_online_algorithm
    """
    count = Nm + 1
    delta = newValue - mean
    mean += delta / count
    if var is None:
        return mean
    var += delta * (newValue - mean)
    return mean, var


def get_masses(masses_dict, ra, dec, m_ini_idx, st_dist_mass, isoch, loga, dm):
    """ """
    # N_obs = cl_dict["N_obs_stars"]
    mass_ini = isoch[m_ini_idx]
    M_obs = mass_ini.sum()
    mass_min, mass_max = mass_ini.min(), mass_ini.max()

    # Select a random IMF sampling array
    Nmets, Nages = len(st_dist_mass), len(st_dist_mass[0])
    i = np.random.randint(Nmets)
    j = np.random.randint(Nages)
    mass_samples = st_dist_mass[i][j]

    # def sampled_inv_cdf(N):
    #     mr = np.random.rand(N)
    #     return inv_cdf(mr)
    # mass_samples = sampled_inv_cdf(500000)

    mass_tot = np.cumsum(mass_samples)

    gamma = 0.62
    t = 10**loga
    qev_t = stellar_evol_mass_loss(loga)
    term1 = (1-qev_t)**gamma
    t0 = minit_LGB05(ra, dec, dm)

    def func_optm(N_max, flag_mass=False):
        M_init_sample = mass_samples[:int(N_max)]
        M_init = mass_tot[int(N_max)]

        # Masks for the min-max mass range
        msk_max = M_init_sample > mass_max
        M_ev = M_init_sample[msk_max].sum()
        M_actual_dyn = M_init - M_ev
        msk_min = M_init_sample < mass_min
        # This is the percentage of mass below the photometric minimum mass limit
        M_ratio = M_init_sample[msk_min].sum() / M_actual_dyn

        term2 = (gamma/(M_init**gamma))*(t/t0)
        M_actual, M_actual_range = 0, 0
        if term1 > term2:
            M_actual = M_init * (term1 - term2)**(1/gamma)
            M_actual_range = M_actual - M_actual*M_ratio

        if flag_mass:
            return M_actual, M_init, abs(M_actual_range - M_obs)
        return abs(M_actual_range - M_obs)

    # Perform grid search
    optimal_param = grid_search_optimal_parameter(func_optm, 100, len(mass_samples))
    M_actual, M_init, mass_diff = func_optm(optimal_param, True)

    # print(round(loga, 2), int(M_actual), int(M_init), int(mass_diff))

    masses_dict['M_actual'].append(M_actual)
    masses_dict['M_init'].append(M_init)

    return masses_dict


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
    """ Fraction of the initial cluster mass (Mini) lost by stellar evolution"""
    a, b, c = 7, 0.26, -1.8
    q_ev = 10**((max(7.1, loga) - a)**b + c)
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
    c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
    lon, lat = c.galactic.l, c.galactic.b
    dist_pc = 10**(.2*(dm+5))
    c = SkyCoord(l=lon, b=lat, distance=dist_pc*u.pc, frame='galactic')
    c.representation_type = 'cylindrical'
    z = c.z.value
    # Estimate R_GC
    gc = c.transform_to(coord.Galactocentric(galcen_distance=8*u.kpc, z_sun=0*u.pc))
    r = np.sqrt(gc.x.value**2+gc.y.value**2)

    # Hernquist potential (bulge)
    # https://galaxiesbook.org/chapters/I-01.-Potential-Theory-and-Spherical-Mass-Distributions.html; Eq 3.58
    # This source above gives a density that does not match the two below formulas:
    # r_B/(r*(1+r/R_B)**3) <-- ?
    # https://docs.galpy.org/en/latest/reference/potentialhernquist.html
    # "note that amp is 2 x [total mass] for the chosen definition of the Two Power Spherical potential"
    # Phi_B_Laplacian = 2*M_B/r_B**3 * 1/((r/r_B)*(1+r/r_B)**3)
    # The above is equivalent to this one by
    # https://academic.oup.com/mnras/article/428/4/2805/992063; Eq 1
    Phi_B_Laplacian = 2*M_B*r_B/(r*(r+r_B)**3)

    # Miyamoto & Nagai potential (disk)
    # https://galaxiesbook.org/chapters/II-01.-Flattened-Mass-Distributions.html#Thickened-disk:-the-Miyamoto-Nagai-model
    # https://articles.adsabs.harvard.edu/pdf/1975PASJ...27..533M
    # https://www.astro.utu.fi/~cflynn/galdyn/lecture4.html
    numerator = M_D * b**2 * (
        a*r**2 + (a + 3*np.sqrt(z**2 + b**2)) * (a + np.sqrt(z**2 + b**2))**2)
    denominator = (b**2 + z**2)**(3/2) * (r**2 + (a + np.sqrt(b**2 + z**2))**2)**(5/2)
    Phi_D_laplacian = numerator / denominator

    # Sanderson potential (dark matter halo)
    A = -M_s / (np.log(2) - 0.5)
    Phi_H_Laplacian = -A/(r*(r+r_s)**2)

    # Ambient density
    rho_amb = (1/(4*np.pi)) * (Phi_B_Laplacian + Phi_D_laplacian + Phi_H_Laplacian)
    t0 = C_env0*(1-epsilon)*10**(-4*gamma)*rho_amb**(-.5)

    return t0
    # Minit = (M_actual**gamma + gamma*(t/t0))**(1/gamma)/(1-q_ev)
