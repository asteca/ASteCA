import logging
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.spatial import KDTree
from astropy.stats import knuth_bin_width, bayesian_blocks
from fast_histogram import histogram2d


def cluster_load(self):
    """ """
    print("Reading and processing cluster data")

    cl_ids = np.array(self.cluster_df[self.id_column])
    mag = np.array(self.cluster_df[self.mag_column])
    colors = [np.array(self.cluster_df[_]) for _ in self.col_column]
    e_mag = np.array(self.cluster_df[self.e_mag_column])
    e_colors = [np.array(self.cluster_df[_]) for _ in self.e_col_column]

    # Obtain bin edges for each dimension, defining a grid.
    bin_edges, ranges, Nbins = bin_edges_f(
        self.bin_method, self.N_mag, self.N_col, mag, colors
    )

    # Obtain histogram for observed cluster.
    hess_diag = []
    for i, col in enumerate(colors):
        # hess_diag.append(
        #     np.histogram2d(
        #         mag, col, bins=[bin_edges[0]] + [bin_edges[i + 1]])[0])

        hess_diag = histogram2d(
            mag,
            col,
            range=[[ranges[0][0], ranges[0][1]], [ranges[i + 1][0], ranges[i + 1][1]]],
            bins=[Nbins[0], Nbins[i + 1]],
        )

    # Flatten array
    cl_histo_f = []
    # probs_f = []
    for i, diag in enumerate(hess_diag):
        cl_histo_f += list(np.array(diag).ravel())
    cl_histo_f = np.array(cl_histo_f)

    # Down sample histogram
    hist_down_samp = np.array(cl_histo_f)
    msk = hist_down_samp > 5
    hist_down_samp[msk] = 5
    # hist_down_samp[~msk] = 0.

    # Index of bins where stars were observed
    cl_z_idx = cl_histo_f != 0

    # Remove all bins where n_i=0 (no observed stars)
    cl_histo_f_z = cl_histo_f[cl_z_idx]

    # Used by the synthetic cluster module
    max_mag_syn = max(mag)
    N_obs_stars = len(mag)
    N_colors = len(colors)
    err_lst = error_distribution(mag, e_mag, e_colors)

    msk = np.isnan(mag) | np.isnan(colors[0])
    mag0 = mag[~msk]
    colors0 = colors[0][~msk]

    cluster_dict = {
        "cl_ids": cl_ids,
        "mag": mag,
        "colors": colors,
        "bin_edges": bin_edges,
        "ranges": ranges,
        "Nbins": Nbins,
        "cl_z_idx": cl_z_idx,
        "cl_histo_f_z": cl_histo_f_z,
        "max_mag_syn": max_mag_syn,
        "N_obs_stars": N_obs_stars,
        "N_colors": N_colors,
        "err_lst": err_lst,
        #
        "hist_down_samp": hist_down_samp,
        "mag0": mag0,
        "colors0": colors0,
    }

    return cluster_dict


def bin_edges_f(bin_method, N_mag, N_col, mag, colors):
    """ """

    bin_edges = []

    if bin_method == "knuth":
        bin_edges.append(
            knuth_bin_width(mag[~np.isnan(mag)], return_bins=True, quiet=True)[1]
        )
        for col in colors:
            bin_edges.append(
                knuth_bin_width(col[~np.isnan(col)], return_bins=True, quiet=True)[1]
            )

    elif bin_method == "fixed":
        # Magnitude
        mag_min, mag_max = np.nanmin(mag), np.nanmax(mag)
        bin_edges.append(np.linspace(mag_min, mag_max, N_mag))
        # Colors
        for col in colors:
            col_min, col_max = np.nanmin(col), np.nanmax(col)
            bin_edges.append(np.linspace(col_min, col_max, N_col))

    elif bin_method == "bayes_blocks":
        bin_edges.append(bayesian_blocks(mag[~np.isnan(mag)]))
        for col in colors:
            bin_edges.append(bayesian_blocks(col[~np.isnan(col)]))

    # Extract ranges and number of bins, used by histogram2d
    ranges, Nbins = [], []
    for be in bin_edges:
        ranges.append([be[0], be[-1]])
        Nbins.append(len(be))

    return bin_edges, ranges, Nbins


def error_distribution(mag, e_mag, e_colors):
    """
    Fit an exponential function to the errors in each photometric dimension,
    using the main magnitude as the x coordinate.
    This data is used to display the error bars, and more importantly, to
    generate the synthetic clusters in the best match module.
    """

    # not_nan_msk = ~(np.isnan(mag) | np.isnan(e_mag) | np.isnan(e_colors[0]))
    # mag, e_mag, e_color = mag[not_nan_msk], e_mag[not_nan_msk], e_colors[0][not_nan_msk]

    # import matplotlib.pyplot as plt
    # from astropy.stats import sigma_clip

    # mag_msk = ~(sigma_clip(mag).mask)
    # e_mag_msk = ~(sigma_clip(e_mag).mask)
    # e_col_msk = ~(sigma_clip(e_color).mask)

    # msk = mag_msk & e_mag_msk
    # mag_clip = mag[msk]
    # e_mag_clip = e_mag[msk]
    # popt_mc, _ = curve_fit(exp_3p, mag_clip, e_mag_clip)
    # y_exp = exp_3p(mag_clip, *popt_mc)
    # plt.scatter(mag, e_mag, alpha=.5)
    # plt.scatter(mag_clip, y_exp)
    # plt.show()

    # msk = mag_msk & e_col_msk
    # mag_clip = mag[msk]
    # e_color_clip = e_color[msk]
    # popt_mc, _ = curve_fit(exp_3p, mag_clip, e_color_clip)
    # y_exp = exp_3p(mag_clip, *popt_mc)
    # plt.scatter(mag, e_color, alpha=.5)
    # plt.scatter(mag_clip, y_exp)
    # plt.show()
    # breakpoint()

    # Mask of not nan values across arrays
    nan_msk = np.isnan(mag) | np.isnan(e_mag)
    for e_col in e_colors:
        nan_msk = nan_msk | np.isnan(e_col)
    not_nan_msk = ~nan_msk
    # Remove nan values
    mag, e_mag = mag[not_nan_msk], e_mag[not_nan_msk]
    e_col_non_nan = []
    for e_col in e_colors:
        e_col_non_nan.append(e_col[not_nan_msk])
    e_colors = e_col_non_nan

    # Left end of magnitude range
    be_m = max(min(mag) + 1.0, np.percentile(mag, 0.5))
    # Width of the intervals in magnitude.
    interv_mag = 0.5
    # Number of intervals.
    delta_mag = mag.max() - be_m
    n_interv = int(round(delta_mag / interv_mag))
    #
    steps_x = np.linspace(be_m - 0.5 * interv_mag, mag.max(), n_interv - 1)

    # Median values for each error array in each magnitude range
    mag_y = []
    for i, e_mc in enumerate([e_mag] + [list(_) for _ in e_colors]):
        x1 = steps_x[0]
        e_mc_medians = []
        for x2 in steps_x[1:]:
            msk = (mag >= x1) & (mag < x2)
            strs_in_range = np.array(e_mc)[msk]
            if len(strs_in_range) > 1:
                e_mc_medians.append(np.median(strs_in_range))
            else:
                # If no stars in interval, use fixed value
                e_mc_medians.append(0.0001)
            x1 = x2
        mag_y.append(e_mc_medians)

    # Make sure that median error values increase with increasing magnitude. This
    # ensures that the 3P exponential fit does not fail
    mag_y_new = []
    for e_arr in mag_y:
        e_arr_new, v_old = [], np.inf
        for i in range(-1, -len(e_arr) - 1, -1):
            if e_arr[i] > v_old:
                e_arr_new.append(v_old)
            else:
                e_arr_new.append(e_arr[i])
            v_old = e_arr[i]
        e_arr_new.reverse()
        mag_y_new.append(e_arr_new)
    mag_y = mag_y_new

    # Mid points in magnitude range
    mag_x = 0.5 * (steps_x[:-1] + steps_x[1:])

    # Fit 3-parameter exponential
    err_lst = []
    for y in mag_y:
        popt_mc, _ = curve_fit(exp_3p, mag_x, y)
        err_lst.append(popt_mc)

    # import matplotlib.pyplot as plt
    # plt.scatter(mag, e_mag, alpha=.5)
    # y_exp = exp_3p(mag_x, *err_lst[0])
    # plt.plot(mag_x, y_exp, c='orange')
    # plt.show()

    # y_exp = exp_3p(mag_x, *err_lst[1])
    # plt.scatter(mag, e_colors[0], alpha=.5)
    # plt.plot(mag_x, y_exp, c='orange')
    # plt.show()
    # breakpoint()

    return err_lst


def exp_3p(x, a, b, c):
    """
    Three-parameters exponential function.

    This function is tied to the 'synth_cluster.add_errors' function.
    """
    return a * np.exp(b * x) + c


def get_masses_binar(my_cluster, synthcl, model_fit, model_std):
    """
    Assign individual masses to the observed cluster's stars, along with binary
    probabilities (if binarity was estimated).

    Estimate the statistics for the mass and binary fractions (not fitted)
    """
    print("Estimating total mass, binary probabilities, and per star masses")
    cl_dict = my_cluster.cluster_dict
    ra, dec = my_cluster.ra, my_cluster.dec
    model_fixed = my_cluster.model_fixed
    m_ini_idx = cl_dict["N_colors"] + 1
    # Extract photometry used in the best fit process
    mags_cols_cl = [cl_dict["mag"]] + [_ for _ in cl_dict["colors"]]
    obs_phot = np.array(mags_cols_cl)

    # Generate random models from the selected solution
    models = ranModels(model_fit, model_std)

    # Identify nans in mag and colors and re-write them as -10
    nan_msk = np.full(obs_phot.shape[1], False)
    for col in obs_phot:
        nan_msk = nan_msk | np.isnan(col)
    obs_phot[:, nan_msk] = -10.0
    obs_phot = obs_phot.T
    N_obs = obs_phot.shape[0]

    # Initiate empty arrays for mean and variance
    st_mass_mean, st_mass_var = np.zeros(N_obs), np.zeros(N_obs)
    st_mass_mean_binar, st_mass_var_binar = np.zeros(N_obs), np.zeros(N_obs)
    prob_binar = np.zeros(N_obs)
    binar_vals = []
    Nm, Nm_binar = 0, 0

    from .synth_cluster import invTrnsfSmpl
    inv_cdf = invTrnsfSmpl(synthcl.IMF_name)

    masses_dict = {'M_actual': [], 'M_init': [], 'M_init_LGB05': []}
    for _, model in enumerate(models):
        # Generate synthetic cluster from the 'model'.
        isoch = synthcl.generate(model, my_cluster, True)
        if not isoch.any():
            continue

        Nm, Nm_binar, st_mass_mean, st_mass_var, st_mass_mean_binar, st_mass_var_binar, binar_vals, prob_binar = xxx(
            Nm, st_mass_mean, st_mass_var, Nm_binar, obs_phot, m_ini_idx, st_mass_mean_binar,
            st_mass_var_binar, prob_binar, binar_vals, isoch)

        # Extract dm and loga
        model_comb = model_fixed | model
        loga, dm = model_comb['loga'], model_comb['dm']
        #
        masses_dict = get_masses(
            masses_dict, cl_dict, ra, dec, m_ini_idx, inv_cdf, isoch, loga, dm)

    # Store standard deviations instead of variances
    st_mass_std = np.sqrt(st_mass_var / Nm)
    st_mass_std_binar = np.sqrt(st_mass_var_binar / max(1, Nm_binar))

    # if synthcl.max_mass < np.median(M_init_arr) + np.std(M_init_arr):
    #     logging.warning(
    #         "The total mass is too close to the 'synth_clusters.max_mass' "
    #         + "parameter. Consider increasing this value."
    #     )

    # Use nan values for stars with initial nan photometric values
    for arr in (
        st_mass_mean,
        st_mass_std,
        st_mass_mean_binar,
        st_mass_std_binar,
        prob_binar,
    ):
        arr[nan_msk] = np.nan

    cl_masses_bfr = pd.DataFrame(
        {
            "ID": cl_dict["cl_ids"],
            "M1": st_mass_mean,
            "M1_std": st_mass_std,
            "M2": st_mass_mean_binar,
            "M2_std": st_mass_std_binar,
            "binar_prob": prob_binar,
        }
    )

    return masses_dict, np.array(binar_vals), cl_masses_bfr


def ranModels(model_fit, model_std, N_models=1000):
    """
    Generate the requested models via sampling a Gaussian centered on the
    selected solution, with standard deviation given by the attached
    uncertainty.

    N_models: number of models to generate (HARDCODED)
    """
    models_ran = {}
    for k, f_val in model_fit.items():
        std = model_std[k]
        models_ran[k] = np.random.normal(f_val, std, N_models)
    # Transpose dict of arrays into list of dicts
    ran_models = [dict(zip(models_ran, t)) for t in zip(*models_ran.values())]

    return ran_models


def xxx(Nm, st_mass_mean, st_mass_var, Nm_binar, obs_phot, m_ini_idx, st_mass_mean_binar, st_mass_var_binar, prob_binar, binar_vals, isoch):
    """
    Estimate the mean and variance for each star via recurrence.
    """
    # Masses, binary mask
    mass_primary = isoch[m_ini_idx]
    mass_secondary = isoch[-1]
    # shape: (N_stars, Ndim)
    photom = isoch[:m_ini_idx].T

    # Binaries have M2>0
    binar_idxs = isoch[-1] > 0.0
    binar_frac = binar_idxs.sum() / isoch.shape[-1]
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

    return Nm, Nm_binar, st_mass_mean, st_mass_var, st_mass_mean_binar, st_mass_var_binar, binar_vals, prob_binar


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


def get_masses(masses_dict, cl_dict, ra, dec, m_ini_idx, inv_cdf, isoch, loga, dm):
    """ """
    N_obs = cl_dict["N_obs_stars"]
    mass_ini = isoch[m_ini_idx]
    M_obs = mass_ini.sum()

    mass_min, mass_max = mass_ini.min(), mass_ini.max()

    def sampled_inv_cdf(N):
        mr = np.random.rand(N)
        return inv_cdf(mr)

    # Sample in chunks until the maximum defined mass is reached. This strategy is
    # a good balance between speed and accuracy
    N_chunk = max(10, int(N_obs * .25))

    qev_t = stellar_evol_mass_loss(loga)
    mass_all = []
    mass_inrange = []

    M_actual, M_phot, M_ev = 0, 0, 0
    N_obs_synth = 0
    while N_obs_synth < N_obs:
        # Sample the IMF
        mass_samples = sampled_inv_cdf(N_chunk)

        mass_all += list(mass_samples)

        # Masks for the min-max mass range
        msk_min = mass_samples > mass_min
        msk_max = mass_samples < mass_max

        # Mass below the minimum mass
        M_phot += mass_samples[~msk_min].sum()
        # Mass below the max mass range
        M_actual += mass_samples[msk_max].sum()
        # Mass above the maximum mass
        M_ev += mass_samples[~msk_max].sum()

        # Number of stars within the min-max mass range
        N_obs_synth += (msk_min & msk_max).sum()

        mass_inrange += list(mass_samples[(msk_min & msk_max)])

    # M_actual ~ M_obs + M_phot

    M_init_ev = M_ev / qev_t
    # M_dyn = M_init - M_actual - M_ev
    print(N_obs_synth - N_obs, int(M_actual), round(qev_t, 3), int(M_ev), int(M_init_ev), int(sum(mass_all)))

    # if M_actual>9000:
    #     breakpoint()

    # import matplotlib.pyplot as plt
    # idx = np.argsort(mass_all)
    # plt.plot(np.array(mass_all)[idx], np.cumsum(np.array(mass_all)[idx]))
    # plt.xscale('log')
    # plt.show()

    # t = 10**loga
    # gamma = 0.62
    # C_env0 = 810e6
    # epsilon = 0.08
    # M_init = sum(mass_all)
    # t0 = gamma*t/((M_init*(1-qev_t))**gamma - M_actual**gamma)
    # rho_amb = ((C_env0*(1-epsilon)*10**(-4*gamma))/t0)**2
    # print("{:.3f}, {:.3f}".format(t0/1e6, rho_amb))

    M_init_LGB05 = minit_LGB05(ra, dec, loga, dm, M_actual, qev_t, M_init_ev)

    # import matplotlib.pyplot as plt
    # plt.hist(mass_all, 50)
    # plt.xscale('log')
    # plt.yscale('log')
    # breakpoint()

    masses_dict['M_actual'].append(M_actual)
    masses_dict['M_init'].append(M_init_ev)
    masses_dict['M_init_LGB05'].append(M_init_LGB05)

    return masses_dict


def stellar_evol_mass_loss(loga):
    """ Fraction of the initial cluster mass (Mini) lost by stellar evolution"""
    a, b, c = 7, 0.26, -1.8
    q_ev = 10**((max(7.1, loga) - a)**b + c)
    return q_ev


def minit_LGB05(ra, dec, loga, dm, M_actual, q_ev, M_init_ev, epsilon=0.08):
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
    t = 10**loga

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

    Minit = (M_actual**gamma + gamma*(t/t0))**(1/gamma)/(1-q_ev)

    # print("{:.3f}, {:.3f}".format(t0/1e6, rho_amb))

    return Minit


def mplot(my_cluster, synthcl, model_fit):
    """ """
    # m_ini_idx = cluster_dict['N_colors'] + 1

    mag_col = my_cluster.mag_column
    col_col = my_cluster.col_column
    cl_dict = my_cluster.cluster_dict

    # Generate synthetic cluster.
    synth_clust = synthcl.generate(model_fit, my_cluster, True)
    # plt.title(f"Lkl dist ={round(final_dist, 3)}")
    # y_edges, x_edges = cluster_dict["bin_edges"]
    # for xe in x_edges:
    #     plt.axvline(xe, c="grey", ls=":")
    # for ye in y_edges:
    #     plt.axhline(ye, c="grey", ls=":")

    # plt.figure(figsize=(8, 8), dpi=300)

    plt.scatter(
        cl_dict["colors"][0],
        cl_dict["mag"],
        c="green",
        alpha=0.5,
        label=f"Observed, N={cl_dict['N_obs_stars']}",
    )

    binar_idx = synth_clust[-1] != 0.0
    x_synth, y_synth = synth_clust[1], synth_clust[0]
    # Single synthetic systems
    plt.scatter(
        x_synth[~binar_idx],
        y_synth[~binar_idx],
        marker="^",
        c="#519ddb",
        alpha=0.5,
        label=f"Single, N={len(x_synth[~binar_idx])}",
    )
    # Binary synthetic systems
    plt.scatter(
        x_synth[binar_idx],
        y_synth[binar_idx],
        marker="v",
        c="#F34C4C",
        alpha=0.5,
        label=f"Binary, N={len(x_synth[binar_idx])}",
    )

    plt.xlabel(col_col)
    plt.ylabel(mag_col)

    plt.gca().invert_yaxis()
    plt.legend()
    plt.show()
    # plt.savefig(f"{fname}.png")
