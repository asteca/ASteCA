
import os
from os.path import join, exists
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.stats import gaussian_kde

from ..structure.king_profile import KingProf
# from ..best_fit.bf_common import varPars
from .. import update_progress
from . import synth_cluster
from . add_errors import getSigmas
from ..out import synth_gen_out
from ..out import make_D0_plot


def main(npd, clp, pd, rt=0.05, synth_CI_rand=False):
    """
    In place for #239
    """

    print("Generating synthetic clusters")

    kcp = clp['KP_conct_par']
    if not np.isnan(kcp):
        rc = rt / (10 ** kcp)
    else:
        rc = np.random.uniform(.5 * rt, rt)

    # if synth_CI_rand is True:
    #     CI = np.random.uniform(.1, .95)
    # else:
    #     CI = clp['cont_index']
    # CI = 0. if np.isnan(CI) else CI
    CI = .8

    # Handle cases where no field region is defined
    if len(clp['field_regions_c']) > 0:
        mags_fl = []
        for fl in clp['field_regions_c']:
            mags_fl += np.ravel(list(zip(*fl))[3]).tolist()
        mags_fl = np.array(mags_fl)
        # Generate mags KDE
        x_grid_fr = np.linspace(min(mags_fl), max(mags_fl), 1000)
        kde = gaussian_kde(mags_fl)
        kdepdf_fr = kde.evaluate(x_grid_fr)
    else:
        kdepdf_fr = None

    # TODO this should come from params_input.dat
    z_vals = (0.0151,)
    a_vals = (8.5,)
    e_vals = (0.5,)
    d_vals = (12.,)
    m_vals = (2000.,)
    b_vals = (0.3,)

    # varIdxs, ndim, ranges = varPars(fundam_params)
    # Take model values from the above list.
    varIdxs = (0, 1, 2, 3, 4, 5)

    models = []
    for z in z_vals:
        for a in a_vals:
            for e in e_vals:
                for d in d_vals:
                    for m in m_vals:
                        for b in b_vals:
                            models.append([z, a, e, d, m, b])

    # models = np.array([z_vals, a_vals, d_vals, e_vals, m_vals, b_vals]).T
    for i, model in enumerate(models):
        # model = np.array([np.random.uniform(*_) for _ in ranges])
        # model_var = model[varIdxs]

        # synth_clust, sigma, extra_pars, isoch_moved, mass_dist, isoch_binar,\
        #     isoch_compl = synth_cluster.main(
        synth_clust = synth_cluster.main(
            pd['fundam_params'], varIdxs, model, clp['completeness'],
            clp['err_lst'], clp['em_float'], clp['max_mag_syn'],
            pd['ext_coefs'], pd['binar_flag'], pd['mean_bin_mr'],
            pd['N_fc'], pd['m_ini_idx'], pd['st_dist_mass'],
            pd['theor_tracks'], pd['err_norm_rand'],
            pd['binar_probs'], pd['ext_unif_rand'], pd['R_V'])

        # Undo transposing performed in add_errors()
        synth_clust = synth_clust.T

        # Get uncertainties
        sigma = []
        for i, popt_mc in enumerate(clp['err_lst']):
            sigma.append(getSigmas(synth_clust[0], popt_mc))

        # sigma = np.array(sigma)
        # Larger errors
        sigma = np.array([sigma[0] * 10, sigma[0] * 5])

        # Generate positional data
        field_dens, cl_dists, x_cl, y_cl, x_fl, y_fl = xyCoords(
            synth_clust.shape[1], CI, rc, rt)

        if kdepdf_fr is not None:
            # Generate field stars' photometry
            synth_field, sigma_field = fldStarsPhot(
                kdepdf_fr, x_grid_fr, clp['max_mag_syn'], synth_clust, sigma,
                len(x_fl))

            # Clip at 'max_mag_syn'
            msk = synth_field[0] < clp['max_mag_syn']
            synth_field, sigma_field = synth_field[:, msk], sigma_field[:, msk]
            x_fl, y_fl = x_fl[msk], y_fl[msk]
        else:
            x_fl, y_fl, synth_field, sigma_field = [
                np.array([]) for _ in range(4)]

        data_file, plot_file = fileName(npd, model)

        # Output data to file.
        synth_gen_out.createFile(
            pd['filters'], pd['colors'], model, synth_clust, sigma,
            x_cl, y_cl, x_fl, y_fl, synth_field, sigma_field, CI, rc, rt,
            data_file)

        # if 'D0' in flag_make_plot:
        #     make_D0_plot.main(
        #         model, isoch_moved, mass_dist, isoch_binar, isoch_compl,
        #         synth_clust, extra_pars, sigma, synth_field, sigma_field, cx,
        #         cy, rc, rt, cl_dists, xmax, ymax, x_cl, y_cl, x_fl, y_fl, CI,
        #         max_mag_syn, flag_make_plot, coords, colors, filters,
        #         plot_file)
        #     print("<<Plots for D0 block created>>")
        # else:
        #     print("<<Skip D0 plot>>")

        update_progress.updt(len(models), i + 1)

    return


def xyCoords(N_clust, CI, rc, rt, cx=0., cy=0.):
    """
    """
    length = rt * 5

    # Estimate number of field stars, given CI, N_clust, and rt
    area_frame = length**2
    area_cl = np.pi * rt**2

    # Generate positions for field stars
    N_fl, field_dens = estimateNfield(N_clust, CI, area_frame, area_cl)
    x_fl = np.random.uniform(-length, length, N_fl)
    y_fl = np.random.uniform(-length, length, N_fl)

    # Sample King's profile with fixed rc, rt values.
    cl_dists = invTrnsfSmpl(rc, rt, N_clust)

    # Generate positions for cluster members, given heir KP distances to the
    # center.
    theta = np.random.uniform(0., 1., N_clust) * 2 * np.pi
    x_cl = cx + cl_dists * np.cos(theta)
    y_cl = cy + cl_dists * np.sin(theta)

    return field_dens, cl_dists, x_cl, y_cl, x_fl, y_fl


def invTrnsfSmpl(rc, rt, N_samp, N_interp=1000):
    """
    Sample King's profile using the inverse CDF method.
    """
    def rKP(r, rc, rt):
        return r * KingProf(r, rc, rt)

    r_0rt = np.linspace(0., rt, N_interp)
    # The CDF is defined as: $F(r)= \int_{r_low}^{r} PDF(r) dr$
    # Sample the CDF
    CDF_samples = []
    for r in r_0rt:
        CDF_samples.append(quad(rKP, 0., r, args=(rc, rt))[0])

    # Normalize CDF
    CDF_samples = np.array(CDF_samples) / CDF_samples[-1]

    # Inverse CDF
    inv_cdf = interp1d(CDF_samples, r_0rt)

    # Sample the inverse CDF
    samples = inv_cdf(np.random.rand(N_samp))

    return samples


def estimateNfield(N_membs, CI, tot_area, cl_area):
    """
    Estimate the total number of field stars that should be generated so
    that the CI is respected.
    """

    # Number of field stars in the cluster area
    if CI == 0.:
        N_field_in_clreg = 0.
    else:
        N_field_in_clreg = N_membs / (1. / CI - 1.)

    # Field stars density
    field_dens = N_field_in_clreg / cl_area

    # Total number of field stars in the entire frame
    N_field = int(field_dens * tot_area)

    return N_field, field_dens


def fldStarsPhot(
        kdepdf_fr, x_grid_fr, max_mag_syn, synth_clust, sigma, N_fl, std_f=5.):
    """

    std_f : number of standard deviations used to perturb the uncertainties
    assigned to the field stars.
    """
    ndim, N_cl = synth_clust.shape

    # *Ordered* synthetic cluster's magnitude
    mags_sort = np.argsort(synth_clust[0])
    synth_clust = synth_clust[:, mags_sort]
    sigma = sigma[:, mags_sort]
    mags_cl = synth_clust[0]

    # # This method does not use the synthetic cluster to generate the field
    # # photometry. The advantage is that it can be controlled with the
    # # 'mag_CI' parameter.
    # #
    # # Exponentiate magnitudes so that stars with larger values are assigned
    # # larger probabilities ans are thus more likely to be sampled below.
    # mag_CI = .5
    # probs = np.exp(mag_CI * synth_clust[0])
    # probs = probs / probs.sum()

    def stretch(x, xmin, xmax):
        return np.interp(x, (x.min(), x.max()), (xmin, xmax))

    def generate_rand_from_pdf(pdf, x_grid, N_fl):
        """
        Source: https://stackoverflow.com/a/29488780/1391441
        """
        cdf = np.cumsum(pdf)
        cdf = cdf / cdf[-1]
        values = np.random.rand(N_fl)
        value_bins = np.searchsorted(cdf, values)
        random_from_cdf = x_grid[value_bins]
        return random_from_cdf

    # Generate N_fl samples from the KDE of the 'mags'.
    mags_kde = generate_rand_from_pdf(kdepdf_fr, x_grid_fr, N_fl)

    # Stretch field regions magnitudes to the synthetic cluster's range.
    mags_kde = stretch(mags_kde, mags_cl.min(), mags_cl.max())
    mags_kde = np.sort(mags_kde)

    # For every mag value sampled from the KDE, find the closest synthetic
    # magnitude value.
    fl_ids = np.searchsorted(mags_cl, mags_kde)
    # Pull elements at the end of the array back one position.
    fl_ids[fl_ids >= len(sigma[0]) - 1] = len(sigma[0]) - 1

    # Generate uncertainties for field stars by randomly perturbing the
    # uncertainties from the (sampled) cluster members.
    sigma_field = sigma[:, fl_ids]
    sigma_field = sigma_field + np.random.normal(0, .1, N_fl) * sigma_field

    # Perturb the magnitudes
    mag_fl = mags_kde + np.random.normal(0., std_f, N_fl) * sigma_field[0]

    # Same for colors
    cols_fl = synth_clust[1:, fl_ids] +\
        np.random.normal(0., std_f, (ndim - 1, N_fl)) * sigma_field[1:]

    # Add an extra perturbation term to the magnitude and colors.
    mag_fl += np.random.normal(0., .5, N_fl)
    cols_fl += np.random.normal(0., .2, (ndim - 1, N_fl))

    # Combine magnitude and color(s) to generate the final field photometry.
    synth_field = np.concatenate(([mag_fl], cols_fl))

    return synth_field, sigma_field


def fileName(npd, model):
    """
    """

    # (z, a) subfolder name
    subfolder = str("{:.7f}".format(model[0])).replace("0.", "") +\
        "_" + str("{:.3f}".format(model[1])).replace(".", "")

    synth_gen_fold = join(npd['output_subdir'], 'synth', subfolder)
    if not exists(synth_gen_fold):
        os.makedirs(synth_gen_fold)

    clust_name = str("{:.3f}".format(model[2])).replace(".", "") +\
        "_" + str("{:.2f}".format(model[3])).replace(".", "") +\
        "_" + str("{:.0f}".format(model[4])).replace(".", "") +\
        "_" + str("{:.2f}".format(model[5])).replace(".", "")

    data_file = join(synth_gen_fold, clust_name + '.dat')
    plot_file = join(synth_gen_fold, clust_name)

    return data_file, plot_file
