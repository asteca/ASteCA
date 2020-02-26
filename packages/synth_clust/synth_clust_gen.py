
import os
from os.path import join, exists
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.stats import gaussian_kde

from ..structure.king_profile import KingProf
# from ..best_fit.bf_common import varPars
from .. import update_progress
from . synth_cluster import properModel
from . import zaWAverage
from . import move_isochrone
from . import cut_max_mag
from . import mass_distribution
from . import mass_interp
from . import binarity
from . import completeness_rm
from . import add_errors

from ..out import synth_gen_out
from ..out import make_D0_plot


def main(
    npd, clp, max_mag_syn, obs_clust, ext_coefs, st_dist_mass, N_fc, err_pars,
    fundam_params, theor_tracks, R_V, m_ini_idx, binar_flag, filters, colors,
    plot_frmt, flag_make_plot, coords, plot_dpi, xmax=2000, ymax=2000,
        synth_CI=False, rt=250., **kwargs):
    """
    In place for #239
    """

    print("Generating synthetic clusters")

    kcp = clp['KP_conct_par']
    if not np.isnan(kcp):
        rc = rt / (10 ** kcp)
    else:
        rc = np.random.uniform(10., rt - 10.)

    if synth_CI is True:
        CI = np.random.uniform(.1, .95)
    else:
        CI = clp['cont_index']

    # Position cluster in the center of the frame
    cx, cy = xmax * .5, ymax * .5

    # TODO handle cases where no field region is defined
    mags_fl = []
    for fl in clp['field_regions_c']:
        mags_fl += np.ravel(list(zip(*fl))[3]).tolist()
    mags_fl = np.array(mags_fl)
    # Generate mags KDE
    x_grid_fr = np.linspace(min(mags_fl), max(mags_fl), 1000)
    kde = gaussian_kde(mags_fl)
    kdepdf_fr = kde.evaluate(x_grid_fr)

    # TODO this should come from params_input.dat
    z_vals = (0.005, 0.0152, 0.03)
    a_vals = (7., 8., 9.)
    e_vals = (1.,)
    d_vals = (14.,)
    m_vals = (500., 2000.)
    b_vals = (0.3,)
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

        isoch_moved, mass_dist, isoch_binar, isoch_compl,\
            (synth_clust, sigma, extra_pars) = synth_cluster(
                fundam_params, varIdxs, model, theor_tracks,
                clp['completeness'], max_mag_syn, st_dist_mass, R_V, ext_coefs,
                N_fc, err_pars, m_ini_idx, binar_flag)

        # Generate positional data
        field_dens, cl_dists, x_cl, y_cl, x_fl, y_fl = xyCoords(
            synth_clust.shape[1], CI, rc, rt, xmax, ymax, cx, cy)

        # Generate field stars' photometry
        synth_field, sigma_field = fldStarsPhot(
            kdepdf_fr, x_grid_fr, max_mag_syn, synth_clust, sigma, len(x_fl))

        # Clip at 'max_mag_syn'
        msk = synth_field[0] < max_mag_syn
        synth_field, sigma_field = synth_field[:, msk], sigma_field[:, msk]
        x_fl, y_fl = x_fl[msk], y_fl[msk]

        data_file, plot_file = fileName(npd, model, plot_frmt)

        # Output data to file.
        synth_gen_out.createFile(
            filters, colors, extra_pars, model, synth_clust, sigma, x_cl,
            y_cl, x_fl, y_fl, synth_field, sigma_field, CI, rc, rt, data_file)

        make_D0_plot.main(
            model, isoch_moved, mass_dist, isoch_binar, isoch_compl,
            synth_clust, extra_pars, sigma, synth_field, sigma_field, cx, cy,
            rc, rt, cl_dists, xmax, ymax, x_cl, y_cl, x_fl, y_fl, CI,
            max_mag_syn, flag_make_plot, coords, colors, filters, plot_dpi,
            plot_file)

        update_progress.updt(len(models), i + 1)

    return


def synth_cluster(
    fundam_params, varIdxs, model, theor_tracks, completeness, max_mag_syn,
        st_dist_mass, R_V, ext_coefs, N_fc, err_pars, m_ini_idx, binar_flag):
    """
    Takes an isochrone and returns a synthetic cluster created according to
    a certain mass distribution.
    """

    # Return proper values for fixed parameters and parameters required
    # for the (z, log(age)) isochrone averaging.
    model_proper, z_model, a_model, ml, mh, al, ah = properModel(
        fundam_params, model, varIdxs)

    # Generate a weighted average isochrone from the (z, log(age)) values in
    # the 'model'.
    isochrone = zaWAverage.main(
        theor_tracks, fundam_params, z_model, a_model, ml, mh, al, ah)

    # Extract parameters
    e, d, M_total, bin_frac = model_proper

    # Move theoretical isochrone using the values 'e' and 'd'.
    isoch_moved = move_isochrone.main(
        isochrone, e, d, R_V, ext_coefs, N_fc, binar_flag, m_ini_idx)

    # Get isochrone minus those stars beyond the magnitude cut.
    isoch_cut = cut_max_mag.main(isoch_moved, max_mag_syn)

    # # Empty list to pass if at some point no stars are left.
    # synth_clust = np.array([])
    # Check for an empty array.
    if isoch_cut.any():
        # Mass distribution to produce a synthetic cluster based on
        # a given IMF and total mass.
        mass_dist = mass_distribution.main(st_dist_mass, M_total)

        # Interpolate masses in mass_dist into the isochrone rejecting those
        # masses that fall outside of the isochrone's mass range.
        # This destroys the order by magnitude.
        isoch_mass = mass_interp.main(isoch_cut, mass_dist, m_ini_idx)

        if isoch_mass.any():
            # Assignment of binarity.
            isoch_binar = binarity.main(isoch_mass, bin_frac, m_ini_idx, N_fc)

            # Completeness limit removal of stars.
            isoch_compl = completeness_rm.main(isoch_binar, completeness)

            if isoch_compl.any():
                # Get errors according to errors distribution.
                synth_clust = add_errors.main(
                    isoch_compl, err_pars, binar_flag, m_ini_idx)

    # if not synth_clust[0].any():
    #     raise ValueError("Synthetic cluster is empty: {}".format(model))

    return isoch_moved, mass_dist, isoch_binar, isoch_compl, synth_clust


def xyCoords(N_clust, CI, rc, rt, xmax, ymax, cx, cy):
    """
    """

    # Estimate number fo field stars, given CI, N_clust, and rt
    area_frame = xmax * ymax
    area_cl = np.pi * rt**2

    # Generate positions for field stars
    N_fl, field_dens = estimateNfield(N_clust, CI, area_frame, area_cl)
    x_fl = np.random.uniform(0., xmax, N_fl)
    y_fl = np.random.uniform(0., ymax, N_fl)

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


def fileName(npd, model, plot_frmt):
    """
    """

    # (z, a) subfolder name
    subfolder = str("{:.4f}".format(model[0])).replace("0.", "") +\
        "_" + str("{:.3f}".format(model[1])).replace(".", "")

    synth_gen_fold = join(npd['output_subdir'], 'synth', subfolder)
    if not exists(synth_gen_fold):
        os.makedirs(synth_gen_fold)

    clust_name = str("{:.3f}".format(model[2])).replace(".", "") +\
        "_" + str("{:.2f}".format(model[3])).replace(".", "") +\
        "_" + str("{:.0f}".format(model[4])).replace(".", "") +\
        "_" + str("{:.2f}".format(model[5])).replace(".", "")

    data_file = join(synth_gen_fold, clust_name + '.dat')
    plot_file = join(synth_gen_fold, clust_name + '.' + plot_frmt)

    return data_file, plot_file
