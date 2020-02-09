
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.stats import gaussian_kde

from ..structure.king_profile import KingProf
from ..best_fit.bf_common import varPars
from . synth_cluster import properModel
from . import zaWAverage
from . import move_isochrone
from . import cut_max_mag
from . import mass_distribution
from . import mass_interp
from . import binarity
from . import completeness_rm
from . import add_errors


def main(
    clp, max_mag_syn, obs_clust, ext_coefs, st_dist_mass, N_fc, err_pars,
    fundam_params, theor_tracks, R_V, m_ini_idx, binar_flag,
        xmax=2000, ymax=2000, synth_CI=False, rt=250., **kwargs):
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

    varIdxs, ndim, ranges = varPars(fundam_params)

    model = np.array([np.random.uniform(*_) for _ in ranges])
    model_var = model[varIdxs]

    isoch_moved, isoch_binar, isoch_compl, (synth_clust, sigma, extra_pars) =\
        synth_cluster(
        fundam_params, varIdxs, model_var, theor_tracks, clp['completeness'],
        max_mag_syn, st_dist_mass, R_V, ext_coefs, N_fc, err_pars, m_ini_idx,
        binar_flag)

    if not synth_clust.any():
        raise ValueError("Synthetic cluster is empty: {}".format(model))

    # Generate positional data
    field_dens, cl_dists, x_cl, y_cl, x_fl, y_fl = xyCoords(
        synth_clust.shape[1], CI, rc, rt, xmax, ymax, cx, cy)

    # Generate field stars' photometry
    synth_field, sigma_field = fldStarsPhot(
        max_mag_syn, synth_clust, sigma, len(x_fl))

    # Clip at 'max_mag_syn'
    msk = synth_field[0] < max_mag_syn
    synth_field, sigma_field = synth_field[:, msk], sigma_field[:, msk]
    x_fl, y_fl = x_fl[msk], y_fl[msk]

    # TODO this needs to be made cleaner
    clp['synth_gen_pars'] = (
        extra_pars, model, isoch_moved, isoch_binar, isoch_compl,
        synth_clust, sigma, cl_dists, x_cl, y_cl, x_fl, y_fl, synth_field,
        sigma_field, field_dens, CI, rc, rt, xmax, ymax, cx, cy)
    return clp


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

    # Empty list to pass if at some point no stars are left.
    synth_clust = np.array([])
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

    return isoch_moved, isoch_binar, isoch_compl, synth_clust


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


def fldStarsPhot(max_mag_syn, synth_clust, sigma, N_fl, std_f=5.):
    """

    std_f : number of standard deviations used to perturb the uncertainties
    assigned to the field stars.
    """
    ndim, N_cl = synth_clust.shape
    mags = np.sort(synth_clust[0])

    # # This method does not use the synthetic luster to generate the field
    # # photometry. The advantage is that it can be controlled with the
    # # 'mac_CI' parameter.
    # #
    # # Exponentiate magnitudes so that stars with larger values are assigned
    # # larger probabilities ans are thus more likely to be sampled below.
    # mag_CI = .5
    # probs = np.exp(mag_CI * synth_clust[0])
    # probs = probs / probs.sum()

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

    # Sample LF's KDE
    x_grid = np.linspace(min(mags), max(mags), 1000)
    kde = gaussian_kde(mags)
    kdepdf = kde.evaluate(x_grid)
    sampled_kde = generate_rand_from_pdf(kdepdf, x_grid, N_fl)

    fl_ids = np.searchsorted(mags, sampled_kde)

    # import matplotlib.pyplot as plt
    # plt.subplot(121)
    # plt.hist(mags, 50, normed=True, alpha=0.5, label='hist')
    # plt.plot(x_grid, kdepdf, color='r', alpha=0.5, lw=3, label='kde')
    # plt.legend()
    # plt.subplot(122)
    # plt.hist(sampled_kde, 50, alpha=0.5, label='from kde')
    # plt.legend()
    # # plt.show()

    # Generate uncertainties for field stars by randomly perturbing the
    # uncertainties from the (sampled) cluster members.
    sigma_field = sigma[:, fl_ids]
    sigma_field = sigma_field + np.random.normal(0, .1, N_fl) * sigma_field

    # Perturb the magnitudes
    mag_fl = synth_clust[0, fl_ids] +\
        np.random.normal(0., std_f, N_fl) * sigma_field[0]
    # Add an extra perturbation term.
    mag_fl = mag_fl + np.random.normal(0., .5, N_fl)

    # Same for colors
    cols_fl = synth_clust[1:, fl_ids] +\
        np.random.normal(0., std_f, (ndim - 1, N_fl)) * sigma_field[1:]
    cols_fl = cols_fl + np.random.normal(0., .2, (ndim - 1, N_fl))

    # Combine magnitude and color(s) to generate the final field photometry.
    synth_field = np.concatenate(([mag_fl], cols_fl))

    return synth_field, sigma_field
