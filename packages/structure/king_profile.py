
import numpy as np
from astropy.stats import circmean
from astropy import units as u
from scipy import spatial
from scipy import stats
import warnings
from .. import update_progress
from ..out import prep_plots
from ..best_fit.bf_common import modeKDE


def main(
    clp, cld_i, kp_flag, kp_nchains, kp_nruns, kp_nburn, rt_max_f, coords,
        **kwargs):
    """
    Bayesian inference over an array of stars' coordinates using a King
    profile model.

    The previously obtained 'clust_rad' value is used to determine the
    range where rt should be fitted.

    The field density value is fixed to the previously estimated value.
    The central density is inferred from the (rc, rt) values and the previously
    estimated number of members ('n_memb_i', which is *also* estimated
    using the field density), OR estimated per step (also using 'fd').

    This means that the value given to the field density has a *very large*
    influence on the final (rc, rt) values.

    I tried setting the field density as as free parameter but the sampler
    tends to maximize its value as much as possible. This is not physically
    reasonable.

    """
    if kp_flag:
        print("Estimating King profile")
        KP_plot, KP_Bys_rc, KP_Bys_rt, KP_Bys_ecc, KP_Bys_theta =\
            fit_King_prof(
                kp_nchains, kp_nruns, kp_nburn, cld_i['x'],
                cld_i['y'], clp['kde_cent'], clp['field_dens'],
                clp['clust_rad'], clp['n_memb_i'], rt_max_f)

        # Obtain number of members and concentration parameter.
        KP_memb_num, KP_conct_par = num_memb_conc_param(
            KP_plot['KP_cent_dens'], KP_Bys_rc[1], KP_Bys_rt[1])

        coord = prep_plots.coord_syst(coords)[0]
        # Set precision of printed values.
        text2 = '{:.1f}, {:.1f}' if coord == 'px' else '{:g}, {:g}'
        text = 'Core & tidal radii obtained: ' + text2 + ' {}'
        print(text.format(KP_Bys_rc[1], KP_Bys_rt[1], coord))

    else:
        print("Skipping King profile fit")
        KP_plot, KP_memb_num, KP_conct_par = {}, np.nan, np.nan
        KP_Bys_rc, KP_Bys_rt, KP_Bys_ecc, KP_Bys_theta =\
            [np.array([np.nan] * 5) for _ in range(4)]

    clp['KP_plot'], clp['KP_Bys_rc'], clp['KP_Bys_rt'],\
        clp['KP_Bys_ecc'], clp['KP_Bys_theta'], clp['KP_memb_num'],\
        clp['KP_conct_par'] = KP_plot, KP_Bys_rc, KP_Bys_rt, KP_Bys_ecc,\
        KP_Bys_theta, KP_memb_num, KP_conct_par
    return clp


def fit_King_prof(
    nchains, nruns, nburn, x, y, cl_cent, field_dens, cl_rad, n_memb_i,
        rt_max_f, N_integ=1000, N_conv=500, tau_stable=0.01):
    """
    N_integ : number of points used in the integration performed in the
    central density function

    If the chain is longer than N_conv times the estimated autocorrelation
    time and if this estimate changed by less than (tau_stable * 100)%

    ** TODO **

    * Generalize to exponential profile? See:
    https://www.aanda.org/articles/aa/abs/2016/02/aa27070-15/aa27070-15.html
    Eq (2),
    https://ned.ipac.caltech.edu/level5/Sept11/Graham/Graham2.html,
    Eq (1)

    * Make the emcee 'move' a general parameter

    * Align the theta with the y axis (North) instead of the x axis (as done
    in Martin et al. 2008; 'A Comprehensive Maximum Likelihood Analysis of
    the Structural Properties of Faint Milky Way Satellites')?

    * Fit the field density?
    ndim = 3 (rc, rt, fd)
    I don't think it can be done...

    """
    from emcee import ensemble
    from emcee import moves

    # HARDCODED ##########################
    # Move used by emcee
    mv = [
        (moves.DESnookerMove(), 0.1),
        (moves.DEMove(), 0.9 * 0.9),
        (moves.DEMove(gamma0=1.0), 0.9 * 0.1),
    ]
    # mv = moves.StretchMove()
    # mv = moves.KDEMove()

    # Steps to store Bayes params.
    KP_steps = int(nruns * .01)

    # Estimated number of members (previously obtained)
    N_memb = n_memb_i
    # Estimate N_memb with each sampler step
    # N_memb = None

    # Select the number of parameters to fit:
    # ndim = 2 fits (rc, rt)
    # ndim = 4 fits (rc, rt, ecc, theta)
    ndim = 2
    # HARDCODED ##########################

    # The tidal radius can not be larger than 'rt_max' times the estimated
    # "optimal" cluster radius. Used as a prior.
    rt_max = rt_max_f * cl_rad
    # Tidal radius array. Used for integrating
    rt_rang = np.linspace(0., rt_max, int(rt_max_f * N_integ))

    # Identify stars inside the cut-out given by the 'rt_max' value. Only these
    # stars will be processed below.
    xy = np.array((x, y)).T
    xy_cent_dist = spatial.distance.cdist([cl_cent], xy)[0]
    msk = xy_cent_dist <= rt_max
    xy_in = xy[msk].T
    r_in = xy_cent_dist[msk]

    # Initial positions for the sampler.
    # Dimensions: rc, rt
    rc_pos0 = np.random.uniform(.05 * rt_max, rt_max, nchains)
    rt_pos0 = np.random.uniform(rc_pos0, rt_max, nchains)
    if ndim == 2:
        pos0 = np.array([rc_pos0, rt_pos0]).T
    elif ndim == 4:
        # Dimensions: rc, rt, ecc, theta
        ecc = np.random.uniform(.0, 1., nchains)
        theta = np.random.uniform(-np.pi / 2., np.pi / 2., nchains)
        pos0 = np.array([rc_pos0, rt_pos0, ecc, theta]).T

    # Define emcee sampler
    args = {
        'ndim': ndim, 'rt_max': rt_max, 'cl_cent': cl_cent, 'fd': field_dens,
        'N_memb': N_memb, 'rt_rang': rt_rang, 'xy_in': xy_in, 'r_in': r_in}
    sampler = ensemble.EnsembleSampler(
        nchains, ndim, lnprob, kwargs=args, moves=mv)

    # Run the sampler hiding some annoying warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        tau_index, autocorr_vals, afs = 0, np.empty(nruns), np.empty(nruns)
        old_tau = np.inf
        for i, (pos, prob, stat) in enumerate(
                sampler.sample(pos0, iterations=nruns)):
            # Every X steps
            if i % KP_steps and i < (nruns - 1):
                continue

            afs[tau_index] = np.mean(sampler.acceptance_fraction)
            tau = sampler.get_autocorr_time(tol=0)
            autocorr_vals[tau_index] = np.mean(tau)
            tau_index += 1

            # Check convergence
            converged = np.all(tau * 100 < sampler.iteration)
            converged &= np.all(np.abs(old_tau - tau) / tau < 0.01)
            if converged:
                print("")
                break
            old_tau = tau
            update_progress.updt(nruns, i + 1)

    # Remove burn-in
    nburn = int(i * nburn)
    samples = sampler.get_chain(discard=nburn, flat=True)

    # Estimate mean, median, 16th, 84th percentiles for each parameter
    rc_m, rt_m = np.mean(samples[:, :2], 0)
    ecc_m, theta_m = 0., 0.
    rc_16, rc_50, rc_84, rt_16, rt_50, rt_84 = np.percentile(
        samples[:, :2], (16, 50, 84), 0).T.flatten()
    ecc_16, ecc_50, ecc_84, theta_16, theta_50, theta_84 =\
        [np.array([np.nan] * 3) for _ in range(6)]
    # Re-write if eccentricity and theta where obtained
    if ndim == 4:
        ecc_m = np.mean(samples[:, 2], 0)
        theta_m = circmean(samples[:, 3] * u.rad).value
        # WARNING: the median and percentiles for theta might are not
        # properly defined. Should use the 'circvar' function instead for
        # the variance, and the probably remove the percentiles
        ecc_16, ecc_50, ecc_84, theta_16, theta_50, theta_84 =\
            np.percentile(samples[:, 2:], (16, 50, 84), 0).T.flatten()

    # Parameters for plotting
    KP_plot = plotParams(
        KP_steps, ndim, sampler, samples, afs, autocorr_vals, tau_index, rc_m,
        rt_m, ecc_m, theta_m, rc_50, rt_50, ecc_50, theta_50, rt_max, cl_cent,
        field_dens, N_memb, xy_in, r_in, rt_rang)

    # Store: 16th, median, 84th, mean, mode
    KP_Bys_rc = np.array([rc_16, rc_50, rc_84, rc_m, KP_plot['KP_mode'][0]])
    KP_Bys_rt = np.array([rt_16, rt_50, rt_84, rt_m, KP_plot['KP_mode'][1]])
    KP_Bys_ecc = np.array([
        ecc_16, ecc_50, ecc_84, ecc_m, KP_plot['KP_mode'][2]])
    KP_Bys_theta = np.array([
        theta_16, theta_50, theta_84, theta_m, KP_plot['KP_mode'][3]])

    return KP_plot, KP_Bys_rc, KP_Bys_rt, KP_Bys_ecc, KP_Bys_theta


def lnprob(pars, ndim, rt_max, cl_cent, fd, N_memb, xy_in, r_in, rt_rang):
    """
    Logarithmic posterior
    """
    if ndim == 2:
        rc, rt = pars
        ecc, theta = 0., 0.
    elif ndim == 4:
        rc, rt, ecc, theta = pars

    # Prior.
    if rt <= rc or rc <= 0. or rt > rt_max or ecc < .0 or ecc > 1. or\
            theta < -np.pi / 2. or theta > np.pi / 2.:
        return -np.inf

    return lnlike(
        (rc, rt, ecc, theta), ndim, rt_max, cl_cent, fd, N_memb, xy_in, r_in,
        rt_rang)


def lnlike(
    pars, ndim, rt_max, cl_cent, fd, N_memb, xy_in, r_in, rt_rang,
        return_dens=False):
    """
    As defined in Pieres et al. (2016)
    """
    rc, rt, ecc, theta = pars

    if ndim == 2:
        # Values outside the tidal radius contribute 'fd' to the likelihood.
        r_in_clip = np.clip(r_in, a_min=0., a_max=rt)
        # N_in_region = r_in.size
        # Evaluate each star in King's profile
        KP = KingProf(r_in_clip, rc, rt)

    elif ndim == 4:
        x, y = xy_in
        # Identify stars inside this ellipse
        in_ellip_msk = inEllipse(xy_in, cl_cent, rt, ecc, theta)
        # N_in_region = in_ellip_msk.sum()

        # General form of the ellipse
        # https://math.stackexchange.com/a/434482/37846
        dx, dy = x - cl_cent[0], y - cl_cent[1]
        x1 = dx * np.cos(theta) + dy * np.sin(theta)
        y1 = dx * np.sin(theta) - dy * np.cos(theta)
        # The 'width' ('a') is used instead of the radius 'r' (sqrt(x**2+y**2))
        # in King's profile, for each star
        a_xy = np.sqrt(x1**2 + y1**2 / (1. - ecc**2))

        # Values outside the ellipse contribute only 'fd' to the likelihood.
        a_xy[~in_ellip_msk] = rt
        KP = KingProf(a_xy, rc, rt)

    # if N_memb is None:
    #     N_memb = NmembEst(ndim, fd, N_in_region, rt_max, rt, ecc)

    # Central density
    k = centDens(N_memb, rt_rang, rc, rt, ecc)
    if return_dens is True:
        return k

    # Likelihood
    li = k * KP + fd

    # Sum of log-likelighood
    sum_log_lkl = np.log(li).sum()

    return sum_log_lkl


# def NmembEst(ndim, fd, N_in_region, rt_max, rt, ecc):
#     """
#     The number of true members is estimated as the total number of stars
#     inside the region, minus the expected number of field stars in
#     the region.
#     """
#     if ndim == 2:
#         N_field_stars = np.pi * rt_max**2 * fd
#     elif ndim == 4:
#         a = rt
#         b = a * np.sqrt(1. - ecc**2)
#         N_field_stars = (np.pi * a * b) * fd

#     N_memb = N_in_region - N_field_stars
#     # Minimum value is 1
#     N_memb = max(1, N_memb)
#     return N_memb


def inEllipse(xy_in, cl_cent, rt, ecc, theta):
    """
    Source: https://stackoverflow.com/a/59718947/1391441
    """
    # Transpose
    xy = xy_in.T

    # The tidal radius 'rt' is made to represent the width ('a')
    # Width (squared)
    a2 = rt**2
    # Height (squared)
    b2 = a2 * (1. - ecc**2)

    # distance between the center and the foci
    foc_dist = np.sqrt(np.abs(b2 - a2))
    # vector from center to one of the foci
    foc_vect = np.array([foc_dist * np.cos(theta), foc_dist * np.sin(theta)])
    # the two foci
    el_foc1 = cl_cent + foc_vect
    el_foc2 = cl_cent - foc_vect

    # For each x,y: calculate z as the sum of the distances to the foci;
    # np.ravel is needed to change the array of arrays (of 1 element) into a
    # single array. Points are exactly on the ellipse when the sum of distances
    # is equal to the width
    z = np.ravel(np.linalg.norm(xy - el_foc1, axis=-1) +
                 np.linalg.norm(xy - el_foc2, axis=-1))

    # Mask that identifies the points inside the ellipse
    in_ellip_msk = z <= 2. * rt  # np.sqrt(max(a2, b2)) * 2.

    return in_ellip_msk


def centDens(N_memb, arr, rc, rt, ecc):
    """
    Central density constant. Integrate up to rt.

    https://math.stackexchange.com/a/1891110/37846
    """
    i = np.searchsorted(arr, rt)
    b = arr[:i] * np.sqrt(1. - ecc**2)
    integ = np.trapz(2. * np.pi * b * KingProf(arr[:i], rc, rt), arr[:i])
    return N_memb / integ


def KingProf(r_in, rc, rt):
    """
    King (1962) profile.
    """
    return ((1. / np.sqrt(1. + (r_in / rc) ** 2)) -
            (1. / np.sqrt(1. + (rt / rc) ** 2))) ** 2


def num_memb_conc_param(cd, rc, rt):
    """
    Calculate approximate number of cluster members with Eq (3) from
    Froebrich et al. (2007); 374, 399-408 and the concentration
    parameter.
    """

    # Approximate number of members.
    x = 1 + (rt / rc) ** 2
    n_c_k = int(round(
        (np.pi * cd * rc ** 2) * (
            np.log(x) - 4 + (4 * np.sqrt(x) + (x - 1)) / x)))
    # Concentration parameter.
    kcp = np.log10(rt / rc)

    return n_c_k, kcp


def plotParams(
    KP_steps, ndim, sampler, samples, afs, autocorr_vals, tau_index, rc_m,
    rt_m, ecc_m, theta_m, rc_50, rt_50, ecc_50, theta_50, rt_max, cl_cent,
        field_dens, N_memb, xy_in, r_in, rt_rang):
    """
    """
    # Mean acceptance fraction
    KP_mean_afs = afs[:tau_index]
    # Autocorrelation
    KP_tau_autocorr = autocorr_vals[:tau_index]
    # Effective sample size
    KP_ESS = samples.shape[0] / np.mean(sampler.get_autocorr_time(tol=0))
    # All samples, shape:(nsteps, nchains, ndim)
    KP_samples = sampler.get_chain()

    # Mode and KDE
    if ndim == 2:
        # This simulates the 'fundam_params and 'varIdxs' arrays.
        fp, vi = [[-np.inf, np.inf], [-np.inf, np.inf]], [0, 1]
        KP_mode, KP_kde = modeKDE(fp, vi, samples.T)
        # Add ecc, theta
        KP_mode += [0., 0.]
        KP_kde += [[], []]
    elif ndim == 4:
        fp, vi = [
            [-np.inf, np.inf], [-np.inf, np.inf], [-np.inf, np.inf],
            [-np.inf, np.inf]], [0, 1, 2, 3]
        KP_mode, KP_kde = modeKDE(fp, vi, samples.T)

    rc_MAD = stats.median_absolute_deviation(samples[:, 0])
    rt_MAD = stats.median_absolute_deviation(samples[:, 1])
    if ndim == 4:
        ecc_MAD = stats.median_absolute_deviation(samples[:, 2])
        theta_MAD = stats.median_absolute_deviation(samples[:, 3])

    # Central density. Use mean values for all the parameters.
    KP_cent_dens = lnlike(
        (rc_m, rt_m, ecc_m, theta_m), ndim, rt_max, cl_cent, field_dens,
        N_memb, xy_in, r_in, rt_rang, True)

    # 16th-84th percentile region for the profile fit
    ecc, theta = 0., 0.
    kpf_yvals = []
    cent_dens_all = []
    for _ in range(1000):
        # Sample rc, rt, ecc, theta. Use the median and MAD.
        rc = np.random.normal(rc_50, 1.4826 * rc_MAD)
        rt = np.random.normal(rt_50, 1.4826 * rt_MAD)
        if ndim == 4:
            ecc = np.random.normal(ecc_50, 1.4826 * ecc_MAD)
            theta = np.random.normal(theta_50, 1.4826 * theta_MAD)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # Obtain central density
            KP_cd_ = lnlike(
                (rc, rt, ecc, theta), ndim, rt_max, cl_cent, field_dens,
                N_memb, xy_in, r_in, rt_rang, True)
            cent_dens_all.append(KP_cd_)
        # Values in y
        kpf_yvals.append(KP_cd_ * KingProf(rt_rang, rc, rt) + field_dens)

    cent_dens_all = np.array(cent_dens_all) / 3600.
    print(np.mean(cent_dens_all), np.median(cent_dens_all), np.std(cent_dens_all))
    kpf_yvals = np.array(kpf_yvals).T
    _16_kp = np.nanpercentile(kpf_yvals, 16, 1)
    _84_kp = np.nanpercentile(kpf_yvals, 84, 1)

    KP_plot = {
        'KP_steps': KP_steps,
        'KP_mean_afs': KP_mean_afs, 'KP_tau_autocorr': KP_tau_autocorr,
        'KP_ESS': KP_ESS, 'KP_samples': KP_samples, 'KP_mode': KP_mode,
        'KP_kde': KP_kde, 'KP_cent_dens': KP_cent_dens,
        '_16_84_rang': rt_rang, '_16_kp': _16_kp, '_84_kp': _84_kp}

    return KP_plot


#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
# All the code below is for testing purposes only

def test():
    """
    """
    seed = np.random.randint(10000)
    print(seed, "\n")
    np.random.seed(seed)

    Nf, Ncl = 2000, 200
    cl_cent = (1000., 1000.)
    frame_rang = 2000.
    xy_field = np.random.uniform(0., frame_rang, (2, Nf))
    fd = Nf / frame_rang**2
    N_integ = 1000

    def neglnprob(pars, ndim, rt_max, cl_cent, fd, Ncl, xy_in, r_in, rt_rang):
        return -1. * lnprob(
            pars, ndim, rt_max, cl_cent, fd, Ncl, xy_in, r_in, rt_rang)

    orig, res = [], []
    for _ in range(10):
        # rt
        width = np.random.uniform(50., 300.)
        ecc = .8 # np.random.uniform(.5, .9)
        theta = 0. # np.random.uniform(-np.pi / 2., np.pi / 2.)
        # rc
        height = width * (1. - ecc**2)
        print(
            ("\nOrig r_c, r_t, ecc, theta: {:.2f}, {:.2f}, {:.2f}," +
             " {:.1f}").format(height, width, ecc, np.rad2deg(theta)))
        orig.append([height, width, ecc, theta])

        N_in_ellip = 0.
        xy_clust = []
        while N_in_ellip < Ncl:
            # Generate stars positions sampling from a King profile
            # rc, rt = height, width
            cl_dists = invTrnsfSmpl(height, width, ecc, theta, 1000)
            phi = np.random.uniform(0., 1., 1000) * 2 * np.pi
            x = cl_cent[0] + cl_dists * np.cos(phi)
            y = cl_cent[1] + cl_dists * np.sin(phi)
            xy = np.array([x, y])

            # Filter stars outside the ellipse
            msk = inEllipse(xy, cl_cent, width, ecc, theta)

            ax = plt.subplot(111)
            plt.scatter(*xy, c='k')
            ellipse = mpatches.Ellipse(
                xy=cl_cent, width=2. * width, height=2. * width,
                angle=0.,
                facecolor='None', edgecolor='black', linewidth=2,
                transform=ax.transData, zorder=6)
            ax.add_patch(ellipse)

            plt.scatter(*xy[:, msk], c='g')
            a2 = width**2
            b2 = a2 * (1. - ecc**2)
            ellipse = mpatches.Ellipse(
                xy=cl_cent, width=2. * np.sqrt(a2), height=2. * np.sqrt(b2),
                angle=np.rad2deg(theta),
                facecolor='None', edgecolor='green', linewidth=2,
                transform=ax.transData, zorder=6)
            ax.add_patch(ellipse)

            plt.show()
            import pdb; pdb.set_trace()  # breakpoint f7c56dc5 //


            xy_clust += list(xy.T[msk])
            N_in_ellip = len(xy_clust)
        xy_clust = np.array(xy_clust)[:Ncl].T

        # Full frame
        xy = np.concatenate([xy_field.T, xy_clust.T])
        rt_max = width * 3
        rt_rang = np.linspace(0., rt_max, N_integ)

        # Stars inside the cut-out given by the 'rt_max' value.
        xy_cent_dist = spatial.distance.cdist([cl_cent], xy)[0]
        msk = xy_cent_dist <= rt_max
        r_in = xy_cent_dist[msk]
        xy_in = xy[msk].T

        ndim = 4
        # rc, rt, ecc, theta
        bounds = ((.05 * rt_max, rt_max), (.05 * rt_max, rt_max),
                  (.5, .9), (-np.pi / 2., np.pi / 2.))
        result = differential_evolution(
            neglnprob, bounds,
            args=(ndim, rt_max, cl_cent, fd, Ncl, xy_in, r_in, rt_rang),
            maxiter=2000, popsize=50, tol=0.00001)
        r_rc, r_rt, r_ecc, r_theta = result.x
        print(
            ("Recov r_c, r_t, ecc, theta: {:.2f}, {:.2f}, {:.2f}, " +
             "{:.1f}").format(r_rc, r_rt, r_ecc, np.rad2deg(r_theta)))
        res.append([r_rc, r_rt, r_ecc, r_theta])

        theta_p = np.pi + theta if theta < 0. else theta
        r_theta_p = np.pi + r_theta if r_theta < 0. else r_theta
        diff1 = np.rad2deg(abs(theta_p - r_theta_p))
        diff2 = np.rad2deg(abs(theta - r_theta))
        print(min(diff1, diff2))

        print(result)
        resPlot(cl_cent, xy_field, xy_clust, height, width, ecc, theta,
                *result.x)

    orig, res = np.array(orig), np.array(res)
    perc_diff = (100. * (orig[:, :2] - res[:, :2]) / orig[:, :2]).T

    plt.subplot(221)
    plt.title("rc")
    plt.scatter(orig[:, 0], perc_diff[0])
    plt.subplot(222)
    plt.title("rt")
    plt.scatter(orig[:, 1], perc_diff[1])
    plt.subplot(223)
    plt.title("ecc")
    plt.scatter(orig[:, 2], orig[:, 2] - res[:, 2])
    plt.subplot(224)
    plt.title("theta")
    plt.scatter(orig[:, 3], np.rad2deg(abs(orig[:, 3] - res[:, 3])))
    plt.show()


def resPlot(
    cl_cent, xy_field, xy_clust, rc, rt, ecc, theta, r_rc, r_rt, r_ecc,
        r_theta):
    """
    """

    ax = plt.subplot(111)
    a2 = rt**2
    b2 = a2 * (1. - ecc**2)
    ellipse = mpatches.Ellipse(
        xy=cl_cent, width=2. * np.sqrt(a2), height=2. * np.sqrt(b2),
        angle=np.rad2deg(theta),
        facecolor='None', edgecolor='green', linewidth=2, ls='--',
        transform=ax.transData, zorder=6)
    ax.add_patch(ellipse)

    a2 = r_rt**2
    b2 = a2 * (1. - r_ecc**2)
    ellipse = mpatches.Ellipse(
        xy=cl_cent, width=2. * np.sqrt(a2), height=2. * np.sqrt(b2),
        angle=np.rad2deg(r_theta),
        facecolor='None', edgecolor='black', linewidth=2,
        transform=ax.transData, zorder=6)
    ax.add_patch(ellipse)

    plt.scatter(*xy_field, c='k', s=15, alpha=.5)
    plt.scatter(*xy_clust, c='g', s=15, zorder=5)
    # plt.scatter(*xy_in, c='r', s=5, zorder=1)
    plt.axvline(1000.)
    plt.axhline(1000.)
    plt.scatter(*cl_cent, marker='x', c='r', s=30, zorder=7)
    plt.xlim(cl_cent[0] - 2. * rt, cl_cent[0] + 2. * rt)
    plt.ylim(cl_cent[1] - 2. * rt, cl_cent[1] + 2. * rt)
    plt.show()


def invTrnsfSmpl(rc, rt, ecc, theta, N_samp, N_interp=1000):
    """
    Sample King's profile using the inverse CDF method.
    """
    from scipy.integrate import quad
    from scipy.interpolate import interp1d

    def rKP(r, rc, rt, ecc, theta):
        x, y = r * np.cos(theta), r * np.sin(theta)
        x1 = x * np.cos(theta) + y * np.sin(theta)
        y1 = x * np.sin(theta) - y * np.cos(theta)
        # The 'width' ('a') is used instead of the radius 'r' (sqrt(x**2+y**2))
        # in King's profile, for each star
        r = np.sqrt(x1**2 + y1**2 / (1. - ecc**2))

        return r * KingProf(r, rc, rt)

    r_0rt = np.linspace(0., rt, N_interp)
    # The CDF is defined as: $F(r)= \int_{r_low}^{r} PDF(r) dr$
    # Sample the CDF
    CDF_samples = []
    for r in r_0rt:
        CDF_samples.append(quad(rKP, 0., r, args=(rc, rt, ecc, theta))[0])

    # Normalize CDF
    CDF_samples = np.array(CDF_samples) / CDF_samples[-1]

    # Inverse CDF
    inv_cdf = interp1d(CDF_samples, r_0rt)

    # Sample the inverse CDF
    samples = inv_cdf(np.random.rand(N_samp))

    return samples


if __name__ == '__main__':
    from scipy.optimize import differential_evolution
    import matplotlib.pyplot as plt
    from matplotlib import patches as mpatches
    test()
