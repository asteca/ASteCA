
import numpy as np
from astropy.stats import circmean
from astropy import units as u
from scipy import spatial
from scipy import stats
try:
    # If this import is not done outside main(), then eval() fails in the
    # definition of the moves
    from emcee import moves
except ImportError:
    pass
import warnings
from .. import update_progress
from ..out import prep_plots
from ..best_fit.bf_common import modeKDE

"""
I attempted to produce a King Profile fit setting the four parameters free:
central density, rc, rt, and field density. Then I use the per-star densities
obtained previously, and compare those with the King densities for each star.

This is, for each star with a distance to center and local density fixed
(r, d):
1. find its King density given some values for its free parameters
(rho0, rc, rt, fd) --> d_K
2. assume a Gaussian distribution for the densities and obtain the likelihood
as ~ exp(-(d-d_K)**2)
3. multiply all exponentials for the full likelihood (apply logarithm to
use a sum of negative squared differences)

I only tried with (ecc, theta) fixed to zero.

What I found was:

1. The number of neighbors used to estimate the local densities can affect the
fit *a lot*.

2. The differential_evolution method finds good fits but tends to overestimate
the rt, and underestimate the rc.

3. Even with reasonable initialization points from the DE, emcee tends to
(rc-->0, rt-->inf) and both the central density and the field density are
almost not taken into account (i.e., their distributions are uniform in the
range given). I don't really know why this happens.

"""


def main(
    clp, cld_i, kp_ndim, kp_nchains, kp_nruns, kp_nburn, kp_emcee_moves,
        rt_max_f, coords, **kwargs):
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
    if kp_ndim in (2, 4):
        print("Estimating King profile")
        KP_plot, KP_Bys_rc, KP_Bys_rt, KP_Bys_ecc, KP_Bys_theta =\
            fit_King_prof(
                kp_ndim, kp_nchains, kp_nruns, kp_nburn, kp_emcee_moves,
                cld_i['x'], cld_i['y'], clp['kde_cent'], clp['fr_dens'],
                clp['field_dens'], clp['clust_rad'], rt_max_f)

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
    ndim, nchains, nruns, nburn, kp_emcee_moves, x, y, cl_cent, fr_dens,
        field_dens, cl_rad, rt_max_f, N_integ=1000, N_conv=500, tau_stable=0.01):
    """
    N_integ : number of points used in the integration performed in the
    central density function

    N_conv, tau_stable: if the chain is longer than N_conv times the
    estimated autocorrelation time and if this estimate changed by less than
    (tau_stable * 100)%

    ndim = 4 fits (rc, rt, fd, rho0,)
    ndim = 6 fits (rc, rt, fd, rho0, ecc, theta)

    ** TODO **

    * Generalize to exponential profile? See:
    https://www.aanda.org/articles/aa/abs/2016/02/aa27070-15/aa27070-15.html
    Eq (2),
    https://ned.ipac.caltech.edu/level5/Sept11/Graham/Graham2.html,
    Eq (1)

    * Align the theta with the y axis (North) instead of the x axis (as done
    in Martin et al. 2008; 'A Comprehensive Maximum Likelihood Analysis of
    the Structural Properties of Faint Milky Way Satellites')?
    """

    from emcee import ensemble
    # Move used by emcee
    mv = [(eval("(moves." + _ + ")")) for _ in kp_emcee_moves]

    # Steps to store Bayes params.
    KP_steps = int(nruns * .01)

    # The tidal radius can not be larger than 'rt_max' times the estimated
    # "optimal" cluster radius. Used as a prior.
    rt_max = rt_max_f * cl_rad
    fd_max = field_dens * 2
    rho0_max = (max(fr_dens) - field_dens) * 2.

    # Tidal radius array. Used for integrating
    rt_rang = np.linspace(0., rt_max, int(rt_max_f * N_integ))

    # Identify stars inside the cut-out given by the 'rt_max' value. Only these
    # stars will be processed below.
    xy = np.array((x, y)).T
    xy_cent_dist = spatial.distance.cdist([cl_cent], xy)[0]
    msk = xy_cent_dist <= rt_max
    xy_in = xy[msk].T
    r_in = xy_cent_dist[msk]
    fr_dens_in = fr_dens[msk]


    from scipy.optimize import differential_evolution

    def nll(pars, *args):
        return -lnprob(pars, *args)
    # # rho0 = (fr_dens_in.max() - field_dens) / (1 - 1. / np.sqrt(1 + 2**2))**2
    bounds = [(0, rt_max), (0, rt_max), (0, fd_max), (0, rho0_max)]
    soln = differential_evolution(nll, bounds, args=(
        cl_cent, ndim, rt_max, fd_max, rho0_max, fr_dens_in, xy_in, r_in))
    rc0, rt0, fd0, rho00 = soln.x
    print(soln)

    # def gridfunc(Nv, rcax, rtax, fdax, rhoax):
    #     gridvals = []
    #     # parvals = []
    #     for rc in rcax:
    #         for rt in rtax:
    #             for fd in fdax:
    #                 for rho in rhoax:
    #                     pars = (rc, rt, fd, rho)
    #                     # parvals.append(pars)
    #                     gridvals.append(-nll(pars, *(
    #                         cl_cent, ndim, rt_max, fd_max, rho0_max,
    #                         fr_dens_in, xy_in, r_in)))
    #     return np.array(gridvals).reshape([Nv, Nv, Nv, Nv])


    # from numpy import unravel_index
    # import matplotlib.pyplot as plt
    # def callgrid(bounds, Nv):
    #     rcax = np.linspace(bounds[0][0], bounds[0][1], Nv)
    #     rtax = np.linspace(bounds[1][0], bounds[1][1], Nv)
    #     fdax = np.linspace(bounds[2][0], bounds[2][1], Nv)
    #     rhoax = np.linspace(bounds[3][0], bounds[3][1], Nv)
    #     gridvals = gridfunc(Nv, rcax, rtax, fdax, rhoax)
    #     max_i = unravel_index(gridvals.argmax(), gridvals.shape)
    #     print(gridvals.max())
    #     print(rcax[max_i[0]], rtax[max_i[1]], fdax[max_i[2]], rhoax[max_i[3]])

    #     return rcax, rtax, fdax, rhoax, gridvals

    # Nv = 10
    # bounds = [(100, 220), (220, 400), (0, fd_max), (0, rho0_max)]
    # rcax, rtax, fdax, rhoax, gridvals = callgrid(bounds, Nv)

    # plt.imshow(gridvals[:,:,5,5],vmin=-.05,vmax=0);plt.colorbar();plt.show()


    # Initial positions for the sampler.
    # Dimensions: rc, rt
    # rc_pos0 = np.random.uniform(.05 * rt_max, rt_max, nchains)
    # rt_pos0 = np.random.uniform(rc_pos0, rt_max, nchains)
    rc_pos0 = rc0 + (rc0 * .1) * np.random.randn(nchains)
    rt_pos0 = rt0 + (rt0 * .1) * np.random.randn(nchains)
    pos0 = [rc_pos0, rt_pos0]
    if ndim == 4:
        # Dimensions: ecc, theta
        # fd_pos0 = np.random.uniform(.0, fd_max, nchains)
        # rho0_pos0 = np.random.uniform(.0, rho0_max, nchains)
        fd_pos0 = fd0 + (fd0 * .1) * np.random.randn(nchains)
        rho0_pos0 = rho00 + (rho00 * .1) * np.random.randn(nchains)
        pos0 += [fd_pos0, rho0_pos0]
    pos0 = np.array(pos0).T

    # Define emcee sampler
    args = {
        'cl_cent': cl_cent, 'ndim': ndim, 'rt_max': rt_max,
        'fd_max': fd_max, 'rho0_max': rho0_max, 'fr_dens_in': fr_dens_in,
        'xy_in': xy_in, 'r_in': r_in}
    sampler = ensemble.EnsembleSampler(
        nchains, ndim, lnprob, kwargs=args, moves=mv)

    # Run the sampler hiding some annoying warnings
    conver_flag = False
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
                conver_flag = True
                break
            old_tau = tau
            # update_progress.updt(nruns, i + 1)

    if conver_flag:
        print("Process converged")

    # Remove burn-in
    nburn = int(i * nburn)
    samples = sampler.get_chain(discard=nburn, flat=True)

    import matplotlib.pyplot as plt

    def plotKP(rc, rt, fd, rho0):
        # rc, rt, fd, rho0 = rc_m, rt_m, fd_m, rho0_m
        r_in_clip = np.clip(r_in, a_min=0., a_max=rt)
        KP = KingProf(r_in_clip, rc, rt)
        li = rho0 * KP + fd
        plt.scatter(r_in, fr_dens_in, c='g', alpha=.5)
        plt.scatter(r_in, li, c='r', alpha=.5)
        plt.show()

    def rho0M(KP_0, rc, rt, fd, rho0=None):
        if rho0 is None:
            rho0 = (KP_0 - fd) / (1 - 1. / np.sqrt(1 + (rt/rc)**2))**2
        pars = (rc, rt, fd, rho0, 0., 0.)
        print(lnlike(pars, ndim, cl_cent, xy_in, r_in, fr_dens_in))
        return rho0

    rc, rt, fd, rho0 = np.mean(samples, 0)
    print('\n Final vals', rc, rt, fd, rho0)
    rho0 = rho0M(np.nan, rc, rt, fd, rho0)
    plotKP(rc, rt, fd, rho0)

    # # KP_0, rc, rt, fd = 1.2e-2, 200, 250, 9.5e-4
    # rc, rt, fd, rho0 = soln.x
    # print(rc, rt, fd, rho0)
    # rho0 = rho0M(np.nan, rc, rt, fd, rho0)
    # plotKP(rc, rt, fd, rho0)

    import corner
    figure = corner.corner(samples)
    plt.show()

    import pdb; pdb.set_trace()  # breakpoint b7998ac2 //


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


def lnprob(
    pars, cl_cent, ndim, rt_max, fd_max, rho0_max, fr_dens_in, xy_in, r_in):
    """
    Logarithmic posterior
    """
    if ndim == 4:
        rc, rt, fd, rho0 = pars
        ecc, theta = 0., 0.
    elif ndim == 6:
        rc, rt, fd, rho0, ecc, theta = pars

    # Prior.
    # ecc < .0 or ecc > 1. or\
    # theta < -np.pi / 2. or theta > np.pi / 2.:
    if rt <= rc or rc <= 0. or rt > rt_max or fd < 0. or fd > fd_max or\
            rho0 < 0. or rho0 > rho0_max:
        return -np.inf

    return lnlike(
        (rc, rt, fd, rho0, ecc, theta), ndim, cl_cent, xy_in, r_in,
        fr_dens_in)


def lnlike(pars, ndim, cl_cent, xy_in, r_in, fr_dens_in, return_dens=False):
    """
    As defined in Pieres et al. (2016)
    """
    rc, rt, fd, rho0, ecc, theta = pars

    if ndim == 4:
        # Values outside the tidal radius contribute 'fd' to the likelihood.
        r_in_clip = np.clip(r_in, a_min=0., a_max=rt)

        # Evaluate each star in King's profile
        KP = KingProf(r_in_clip, rc, rt)

    elif ndim == 6:
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

        # Values outside the ellipse contribute 'fd' to the likelihood.
        a_xy[~in_ellip_msk] = rt
        KP = KingProf(a_xy, rc, rt)

    # Likelihood
    li = rho0 * KP + fd
    # Density delta
    dens = np.square(fr_dens_in - li)
    # Sum of log-likelihood
    sum_log_lkl = dens.sum()
    if np.random.randint(100) == 50:
        print(rt, rho0, fd, sum_log_lkl)
    #     # print(rc, rt, fd, rho0)
    #     import matplotlib.pyplot as plt
    #     plt.title("{:.5f}, {:.0f}".format(sum_log_lkl, rt))
    #     plt.scatter(r_in_clip, dens, alpha=.5)
    #     plt.ylim(0., max(dens))
    # #     plt.scatter(r_in_clip,li,c='r')
    # #     plt.scatter(r_in_clip,fr_dens_in,c='g')
    #     plt.show()

    # sum_log_lkl = -np.inf if np.isnan(sum_log_lkl) else sum_log_lkl
    return sum_log_lkl


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
    z = np.ravel(np.linalg.norm(xy - el_foc1, axis=-1)
                 + np.linalg.norm(xy - el_foc2, axis=-1))

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
    if ndim in (2, 3):
        # This simulates the 'fundam_params and 'varIdxs' arrays.
        if ndim == 2:
            fp, vi = [[-np.inf, np.inf], [-np.inf, np.inf]], [0, 1]
        else:
            fp, vi = [
                [-np.inf, np.inf], [-np.inf, np.inf], [-np.inf, np.inf]],\
                [0, 1, 2]
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
