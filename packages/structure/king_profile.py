
import numpy as np
from astropy.stats import circmean, circstd
from astropy import units as u
from scipy import stats
import warnings
from .. import update_progress
from ..best_fit.bf_common import modeKDE


def main(
    clp, cld, kp_ndim, kp_nchains, kp_nruns, kp_nburn,
        rt_max_f, **kwargs):
    """
    Bayesian inference over an array of stars' coordinates using a King
    profile model.

    The previously obtained 'clust_rad' value is used to determine the
    range where rt should be fitted.

    The field density value is fixed to the previously estimated value.
    The central density is inferred from the (rc, rt) values and the previously
    estimated number of members ('n_memb', which is *also* estimated
    using the field density).

    This means that the value given to the field density has a *very large*
    influence on the final (rc, rt) values.

    IMPORTANT: the function will not work properly if the cluster's area
    is cropped; i.e.: near a frame's border.

    """
    if kp_ndim in (2, 4):
        print("Estimating King profile")
        KP_plot, KP_Bys_rc, KP_Bys_rt, KP_Bys_ell, KP_Bys_theta =\
            fit_King_prof(
                kp_ndim, kp_nchains, kp_nruns, kp_nburn,
                clp['kde_cent'], clp['xy_filtered'], clp['xy_cent_dist'],
                clp['field_dens'], clp['clust_rad'], clp['n_memb'], rt_max_f)

        # Obtain number of members and concentration parameter.
        if KP_Bys_rt['median'] > 0.:
            arr = np.linspace(0, KP_Bys_rt['median'], 1000)
            KP_memb_num = KP_memb_x(
                KP_plot['KP_cent_dens'], KP_Bys_rc['median'],
                KP_Bys_rt['median'], KP_Bys_ell['median'], arr,
                KP_Bys_rt['median'])
            KP_conct_par = np.log10(KP_Bys_rt['median'] / KP_Bys_rc['median'])
        else:
            KP_memb_num, KP_conct_par = np.nan, np.nan

        # Set precision of printed values.
        print('Core & tidal radii obtained: {:g}, {:g} deg'.format(
            KP_Bys_rc['median'], KP_Bys_rt['median']))

    else:
        print("Skipping King profile fit")
        KP_plot, KP_memb_num, KP_conct_par = {}, np.nan, np.nan
        KP_Bys_rc, KP_Bys_rt, KP_Bys_ell, KP_Bys_theta =\
            [{'16th': np.nan, 'median': np.nan, '84th': np.nan, 'mean': np.nan,
              'std': np.nan, 'mode': np.nan} for _ in range(4)]

    clp['KP_plot'], clp['KP_Bys_rc'], clp['KP_Bys_rt'],\
        clp['KP_Bys_ell'], clp['KP_Bys_theta'], clp['KP_memb_num'],\
        clp['KP_conct_par'] = KP_plot, KP_Bys_rc, KP_Bys_rt, KP_Bys_ell,\
        KP_Bys_theta, KP_memb_num, KP_conct_par
    return clp


def fit_King_prof(
    ndim, nchains, nruns, nburn, cl_cent, xy_filtered,
    xy_cent_dist, field_dens, cl_rad, n_memb, rt_max_f, N_integ=1000,
        N_conv=500, tau_stable=0.01):
    """
    N_integ : number of points used in the integration performed in the
    central density function

    N_conv, tau_stable: if the chain is longer than N_conv times the
    estimated autocorrelation time and if this estimate changed by less than
    (tau_stable * 100)%

    ndim = 1 fits (rc)                 <-- Not implemented
    ndim = 2 fits (rc, rt)
    ndim = 3 fits (rc, rt, fd)         <-- Not implemented
    ndim = 4 fits (rc, rt, ell, theta)

    ## About the 'N_memb' and 'field_dens' parameters

    * 'N_memb' fixed + 'fd' fixed: best results
    * 'N_memb' free  + 'fd' fixed: N_memb grows and very large tidal radii
       are favored
    * 'N_memb' free  + 'fd' free : the field density tends to the maximum
      allowed, which pulls the tidal radius to lower values, The number of
      members also drops to small values.
    * 'N_memb' fixed  + 'fd' free: the field density tends to the maximum
      allowed, which pulls the tidal radius to lower values


    ** TODO **

    * Fit using the RDP and a chi-square (Kuepper et al. 2010, Tarricq et
    al. 2021)

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
    # DESnookerMove(), 0.1; DEMove(), 0.9 * 0.9; DEMove(gamma0=1.0), 0.9 * 0.1
    # mv = [(eval("(moves." + _ + ")")) for _ in kp_emcee_moves]

    # Steps to store Bayes params.
    KP_steps = max(1, int(nruns * .01))

    # Fixed number of members (previously estimated)
    N_memb = 1 * n_memb

    # The tidal radius can not be larger than 'rt_max' times the estimated
    # "optimal" cluster radius. Used as a prior.
    rt_max = rt_max_f * cl_rad
    # Tidal radius array. Used for integrating
    rt_rang = np.linspace(0., rt_max, int(rt_max_f * N_integ))

    # Identify stars inside the cut-out given by the 'rt_max' value. Only these
    # stars will be processed below.
    msk = xy_cent_dist <= rt_max
    xy_in = xy_filtered[msk].T
    r_in = xy_cent_dist[msk]

    # Initial positions for the sampler.
    # Dimensions: rc, rt
    rc_pos0 = np.random.uniform(.05 * rt_max, rt_max, nchains)
    rt_pos0 = np.random.uniform(rc_pos0, rt_max, nchains)
    pos0 = [rc_pos0, rt_pos0]
    if ndim == 4:
        # Dimensions: ell, theta
        ell = np.random.uniform(.0, 1., nchains)
        theta = np.random.uniform(-np.pi / 2., np.pi / 2., nchains)
        pos0 += [ell, theta]
    pos0 = np.array(pos0).T

    # Define emcee sampler
    args = {
        'ndim': ndim, 'rt_max': rt_max, 'cl_cent': cl_cent, 'fd': field_dens,
        # 'fd_max': fd_max,
        'N_memb': N_memb, 'rt_rang': rt_rang, 'xy_in': xy_in, 'r_in': r_in}
    sampler = ensemble.EnsembleSampler(nchains, ndim, lnprob, kwargs=args)

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
            update_progress.updt(nruns, i + 1)

    if conver_flag:
        print("Process converged")

    # Remove burn-in
    nburn = int(i * nburn)
    samples = sampler.get_chain(discard=nburn, flat=True)

    # Estimate mean, median, 16th, 84th percentiles for each parameter
    if ndim in (2, 4):
        rc_m, rt_m = np.mean(samples[:, :2], 0)
        rc_std, rt_std = np.std(samples[:, :2], 0)
        rc_16, rc_50, rc_84, rt_16, rt_50, rt_84 = np.percentile(
            samples[:, :2], (16, 50, 84), 0).T.flatten()
    ell_m, ell_std, theta_m, theta_std = 0., 0., 0., 0
    ell_16, ell_50, ell_84, theta_std = [np.nan] * 4
    # Re-write if ellipticity and theta where obtained
    if ndim == 4:
        ell_m = np.mean(samples[:, 2], 0)
        ell_std = np.std(samples[:, 2], 0)
        ell_16, ell_50, ell_84 = np.percentile(
            samples[:, 2], (16, 50, 84), 0).T.flatten()
        theta_m = circmean(samples[:, 3] * u.rad).value
        theta_std = circstd(samples[:, 3] * u.rad).value

    # Parameters for plotting
    KP_plot = plotParams(
        KP_steps, ndim, sampler, samples, afs, autocorr_vals, tau_index, rc_m,
        rt_m, ell_m, theta_m, theta_std, rc_50, rt_50, ell_50, rt_max, cl_cent,
        field_dens, N_memb, xy_in, r_in, rt_rang)

    # Store: 16th, median, 84th, mean, STDDEV, mode
    KP_Bys_rc = {'16th': rc_16, 'median': rc_50, '84th': rc_84, 'mean': rc_m,
                 'std': rc_std, 'mode': KP_plot['KP_mode'][0]}
    KP_Bys_rt = {'16th': rt_16, 'median': rt_50, '84th': rt_84, 'mean': rt_m,
                 'std': rt_std, 'mode': KP_plot['KP_mode'][1]}
    KP_Bys_ell = {'16th': ell_16, 'median': ell_50, '84th': ell_84,
                  'mean': ell_m, 'std': ell_std, 'mode': KP_plot['KP_mode'][2]}
    KP_Bys_theta = {'16th': np.nan, 'median': np.nan, '84th': np.nan,
                    'mean': theta_m, 'std': theta_std,
                    'mode': KP_plot['KP_mode'][3]}

    return KP_plot, KP_Bys_rc, KP_Bys_rt, KP_Bys_ell, KP_Bys_theta


def lnprob(
        pars, ndim, rt_max, cl_cent, fd, N_memb, rt_rang, xy_in, r_in):
    """
    Logarithmic posterior
    """
    if ndim == 2:
        rc, rt = pars
        ell, theta = 0., 0.
    elif ndim == 4:
        rc, rt, ell, theta = pars

    # Prior.
    if rt <= rc or rc <= 0. or rt > rt_max or ell < .0 or ell >= 1. or\
            theta < -np.pi / 2. or theta > np.pi / 2.:
        return -np.inf

    return lnlike(
        (rc, rt, ell, theta), ndim, rt_max, cl_cent, fd, N_memb, xy_in, r_in,
        rt_rang)


def lnlike(
    pars, ndim, rt_max, cl_cent, fd, N_memb, xy_in, r_in, rt_rang,
        return_dens=False):
    """
    As defined in Pieres et al. (2016)
    """
    rc, rt, ell, theta = pars

    if ndim in (1, 2, 3):
        # Values outside the tidal radius contribute 'fd' to the likelihood.
        r_in_clip = np.clip(r_in, a_min=0., a_max=rt)

        # Evaluate each star in King's profile
        KP = KingProf(r_in_clip, rc, rt)
        # N_in_region = (r_in <= rt).sum()

    elif ndim == 4:
        x, y = xy_in
        # Identify stars inside this ellipse
        in_ellip_msk = inEllipse(xy_in, cl_cent, rt, ell, theta)
        # N_in_region = in_ellip_msk.sum()

        # General form of the ellipse
        # https://math.stackexchange.com/a/434482/37846
        dx, dy = x - cl_cent[0], y - cl_cent[1]
        x1 = dx * np.cos(theta) + dy * np.sin(theta)
        y1 = dx * np.sin(theta) - dy * np.cos(theta)
        # The 'width' ('a') is used instead of the radius 'r' (sqrt(x**2+y**2))
        # in King's profile, for each star
        a_xy = np.sqrt(x1**2 + y1**2 / (1. - ell)**2)

        # Values outside the ellipse contribute 'fd' to the likelihood.
        a_xy[~in_ellip_msk] = rt
        r_in_clip = a_xy

        KP = KingProf(r_in_clip, rc, rt)

    # Central density
    rho_0 = centDens(N_memb, rt_rang, rc, rt, ell)

    # Testing to avoid using 'N_memb'. A small number must be added to 'fd'
    # and still it does not work as well as the above method
    # rho_0 = (N_in_region - np.pi * rt**2 * fd * (1 + 0.025)) / KP_memb_x(1, rc, rt, rt)

    if return_dens is True:
        return rho_0

    # Likelihood
    # 'fd' IS necessary
    li = rho_0 * KP + fd
    # Sum of log-likelihood
    sum_log_lkl = np.log(li).sum()

    sum_log_lkl = -np.inf if np.isnan(sum_log_lkl) else sum_log_lkl
    return sum_log_lkl


def inEllipse(xy_in, cl_cent, rt, ell, theta):
    """
    Source: https://stackoverflow.com/a/59718947/1391441
    """
    # Transpose
    xy = xy_in.T

    # The tidal radius 'rt' is made to represent the width ('a')
    # Width (squared)
    a2 = rt**2
    # Height (squared)
    b2 = a2 * (1. - ell)**2

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


def centDens(N_memb, arr, rc, rt, ell):
    """
    Central density constant. Integrate up to rt.

    https://math.stackexchange.com/a/1891110/37846

    The number of members is estimated with:

    N_memb = cd * integ
    """
    # Avoid dividing by zero
    integ = max(0.0000001, KP_integ_x(rc, rt, ell, arr, rt))
    return N_memb / integ

    # This does not work, r_t grows to the maximum allowed value
    # cent_dens = some value
    # A = 1 / np.sqrt(1 + (rt / rc)**2)
    # rho_0 = cent_dens / (1 - A)**2
    # return rho_0


def KingProf(r_in, rc, rt):
    """
    King (1962) profile.
    """
    return ((1. / np.sqrt(1. + (r_in / rc) ** 2))
            - (1. / np.sqrt(1. + (rt / rc) ** 2))) ** 2


def KP_integ_x(rc, rt, ell, arr, rad):
    """
    Integral of King's profile up to rad. It must be understood as:

    N_memb = cd * integ

    where 'integ' is:

    x = 1 + (rt / rc) ** 2
    integ = (np.pi * rc ** 2) * (
            np.log(x) - 4 + (4 * np.sqrt(x) + (x - 1)) / x)

    i.e., the integral of 2*pi*r*KP divided by 'cd'. Below, 'integ' is
    equivalent to the one given above, but generalized to an ellipse.
    """
    i = np.searchsorted(arr, rad)
    b = arr[:i] * (1. - max(0, ell))  # The max() catches 'nan'
    integ = np.trapz(2. * np.pi * b * KingProf(arr[:i], rc, rt), arr[:i])
    return integ


def KP_memb_x(cd, rc, rt, ell, arr, rad):
    """
    Calculate approximate number of cluster members for a given King Profile,
    from 0 up to x.

    Used by packages/out/mp_radius.py
    """
    integ = KP_integ_x(rc, rt, ell, arr, rad)
    N_memb_x = int(cd * integ)

    return N_memb_x


def plotParams(
    KP_steps, ndim, sampler, samples, afs, autocorr_vals, tau_index, rc_m,
    rt_m, ell_m, theta_m, theta_std, rc_50, rt_50, ell_50, rt_max, cl_cent,
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
        # Add ell, theta
        KP_mode += [0., 0.]
        KP_kde += [[], []]
    elif ndim == 4:
        fp, vi = [
            [-np.inf, np.inf], [-np.inf, np.inf], [-np.inf, np.inf],
            [-np.inf, np.inf]], [0, 1, 2, 3]
        KP_mode, KP_kde = modeKDE(fp, vi, samples.T)

    rc_MAD = stats.median_abs_deviation(samples[:, 0])
    rt_MAD = np.inf
    if ndim != 1:
        rt_MAD = stats.median_abs_deviation(samples[:, 1])
    if ndim == 4:
        ell_MAD = stats.median_abs_deviation(samples[:, 2])

    # Central density. Use mean values for all the parameters.
    KP_cent_dens = lnlike(
        (rc_m, rt_m, ell_m, theta_m), ndim, rt_max, cl_cent, field_dens,
        N_memb, xy_in, r_in, rt_rang, True)

    # 16th-84th percentile region for the profile fit
    ell, theta = 0., 0.
    kpf_yvals, cd_rc_rt_ell_sampled = [], []

    # Sample rc, rt, ell, theta. Use the median and MAD.
    N = 1000
    rc_s = np.random.normal(rc_50, 1.4826 * rc_MAD, N)
    rt_s = np.random.normal(rt_50, 1.4826 * rt_MAD, N)
    if ndim == 4:
        ell_s = np.random.normal(ell_50, 1.4826 * ell_MAD, N)
        ell_s = np.clip(ell_s, a_min=0, a_max=.99)
        theta_s = np.random.normal(theta_m, theta_std, N)

    for _ in range(N):
        rc, rt, ell = rc_s[_], rt_s[_], np.nan
        if ndim == 4:
            ell, theta = ell_s[_], theta_s[_]

        # with warnings.catch_warnings():
        #     warnings.simplefilter("ignore")
        # Obtain central density
        KP_cd_ = lnlike(
            (rc, rt, ell, theta), ndim, rt_max, cl_cent, field_dens,
            N_memb, xy_in, r_in, rt_rang, True)
        cd_rc_rt_ell_sampled.append([KP_cd_, rc, rt, ell])

        # Values in y
        kpf_yvals.append(KP_cd_ * KingProf(rt_rang, rc, rt) + field_dens)

    cd_rc_rt_ell_sampled = np.array(cd_rc_rt_ell_sampled).T

    kpf_yvals = np.array(kpf_yvals).T
    _16_kp = np.nanpercentile(kpf_yvals, 16, 1)
    _84_kp = np.nanpercentile(kpf_yvals, 84, 1)

    KP_plot = {
        'KP_steps': KP_steps,
        'KP_mean_afs': KP_mean_afs, 'KP_tau_autocorr': KP_tau_autocorr,
        'KP_ESS': KP_ESS, 'KP_samples': KP_samples, 'KP_mode': KP_mode,
        'KP_kde': KP_kde, 'KP_cent_dens': KP_cent_dens,
        '_16_84_rang': rt_rang, '_16_kp': _16_kp, '_84_kp': _84_kp,
        'cd_rc_rt_ell_sampled': cd_rc_rt_ell_sampled}

    return KP_plot
