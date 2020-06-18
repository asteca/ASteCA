
import numpy as np
from astropy.stats import circmean
from astropy import units as u
from scipy import spatial
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

    # Check flag to run or skip.
    if kp_flag:
        print("Estimating King profile")
        KP_cd, KP_steps, KP_mean_afs, KP_tau_autocorr, KP_ESS, KP_samples,\
            KP_Bys_rc, KP_Bys_rt, KP_Bys_ecc, KP_Bys_theta, KP_Bayes_kde =\
            fit_King_prof(
                kp_nchains, kp_nruns, kp_nburn, cld_i['x'], cld_i['y'],
                clp['kde_cent'], clp['field_dens'],
                clp['clust_rad'], clp['n_memb_i'], rt_max_f)

        if not np.isnan([KP_cd, KP_Bys_rc[1], KP_Bys_rt[1]]).any():
            # Obtain number of members and concentration parameter.
            KP_memb_num, KP_conct_par = num_memb_conc_param(
                KP_cd, KP_Bys_rc[1], KP_Bys_rt[1])

            # Print results.
            coord = prep_plots.coord_syst(coords)[0]
            # Set precision of printed values.
            text2 = '{:.1f}, {:.1f}' if coord == 'px' else '{:g}, {:g}'
            text = 'Core & tidal radii obtained: ' + text2 + ' {}'
            print(text.format(KP_Bys_rc[1], KP_Bys_rt[1], coord))
        else:
            KP_memb_num, KP_conct_par = np.nan, np.nan

    else:
        print("Skipping King profile fit")
        KP_cd, KP_steps, KP_mean_afs, KP_tau_autocorr, KP_ESS, KP_samples,\
            KP_memb_num, KP_conct_par = [np.nan] * 8
        KP_Bys_rc, KP_Bys_rt, KP_Bys_ecc, KP_Bys_theta =\
            [np.array([np.nan] * 5) for _ in range(4)]
        KP_Bayes_kde = np.array([])

    clp['KP_cent_dens'], clp['KP_steps'], clp['KP_mean_afs'],\
        clp['KP_tau_autocorr'], clp['KP_ESS'], clp['KP_samples'],\
        clp['KP_Bys_rc'], clp['KP_Bys_rt'], clp['KP_Bys_ecc'],\
        clp['KP_Bys_theta'], clp['KP_Bayes_kde'],\
        clp['KP_memb_num'], clp['KP_conct_par'] = KP_cd, KP_steps,\
        KP_mean_afs, KP_tau_autocorr, KP_ESS, KP_samples, KP_Bys_rc,\
        KP_Bys_rt, KP_Bys_ecc, KP_Bys_theta, KP_Bayes_kde, KP_memb_num,\
        KP_conct_par
    return clp


def fit_King_prof(
    nchains, nruns, nburn, x, y, cl_cent, field_dens, cl_rad, n_memb_i,
        rt_max_f, N_integ=1000, N_conv=1000, tau_stable=0.05):
    """
    rt_max_f: factor that caps the maximum tidal radius, given the previously
    estimated cluster radius.
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

    # Field density and estimated number of members (previously
    # obtained)
    fd = field_dens

    # Fix N_memb
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

    # Initial positions for the sampler.
    if ndim == 2:
        # Dimensions: rc, rt
        rc_pos0 = np.random.uniform(.05 * rt_max, rt_max, nchains)
        rt_pos0 = np.random.uniform(rc_pos0, rt_max, nchains)
        pos0 = np.array([rc_pos0, rt_pos0]).T
    elif ndim == 4:
        # Dimensions: rc, rt, ecc, theta
        rc_pos0 = np.random.uniform(.05 * rt_max, rt_max, nchains)
        rt_pos0 = np.random.uniform(rc_pos0, rt_max, nchains)
        ecc = np.random.uniform(0., 1., nchains)
        theta = np.random.uniform(0., np.pi, nchains)
        pos0 = np.array([rc_pos0, rt_pos0, ecc, theta]).T

    # Identify stars inside the cut-out given by the 'rt_max' value. Only these
    # stars will be processed below.
    xy = np.array((x, y)).T
    xy_cent_dist = spatial.distance.cdist([cl_cent], xy)[0]
    msk = xy_cent_dist <= rt_max
    xy_in = xy[msk].T
    r_in = xy_cent_dist[msk]

    args = {
        'ndim': ndim, 'rt_max': rt_max, 'cl_cent': cl_cent, 'fd': fd,
        'N_memb': N_memb, 'rt_rang': rt_rang, 'xy_in': xy_in, 'r_in': r_in}

    # emcee sampler
    sampler = ensemble.EnsembleSampler(
        nchains, ndim, lnprob, kwargs=args, moves=mv)

    # Run the smpler hiding some warnings
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
            tau = np.mean(sampler.get_autocorr_time(tol=0))
            autocorr_vals[tau_index] = tau
            tau_index += 1

            # Check convergence
            converged = tau * (N_conv / nchains) < i * nburn
            converged &= np.abs(old_tau - tau) / tau < tau_stable
            if converged:
                print("")
                break
            old_tau = tau

            update_progress.updt(nruns, i + 1)

        KP_mean_afs = afs[:tau_index]
        KP_tau_autocorr = autocorr_vals[:tau_index]

        # Remove burn-in
        nburn = int(i * nburn)
        samples = sampler.get_chain(discard=nburn, flat=True)

        # Extract mean, median, mode, 16th, 84th percentiles for each parameter
        rc, rt = np.mean(samples[:, :2], 0)
        rc_16, rc_50, rc_84, rt_16, rt_50, rt_84 = np.percentile(
            samples[:, :2], (16, 50, 84), 0).T.flatten()

        if ndim == 2:
            ecc, theta = 0., 0.
            ecc_16, ecc_50, ecc_84, theta_16, theta_50, theta_84 =\
                [np.array([np.nan] * 3) for _ in range(6)]
            # Mode and KDE to plot
            # This simulates the 'fundam_params and 'varIdxs' arrays.
            fp, vi = [[-np.inf, np.inf], [-np.inf, np.inf]], [0, 1]
            KP_Bys_mode, KP_Bayes_kde = modeKDE(fp, vi, samples.T)
            KP_Bys_mode += [0., 0.]
            KP_Bayes_kde += [[], []]

        elif ndim == 4:
            ecc = np.mean(samples[:, 2], 0)
            theta = circmean(samples[:, 3] * u.rad).value
            # Beware: the median and percentiles for theta might not be
            # properly defined.
            ecc_16, ecc_50, ecc_84, theta_16, theta_50, theta_84 =\
                np.percentile(samples[:, 2:], (16, 50, 84), 0).T.flatten()
            # Estimate the mode
            fp, vi = [
                [-np.inf, np.inf], [-np.inf, np.inf], [-np.inf, np.inf],
                [-np.inf, np.inf]], [0, 1, 2, 3]
            KP_Bys_mode, KP_Bayes_kde = modeKDE(fp, vi, samples.T)

        # Store: 16th, median, 84th, mean, mode
        KP_Bys_rc = np.array([rc_16, rc_50, rc_84, rc, KP_Bys_mode[0]])
        KP_Bys_rt = np.array([rt_16, rt_50, rt_84, rt, KP_Bys_mode[1]])
        KP_Bys_ecc = np.array([ecc_16, ecc_50, ecc_84, ecc, KP_Bys_mode[2]])

        KP_Bys_theta = np.array([
            theta_16, theta_50, theta_84, theta, KP_Bys_mode[3]])

        # Effective sample size
        KP_ESS = samples.shape[0] / np.mean(sampler.get_autocorr_time(tol=0))

        # For plotting, (nsteps, nchains, ndim)
        KP_samples = sampler.get_chain()

    # Central density, for plotting. Use mean values for all the parameters.
    KP_cd = lnlike(
        (rc, rt, ecc, theta), ndim, rt_max, cl_cent, fd, N_memb, xy_in, r_in,
        rt_rang, True)

    return KP_cd, KP_steps, KP_mean_afs, KP_tau_autocorr, KP_ESS, KP_samples,\
        KP_Bys_rc, KP_Bys_rt, KP_Bys_ecc, KP_Bys_theta, KP_Bayes_kde


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
    if rt <= rc or rc <= 0. or rt > rt_max or ecc < 0. or ecc > 1. or\
            theta < 0. or theta > np.pi:
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
        r_in = np.clip(r_in, a_min=0., a_max=rt)
        N_in_region = r_in.size
        # Evaluate each star in King's profile
        KP = KingProf(r_in, rc, rt)

    elif ndim == 4:
        x, y = xy_in
        # Identify stars inside this ellipse
        in_ellip_msk = inEllipse(xy_in, cl_cent, rt, ecc, theta)
        N_in_region = in_ellip_msk.sum()

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

    if N_memb is None:
        N_memb = NmembEst(ndim, fd, N_in_region, rt_max, rt, ecc)

    # Central density
    k = centDens(N_memb, rt_rang, rc, rt, ecc)
    if return_dens is True:
        return k

    # Likelihood
    li = k * KP + fd

    # Sum of log-likelighood
    sum_log_lkl = np.log(li).sum()

    return sum_log_lkl


def NmembEst(ndim, fd, N_in_region, rt_max, rt, ecc):
    """
    The number of true members is estimated as the total number of stars
    inside the region, minus the expected number of field stars in
    the region.
    """
    if ndim == 2:
        N_field_stars = fd * rt_max**2
    elif ndim == 4:
        a = rt
        b = a * np.sqrt(1. - ecc**2)
        N_field_stars = (np.pi * a * b) * fd

    N_memb = N_in_region - N_field_stars
    # Minimum value is 1
    N_memb = max(1, N_memb)
    return N_memb


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


#############################################################################
# All the code below is for testing purposes only

# def NmembEst(fd, r_in, rt_max):
#     """
#     """
#     N_memb = r_in.size - fd * rt_max**2

#     return N_memb


# def rMmembN(fd, rt_max, rt, xy, area_tot, cxy, circarea_pars):
#     """
#     The field density is kept fixed.
#     """
#     # Distances to the estimated center of the cluster.
#     xy_cent_dist = spatial.distance.cdist([cxy], xy)[0]
#     msk = xy_cent_dist <= rt_max
#     r_in = xy_cent_dist[msk]

#     # # TODO areas that spill outside of frame
#     # area_tr_max = np.pi * rt_max**2
#     # area_out_rt_max = area_tot - area_tr_max

#     # # New field density
#     # fd = np.sum(~msk) / area_out_rt_max

#     area_cl = np.pi * rt**2

#     # Handle areas that spill outside of frame
#     fr_area = 1.
#     dxy = circarea_pars[0]
#     if rt > dxy:
#         x0, x1, y0, y1 = circarea_pars[1:5]
#         fr_area = circFrac((cxy), rt, x0, x1, y0, y1, *circarea_pars[5:])
#     area_cl *= fr_area

#     msk = xy_cent_dist <= rt
#     N_memb = max(1, msk.sum() - fd * area_cl)

#     return r_in, N_memb


    # def temp(rc, rt, ecc, theta, sampler=None):
    #     import matplotlib.pyplot as plt
    #     from matplotlib import patches as mpatches

    #     loglkl = lnlike((rc, rt, ecc, theta), rt_rang, fd, N_memb, xy_in, cl_cent)
    #     print(loglkl)

    #     Nmemb_ellip, msk = NmembEllipse(xy_in, cl_cent, rt, ecc, theta, fd)
    #     print(N_memb, Nmemb_ellip)

    #     if sampler is not None:
    #         samples = sampler.get_chain(discard=nburn)
    #         plt.subplot(231)
    #         plt.plot(samples[:, :, 1])
    #         plt.subplot(232)
    #         plt.plot(samples[:, :, 2])
    #         plt.subplot(233)
    #         plt.plot(samples[:, :, 3])

    #         samples = sampler.get_chain(discard=nburn, flat=True)
    #         plt.subplot(234)
    #         plt.hist(samples.T[2], 50)
    #         plt.subplot(235)
    #         plt.hist(samples.T[3], 50)

    #         ax = plt.subplot(236)
    #     else:
    #         ax = plt.subplot(111)

    #     a2 = rt**2
    #     b2 = a2 * (1. - ecc**2)
    #     ellipse = mpatches.Ellipse(
    #         xy=cl_cent, width=2. * np.sqrt(a2), height=2. * np.sqrt(b2),
    #         angle=np.rad2deg(theta),
    #         facecolor='None', edgecolor='black', linewidth=2,
    #         transform=ax.transData)
    #     ax.add_patch(ellipse)
    #     plt.scatter(*xy_in[:, msk], c='g', s=3)
    #     plt.scatter(*xy_in[:, ~msk], c='r', s=3)
    #     plt.show()


# def test():
#     """
#     """
#     import matplotlib.pyplot as plt
#     from matplotlib import patches as mpatches
#     from scipy.optimize import differential_evolution

#     seed = np.random.randint(10000)
#     print(seed, "\n")
#     np.random.seed(seed)

#     Nf, Ncl = 2000, 200
#     cl_cent = (1000., 1000.)
#     width = 300.
#     ecc, theta = np.random.uniform(), np.random.uniform(0., 180.)
#     height = width * (1. - ecc**2)

#     N_in_ellip = 0.
#     xy_clust = []
#     while N_in_ellip < Ncl:
#         xy = np.random.uniform(
#             cl_cent[0] - width, cl_cent[1] + width, (2, 1000))
#         x, y = xy
#         dx, dy = x - cl_cent[0], y - cl_cent[1]
#         x1 = dx * np.cos(theta) + dy * np.sin(theta)
#         y1 = dx * np.sin(theta) - dy * np.cos(theta)
#         dist = (x1 / width)**2 + (y1 / height)**2
#         msk = dist <= 1.
#         xy_clust += list(xy.T[msk])
#         N_in_ellip = len(xy_clust)
#     xy_clust = np.array(xy_clust)[:Ncl].T

#     # rt = width
#     # rc = rt / 4.
#     # cl_dists = invTrnsfSmpl(rc, rt, Ncl)
#     # # Generate positions for cluster members, given heir KP distances to the
#     # # center.
#     # phi = np.random.uniform(0., 1., Ncl) * 2 * np.pi
#     # x_cl = cl_cent[0] + cl_dists * np.cos(phi)
#     # y_cl = cl_cent[1] + cl_dists * np.sin(phi)
#     # xy_clust = np.array([x_cl, y_cl])

#     frame_rang = 2000.
#     xy_field = np.random.uniform(0., frame_rang, (2, Nf))
#     xy = np.concatenate([xy_field.T, xy_clust.T])

#     rt_max = width * 2
#     N_integ = 1000
#     rt_rang = np.linspace(0., rt_max, N_integ)

#     fd = Nf / frame_rang**2

#     # Stars inside the cut-out given by the 'rt_max' value.
#     xy_cent_dist = spatial.distance.cdist([cl_cent], xy)[0]
#     msk = xy_cent_dist <= rt_max
#     # r_in = xy_cent_dist[msk]
#     xy_in = xy[msk].T

#     def temp(rc, rt, ecc=0., theta=0.):
#         # loglkl = lnlike(
#         #     (rc, rt, ecc, theta), rt_rang, fd, Ncl, xy_in, cl_cent)
#         # print(loglkl)

#         N_memb, msk = NmembEllipse(xy_in, cl_cent, rt, ecc, theta, fd)
#         # print("N_memb: {}".format(N_memb))

#         ax = plt.subplot(111)
#         a2 = rt**2
#         b2 = a2 * (1. - ecc**2)
#         ellipse = mpatches.Ellipse(
#             xy=cl_cent, width=2. * np.sqrt(a2), height=2. * np.sqrt(b2),
#             angle=np.rad2deg(theta),
#             facecolor='None', edgecolor='black', linewidth=2,
#             transform=ax.transData)
#         ax.add_patch(ellipse)
#         plt.scatter(*xy_field, c='k', s=3)
#         plt.scatter(*xy_clust, c='g', s=15, marker='s', zorder=5)
#         plt.scatter(*xy_in[:, msk], c='r', s=5, zorder=1)
#         plt.scatter(*xy_in[:, ~msk], c='r', s=5, zorder=1)
#         plt.axvline(1000.)
#         plt.axhline(1000.)
#         plt.show()

#     ndim = 4

#     bounds = ((.05 * rt_max, rt_max), (.05 * rt_max, rt_max),
#               (0., 1.), (0., np.pi))
#     result = differential_evolution(
#         lnprob, bounds,
#         args=(ndim, rt_rang, fd, Ncl, xy_in, cl_cent, rt_max),
#         maxiter=2000, popsize=25, tol=0.00001)
#     r_rc, r_rt, r_ecc, r_theta = result.x
#     print("{:.2f}, {:.2f}, {:.1f}".format(width, ecc, theta))
#     print("{:.2f}, {:.2f}, {:.1f}".format(r_rt, r_ecc, np.rad2deg(r_theta)))
#     print(result)
#     temp(*result.x)


# def invTrnsfSmpl(rc, rt, N_samp, N_interp=1000):
#     """
#     Sample King's profile using the inverse CDF method.
#     """
#     from scipy.integrate import quad
#     from scipy.interpolate import interp1d

#     def rKP(r, rc, rt):
#         return r * KingProf(r, rc, rt)

#     r_0rt = np.linspace(0., rt, N_interp)
#     # The CDF is defined as: $F(r)= \int_{r_low}^{r} PDF(r) dr$
#     # Sample the CDF
#     CDF_samples = []
#     for r in r_0rt:
#         CDF_samples.append(quad(rKP, 0., r, args=(rc, rt))[0])

#     # Normalize CDF
#     CDF_samples = np.array(CDF_samples) / CDF_samples[-1]

#     # Inverse CDF
#     inv_cdf = interp1d(CDF_samples, r_0rt)

#     # Sample the inverse CDF
#     samples = inv_cdf(np.random.rand(N_samp))

#     return samples


# if __name__ == '__main__':
#     test()
