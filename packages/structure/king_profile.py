
import numpy as np
from scipy import spatial
import warnings
from .. import update_progress
from ..out import prep_plots
from ..aux_funcs import circFrac


def main(
    clp, cld_i, kp_flag, kp_nchains, kp_nruns, kp_nburn, rt_max_f, coords,
        **kwargs):
    """
    Fit a King profile given an array of stars. The field density value is
    fixed and the core radius, tidal radius and maximum central density are
    fitted.
    """

    # Check flag to run or skip.
    if kp_flag:
        print("Estimating King profile")
        KP_cd, KP_steps, KP_mean_afs, KP_tau_autocorr, KP_ESS, KP_samples,\
            KP_Bys_rc, KP_Bys_rt = fit_King_prof(
                kp_nchains, kp_nruns, kp_nburn, cld_i['x'], cld_i['y'],
                clp['kde_cent'], clp['field_dens'],
                clp['clust_rad'], clp['n_memb_i'], clp['N_MC'],
                clp['rand_01_MC'], clp['cos_t'], clp['sin_t'], rt_max_f)

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
        KP_Bys_rc, KP_Bys_rt = np.array([np.nan, np.nan, np.nan]),\
            np.array([np.nan, np.nan, np.nan])

    clp['KP_cent_dens'], clp['KP_steps'], clp['KP_mean_afs'],\
        clp['KP_tau_autocorr'], clp['KP_ESS'], clp['KP_samples'],\
        clp['KP_Bys_rc'], clp['KP_Bys_rt'], clp['KP_memb_num'],\
        clp['KP_conct_par'] = KP_cd, KP_steps, KP_mean_afs, KP_tau_autocorr,\
        KP_ESS, KP_samples, KP_Bys_rc, KP_Bys_rt, KP_memb_num, KP_conct_par
    return clp


def fit_King_prof(
    nchains, nruns, nburn, x, y, cl_cent, fd, cl_rad, n_memb_i, N_MC,
    rand_01_MC, cos_t, sin_t, rt_max_f, N_integ=1000, KP_fit_cent=False,
        N_conv=1000, tau_stable=0.05):
    """
    Fit central density, core radius and tidal radius, using a fixed
    value for the field density (fd).

    rt_max_f: factor that caps the maximum tidal radius, given the previously
    estimated cluster radius.

    HARDCODED:

    N_integ: number of values used in the tidal radius array.
    KP_fit_cent : fit center coords or used fixed values.
    N_conv
    tau_stable
    """

    from emcee import ensemble
    from emcee import moves

    # HARDCODED
    # mv = moves.StretchMove()
    # mv = moves.KDEMove()
    mv = [
        (moves.DESnookerMove(), 0.1),
        (moves.DEMove(), 0.9 * 0.9),
        (moves.DEMove(gamma0=1.0), 0.9 * 0.1),
    ]

    # HARDCODED
    # Steps to store Bayes params.
    KP_steps = int(nruns * .01)

    xy = np.array((x, y)).T
    area_tot = np.ptp(x) * np.ptp(y)
    cx, cy = cl_cent

    # Frame limits
    x0, x1 = min(xy.T[0]), max(xy.T[0])
    y0, y1 = min(xy.T[1]), max(xy.T[1])
    dx0, dx1 = abs(cl_cent[0] - x0), abs(cl_cent[0] - x1)
    dy0, dy1 = abs(cl_cent[1] - y0), abs(cl_cent[1] - y1)
    dxy = min(dx0, dx1, dy0, dy1)

    # Parameters used to estimate the fraction of cluster area inside the frame
    # if it spills outside of it.
    circarea_pars = (dxy, x0, x1, y0, y1, N_MC, rand_01_MC, cos_t, sin_t)

    # This assumes that the tidal radius can not be larger than 'rt_max' times
    # the estimated cluster radius.
    rt_max = rt_max_f * cl_rad
    # Tidal radius array. Used for integrating
    rt_rang = np.linspace(0., rt_max, int(rt_max_f * N_integ))

    if KP_fit_cent is True:
        # Dimensions: cent_x, cent_y, rc, rt
        ndim = 4
        # Initial guesses.
        cx_pos0 = np.random.uniform(
            cx - .9 * cl_rad, cx + .9 * cl_rad, nchains)
        cy_pos0 = np.random.uniform(
            cy - .9 * cl_rad, cy + .9 * cl_rad, nchains)
        rc_pos0 = np.random.uniform(.05 * rt_max, rt_max, nchains)
        rt_pos0 = np.random.uniform(rc_pos0, rt_max, nchains)
        pos0 = np.array([cx_pos0, cy_pos0, rc_pos0, rt_pos0]).T
        r_in = []

    else:
        # Dimensions: rc, rt
        ndim = 2
        # Initial guesses.
        rc_pos0 = np.random.uniform(.05 * rt_max, rt_max, nchains)
        rt_pos0 = np.random.uniform(rc_pos0, rt_max, nchains)
        pos0 = np.array([rc_pos0, rt_pos0]).T

        # Distances to the center of the cluster.
        xy_cent_dist = spatial.distance.cdist([cl_cent], xy)[0]
        msk = xy_cent_dist <= rt_max
        r_in = xy_cent_dist[msk]

    args = {
        'ndim': ndim, 'rt_max': rt_max, 'rt_rang': rt_rang, 'fd': fd,
        'N_memb': n_memb_i, 'r_in': r_in, 'cl_cent': cl_cent, 'cl_rad': cl_rad,
        'xy': xy, 'area_tot': area_tot, 'circarea_pars': circarea_pars}

    # emcee sampler
    sampler = ensemble.EnsembleSampler(
        nchains, ndim, lnprob, kwargs=args, moves=mv)

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
        nburn = int(i * nburn)
        samples = sampler.get_chain(discard=nburn, flat=True)

        if ndim == 4:
            cx, cy, rc, rt = np.mean(samples, 0)
            print(cx, cy)
            e_cx, e_cy, e_rc, e_rt = np.percentile(samples, (16, 84), 0).T
        elif ndim == 2:
            rc, rt = np.mean(samples, 0)
            e_rc_16, e_rc_84, e_rt_16, e_rt_84 = np.percentile(
                samples, (16, 84), 0).T.flatten()
            KP_Bys_rc = np.array([e_rc_16, rc, e_rc_84])
            KP_Bys_rt = np.array([e_rt_16, rt, e_rt_84])

        # Effective sample size
        KP_ESS = samples.shape[0] / np.mean(sampler.get_autocorr_time(tol=0))

        # For plotting, (nsteps, nchains, ndim)
        KP_samples = sampler.get_chain()

    # Central density, for plotting.
    if ndim == 2:
        N_memb = n_memb_i
    else:
        N_memb = rMmembN(
            fd, rt_max, rt, xy, area_tot, (cx, cy), circarea_pars)[-1]
    KP_cd = centDens(N_memb, rc, rt, rt_rang)

    return KP_cd, KP_steps, KP_mean_afs, KP_tau_autocorr, KP_ESS, KP_samples,\
        KP_Bys_rc, KP_Bys_rt


def lnprob(
    pars, ndim, rt_max, rt_rang, fd, N_memb, r_in, cl_cent, cl_rad, xy,
        area_tot, circarea_pars):
    """
    Logarithmic posterior
    """
    if ndim == 2:
        rc, rt = pars
        # Prior.
        if rt <= rc or rc <= 0. or rt > rt_max:
            return -np.inf
    if ndim == 4:
        cx, cy, rc, rt = pars
        # Prior.
        if cx < cl_cent[0] - .7 * cl_rad or cx > cl_cent[0] + .7 * cl_rad or\
            cy < cl_cent[1] - .7 * cl_rad or cy > cl_cent[1] + .7 * cl_rad or\
                rt <= rc or rc <= 0. or rt > rt_max:
            return -np.inf

    return lnlike(
        pars, ndim, rt_max, rt_rang, fd, N_memb, r_in, xy, area_tot,
        circarea_pars)


def lnlike(
    pars, ndim, rt_max, rt_rang, fd, N_memb, r_in, xy, area_tot,
        circarea_pars):
    """
    As defined in Pieres et al. (2016)
    """

    if ndim == 2:
        rc, rt = pars
    if ndim == 4:
        cx, cy, rc, rt = pars
        r_in, N_memb = rMmembN(
            fd, rt_max, rt, xy, area_tot, (cx, cy), circarea_pars)

    # Values outside the tidal radius contribute 'fd'.
    r_in = np.clip(r_in, a_min=0., a_max=rt)

    li = centDens(N_memb, rc, rt, rt_rang) * KingProf(r_in, rc, rt) + fd

    return np.log(li).sum()


def rMmembN(fd, rt_max, rt, xy, area_tot, cxy, circarea_pars):
    """
    The field density is kept fixed.
    """
    # Distances to the estimated center of the cluster.
    xy_cent_dist = spatial.distance.cdist([cxy], xy)[0]
    msk = xy_cent_dist <= rt_max
    r_in = xy_cent_dist[msk]

    # # TODO areas that spill outside of frame
    # area_tr_max = np.pi * rt_max**2
    # area_out_rt_max = area_tot - area_tr_max

    # # New field density
    # fd = np.sum(~msk) / area_out_rt_max

    area_cl = np.pi * rt**2

    # Handle areas that spill outside of frame
    fr_area = 1.
    dxy = circarea_pars[0]
    if rt > dxy:
        x0, x1, y0, y1 = circarea_pars[1:5]
        fr_area = circFrac((cxy), rt, x0, x1, y0, y1, *circarea_pars[5:])
    area_cl *= fr_area

    msk = xy_cent_dist <= rt
    N_memb = max(1, msk.sum() - fd * area_cl)

    return r_in, N_memb


def centDens(N_memb, rc, rt, arr):
    """
    Central density constant. Integrate up to rt.
    """
    i = np.searchsorted(arr, rt)
    return N_memb / np.trapz(
        2. * np.pi * arr[:i] * KingProf(arr[:i], rc, rt), arr[:i])


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
