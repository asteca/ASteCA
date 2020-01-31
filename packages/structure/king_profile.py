
import numpy as np
from scipy import spatial
import warnings
from .. import update_progress
from ..out import prep_plots


def main(
    clp, cld_i, inst_packgs_lst, kp_flag, kp_nchains, kp_nruns, rt_max_f,
        coords, **kwargs):
    """
    Fit a King profile given an array of stars. The field density value is
    fixed and the core radius, tidal radius and maximum central density are
    fitted.
    """

    # Initial dummy values, used if function is skipped.
    cd, rc, e_rc, rt, e_rt, mean_afs_KP, tau_autocorr_KP, N_ess_KP, n_c_k,\
        kcp = [np.nan] * 10

    # Check flag to run or skip.
    if kp_flag:
        if 'emcee' not in inst_packgs_lst:
            print("  WARNING: 'emcee' is not installed. Can not fit" +
                  " King's profile")

        else:
            print("Estimating King profile")
            cd, rc, e_rc, rt, e_rt, mean_afs_KP, tau_autocorr_KP, N_ess_KP =\
                fit_3P_King_prof(
                    cld_i['x'], cld_i['y'], clp['kde_cent'], clp['field_dens'],
                    clp['clust_rad'], clp['n_memb_i'], kp_nchains, kp_nruns,
                    rt_max_f)

            if not np.isnan([cd, rt, rc]).any():
                # Obtain number of members and concentration parameter.
                n_c_k, kcp = num_memb_conc_param(cd, rt, rc)

                # Print results.
                coord = prep_plots.coord_syst(coords)[0]
                # Set precision of printed values.
                text2 = '{:.1f}, {:.1f}' if coord == 'px' else '{:g}, {:g}'
                text = 'Core & tidal radii obtained: ' + text2 + ' {}'
                print(text.format(rc, rt, coord))
            else:
                n_c_k, kcp = np.nan, np.nan

    clp['KP_cent_dens'], clp['core_rad'], clp['e_core'], clp['tidal_rad'],\
        clp['e_tidal'], clp['mean_afs_KP'], clp['tau_autocorr_KP'],\
        clp['N_ess_KP'], clp['KP_memb_num'], clp['KP_conct_par'] =\
        cd, rc, e_rc, rt, e_rt, mean_afs_KP, tau_autocorr_KP, N_ess_KP,\
        n_c_k, kcp
    return clp


def fit_3P_King_prof(
    x, y, cl_cent, fd, cl_rad, n_memb_i, nchains, nruns, rt_max_f,
        N_integ=1000):
    """
    Fit central density, core radius and tidal radius, using a fixed
    value for the field density (fd).

    rt_max_f: factor that caps the maximum tidal radius, given the previously
    estimated cluster radius.

    HARDCODED:

    N_integ: number of values used in the tidal radius array.
    """
    from emcee import ensemble

    # Distances to the center of the cluster.
    xy_cent_dist = spatial.distance.cdist([cl_cent], np.array((x, y)).T)[0]

    # This assumes that the tidal radius can not be larger than 'rt_max' times
    # the estimated cluster radius.
    rt_max = rt_max_f * cl_rad
    msk = xy_cent_dist <= rt_max
    r_in = xy_cent_dist[msk]
    # Tidal radius array. Used for integrating
    rt_rang = np.linspace(0., rt_max, int(rt_max_f * N_integ))

    # Two dimensions: rt, rc
    ndim = 2
    # Initial guesses.
    rc_pos0 = np.random.uniform(.05 * rt_max, rt_max, nchains)
    rt_pos0 = np.random.uniform(rc_pos0, rt_max, nchains)
    pos0 = np.array([rc_pos0, rt_pos0]).T

    # emcee sampler
    sampler = ensemble.EnsembleSampler(
        nchains, ndim, lnprob, args=(rt_max, rt_rang, r_in, n_memb_i, fd))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        st_check = int(nruns * .1)
        tau_index, autocorr_vals, afs = 0, np.empty(nruns), np.empty(nruns)
        for i, (pos, prob, stat) in enumerate(
                sampler.sample(pos0, iterations=nruns)):

            # Every X steps
            if i % st_check and i < (nruns - 1):
                continue

            afs[tau_index] = np.mean(sampler.acceptance_fraction)
            tau = sampler.get_autocorr_time(tol=0)
            autocorr_vals[tau_index] = np.mean(tau)
            tau_index += 1

            update_progress.updt(nruns, i + 1)

        mean_afs = afs[:tau_index]
        tau_autocorr = autocorr_vals[:tau_index]
        # Remove burn-in (25% of chain)
        nburn = int(i * .25)
        samples = sampler.get_chain(discard=nburn, flat=True)

        rc, rt = np.mean(samples, 0)
        e_rc, e_rt = np.percentile(samples, (16, 84), 0)
        tau_mean = np.mean(sampler.get_autocorr_time(tol=0))
        # Effective sample size
        N_ess = samples.shape[0] / tau_mean

    cd = centDens(n_memb_i, (rc, rt), rt_rang)

    return cd, rc, e_rc, rt, e_rt, mean_afs, tau_autocorr, N_ess


def lnprob(rct, rt_max, rt_rang, r_in, N_memb, fd):
    """
    Logarithmic posterior
    """
    # Prior.
    if rct[1] <= rct[0] or rct[0] <= 0. or rct[1] > rt_max:
        return -np.inf

    return lnlike(rct, rt_rang, r_in, N_memb, fd)


def lnlike(rct, rt_rang, r_in, N_memb, fd):
    """
    As defined in Pieres et al. (2016)
    """

    # Values outside the tidal radius contribute 'fd'.
    r_in = np.clip(r_in, a_min=0., a_max=rct[1])

    li = centDens(N_memb, rct, rt_rang) * KingProf(rct, r_in) + fd

    return np.log(li).sum()


def centDens(N_memb, rct, arr):
    """
    Central density constant. Integrate up to rt (=rct[1])
    """
    i = np.searchsorted(arr, rct[1])
    return N_memb / np.trapz(
        2. * np.pi * arr[:i] * KingProf(rct, arr[:i]), arr[:i])


def KingProf(rct, r_in):
    """
    King (1962) profile.
    """
    return ((1. / np.sqrt(1. + (r_in / rct[0]) ** 2)) -
            (1. / np.sqrt(1. + (rct[1] / rct[0]) ** 2))) ** 2


def num_memb_conc_param(cd, rt, rc):
    """
    If 3-P King profile converged, ie: the tidal radius was found,
    calculate approximate number of cluster members with Eq (3) from
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
