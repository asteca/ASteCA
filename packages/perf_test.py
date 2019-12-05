
import time as t
import pickle
import sys
import numpy as np
import matplotlib.pyplot as plt

from packages.best_fit.bf_common import varPars
from packages.synth_clust import zaWAverage
from packages.synth_clust.synth_cluster import properModel
from packages.synth_clust import move_isochrone
from packages.synth_clust import cut_max_mag
from packages.synth_clust import mass_distribution
from packages.synth_clust import mass_interp
from packages.synth_clust import binarity
from packages.synth_clust import completeness_rm
from packages.synth_clust import add_errors
from packages.best_fit import likelihood


def main(model, lkl_method, obs_clust, fundam_params, synthcl_args, sel):
    """
    1. This file needs to be in the top folder
    2. The 'perf_test.pickle' is created by 'best_fit_synth_cl'
    3. The input 'perf_test.pickle' should also be present in the top folder
    """

    varIdxs, ndim, ranges = varPars(fundam_params)

    theor_tracks, e_max, err_lst, completeness, max_mag_syn, st_dist_mass,\
        R_V, ext_coefs, N_fc, m_ini, cmpl_rnd, err_rnd = synthcl_args

    s = t.time()
    model_proper, z_model, a_model, ml, mh, al, ah = properModel(
        fundam_params, model, varIdxs)
    t0 = t.time() - s

    s = t.time()
    isochrone = zaWAverage.main(
        theor_tracks, fundam_params, z_model, a_model, ml, mh, al, ah)
    t1 = t.time() - s

    # Generate synthetic cluster.
    e, d, M_total, bin_frac = model_proper
    s = t.time()
    isoch_moved = move_isochrone.main(isochrone, e, d, R_V, ext_coefs, N_fc)
    t2 = t.time() - s
    s = t.time()
    isoch_cut = cut_max_mag.main(isoch_moved, max_mag_syn)
    t3 = t.time() - s

    synth_clust = []
    if isoch_cut.any():
        s = t.time()
        mass_dist = mass_distribution.main(st_dist_mass, M_total)
        t4 = t.time() - s
        s = t.time()
        isoch_mass = mass_interp.main(isoch_cut, mass_dist, m_ini)
        t5 = t.time() - s
        if isoch_mass.any():
            s = t.time()
            isoch_binar = binarity.main(isoch_mass, bin_frac, m_ini, N_fc)
            t6 = t.time() - s
            s = t.time()
            isoch_compl = completeness_rm.main(
                isoch_binar, completeness, cmpl_rnd)
            t7 = t.time() - s
            if isoch_compl.any():
                s = t.time()
                synth_clust = add_errors.main(
                    isoch_compl, err_lst, e_max, m_ini, err_rnd)
                t8 = t.time() - s

    s = t.time()
    likelihood.main(lkl_method, synth_clust, obs_clust)
    t9 = t.time() - s

    return t0, t1, t2, t3, t4, t5, t6, t7, t8, t9


if __name__ == '__main__':

    # from synth_clust import imf
    mpath = sys.path[0].replace('synth_clust', '').replace('packages/', '')

    print("Reading data")
    with open(mpath + '/perf_test.pickle', 'rb') as f:
        obs_clust, fundam_params, theor_tracks, lkl_method, R_V, e_max,\
            err_lst, completeness, max_mag_syn, st_dist_mass, ext_coefs, N_fc,\
            m_ini, cmpl_rnd, err_rnd = pickle.load(f)

    # Pack synthetic cluster arguments.
    synthcl_args = [
        theor_tracks, e_max, err_lst, completeness, max_mag_syn, st_dist_mass,
        R_V, ext_coefs, N_fc, m_ini, cmpl_rnd, err_rnd]

    print("Running")
    np.random.seed(12345)

    t_new, t_old = 0., 0.
    times_all_new, times_all_old = [], []

    max_time, elapsed, start, N_tot = 600., 0., t.time(), 0
    while elapsed < max_time:

        model = []
        for p in fundam_params:
            model.append(np.random.uniform(min(p), max(p)))

        if np.random.randint(1000) == 500:
            print(("{}, {:.1f}| {:.4f}, {:.3f}, {:.2f}, {:.2f}, {:.0f}, "
                   "{:.1f}").format(N_tot, elapsed, *model))

        s = t.time()
        sel = 'new'
        times_all_new.append(main(
            model, lkl_method, obs_clust, fundam_params, synthcl_args,
            sel))
        t_new += t.time() - s

        s = t.time()
        sel = 'old'
        times_all_old.append(main(
            model, lkl_method, obs_clust, fundam_params, synthcl_args,
            sel))
        t_old += t.time() - s

        N_tot += 1
        elapsed += t.time() - start
        start = t.time()
        if elapsed >= max_time:
            break

    times_norm_old = 100. * np.sum(times_all_old, 0) / t_old
    times_norm_new = 100. * np.sum(times_all_new, 0) / t_new
    cols = [
        'propModel', 'zaWAvrg', 'move', 'cut', 'M_dst', 'M_intrp', 'binar',
        'complete', 'errors', lkl_method]

    for i, c in enumerate(cols):
        print("{}: {}, {}".format(c, times_norm_old[i], times_norm_new[i]))

    fig = plt.figure(figsize=(20, 10))
    plt.title("N={}, t={:.2f}".format(N_tot, t_old), fontsize=18)
    plt.grid(zorder=0)
    plt.bar(cols, times_norm_old, zorder=4)
    # Text annotations
    x, y = np.arange(len(times_norm_old)), np.round(times_norm_old, 2)
    up = max(y) * .03
    plt.ylim(0, max(y) + 4 * up)
    for xi, yi, l in zip(*[x, y, list(map(str, y))]):
        plt.text(xi - len(l) * .02, yi + up, l, fontsize=18,
                 bbox=dict(facecolor='w', edgecolor='w', alpha=.5))
    plt.ylabel("% of time used", fontsize=18)
    plt.xticks(rotation=45, fontsize=18)
    plt.yticks(fontsize=18)
    fig.tight_layout()
    plt.savefig("perf_test_old.png", dpi=150)

    fig = plt.figure(figsize=(20, 10))
    plt.title("N={}, t={:.2f}".format(N_tot, t_new), fontsize=18)
    plt.grid(zorder=0)
    plt.bar(cols, times_norm_new, zorder=4)
    # Text annotations
    x, y = np.arange(len(times_norm_new)), np.round(times_norm_new, 1)
    up = max(y) * .03
    plt.ylim(0, max(y) + 4 * up)
    for xi, yi, l in zip(*[x, y, list(map(str, y))]):
        plt.text(xi - len(l) * .02, yi + up, l, fontsize=18,
                 bbox=dict(facecolor='w', edgecolor='w', alpha=.5))
    plt.ylabel("% of time used", fontsize=18)
    plt.xticks(rotation=45, fontsize=18)
    plt.yticks(fontsize=18)
    fig.tight_layout()
    plt.savefig("perf_test_new.png", dpi=150)
