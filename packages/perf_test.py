
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

"""

1. **This file needs to be in the top folder of the repo**

2. The 'perf_test.pickle' is created by 'best_fit_synth_cl'

3. The input 'perf_test.pickle' should also be present in the top folder

"""


def main(model, lkl_method, obs_clust, synthcl_args):
    """
    """
    fundam_params, completeness, err_lst, em_float,\
        max_mag_syn, ext_coefs, binar_flag, mean_bin_mr, N_fc, m_ini_idx,\
        st_dist_mass, theor_tracks, err_norm_rand, binar_probs, ext_unif_rand,\
        R_V = synthcl_args

    varIdxs, ndim, ranges = varPars(fundam_params)
    model = np.array(model)[varIdxs]

    t0, t1, t2, t3, t4, t5, t6, t7, t8, t9 = [0.] * 10

    s = t.time()
    model_proper, z_model, a_model, ml, mh, al, ah = properModel(
        fundam_params, model, varIdxs)
    t0 = t.time() - s

    s = t.time()
    isochrone = zaWAverage.main(
        theor_tracks, fundam_params, binar_flag, m_ini_idx, z_model, a_model,
        ml, mh, al, ah)
    t1 = t.time() - s

    # Generate synthetic cluster.
    e, d, M_total, bin_frac = model_proper

    s = t.time()
    isoch_moved = move_isochrone.main(
        isochrone, e, d, R_V, ext_coefs, N_fc, ext_unif_rand[ml], m_ini_idx,
        binar_flag)
    t2 = t.time() - s

    s = t.time()
    isoch_cut = cut_max_mag.main(isoch_moved, max_mag_syn)
    t3 = t.time() - s

    synth_clust = np.array([])
    if isoch_cut.any():

        s = t.time()
        mass_dist = mass_distribution.main(
            st_dist_mass[ml], mean_bin_mr, bin_frac, M_total)
        t4 = t.time() - s

        s = t.time()
        isoch_mass = mass_interp.main(isoch_cut, m_ini_idx, mass_dist)
        t5 = t.time() - s

        if isoch_mass.any():
            s = t.time()
            isoch_binar = binarity.main(
                isoch_mass, bin_frac, m_ini_idx, N_fc, binar_probs[ml])
            t6 = t.time() - s

            s = t.time()
            isoch_compl = completeness_rm.main(isoch_binar, completeness)
            t7 = t.time() - s

            if isoch_compl.any():
                s = t.time()
                synth_clust = add_errors.main(
                    isoch_compl, err_lst, em_float, err_norm_rand[ml])
                t8 = t.time() - s

    s = t.time()
    # Transposing is necessary for np.histogramdd() in the likelihood
    synth_clust = synth_clust[:sum(N_fc)].T
    _ = likelihood.main(lkl_method, synth_clust, obs_clust)
    t9 = t.time() - s

    return t0, t1, t2, t3, t4, t5, t6, t7, t8, t9


if __name__ == '__main__':

    # from synth_clust import imf
    mpath = sys.path[0].replace('synth_clust', '').replace('packages/', '')

    print("Reading data")
    with open(mpath + '/perf_test.pickle', 'rb') as f:
        fundam_params, completeness, err_lst, em_float, max_mag_syn,\
            ext_coefs, binar_flag, mean_bin_mr, N_fc, m_ini_idx,\
            st_dist_mass, theor_tracks, err_norm_rand, binar_probs,\
            ext_unif_rand, R_V, lkl_method, obs_clust = pickle.load(f)

    # Pack synthetic cluster arguments.
    synthcl_args = [
        fundam_params, completeness, err_lst, em_float, max_mag_syn, ext_coefs,
        binar_flag, mean_bin_mr, N_fc, m_ini_idx, st_dist_mass, theor_tracks,
        err_norm_rand, binar_probs, ext_unif_rand, R_V]

    print("Running")
    np.random.seed(12345)

    times_all = []
    max_time, elapsed, start, N_tot, models_sec = 600., 0., t.time(), 0, []
    while elapsed < max_time:

        model = []
        for p in fundam_params:
            model.append(np.random.uniform(min(p), max(p)))

        data = main(model, lkl_method, obs_clust, synthcl_args)
        times_all.append(data)

        N_tot += 1
        elapsed += t.time() - start
        models_sec.append(N_tot / elapsed)
        if np.random.randint(1000) == 500:
            print(("{}, {:.1f} | {:.0f}").format(
                N_tot, elapsed, models_sec[-1]))

        start = t.time()
        if elapsed >= max_time:
            break

    times_norm = 100. * np.sum(times_all, 0) / elapsed
    cols = [
        'propModel', 'zaWAvrg', 'move', 'cut', 'M_dst', 'M_interp', 'binar',
        'complete', 'errors', lkl_method]

    for i, c in enumerate(cols):
        print("{}: {}".format(c, times_norm[i]))

    fig, axs = plt.subplots(2, 1, figsize=(20, 20))

    M_min, M_max = fundam_params[-2]
    axs[0].set_title(
        "Exclude initial 5% of runs | M=[{:.0f}, {:.0f}]".format(M_min, M_max),
        fontsize=18)
    axs[0].plot(models_sec[-int(.95 * len(models_sec)):])
    axs[0].grid(True, zorder=0)
    axs[0].xaxis.set_tick_params(labelsize=18)
    axs[0].yaxis.set_tick_params(labelsize=18)

    axs[1].set_title("N={}, t={:.2f} -> [{:.0f} m/s]".format(
        N_tot, elapsed, N_tot / elapsed), fontsize=18)
    plt.grid(zorder=0)
    plt.bar(cols, times_norm, zorder=4)
    # Text annotations
    x, y = np.arange(len(times_norm)), np.round(times_norm, 2)
    up = max(y) * .03
    plt.ylim(0, max(y) + 4 * up)
    for xi, yi, l in zip(*[x, y, list(map(str, y))]):
        plt.text(xi - len(l) * .02, yi + up, l, fontsize=18,
                 bbox=dict(facecolor='w', edgecolor='w', alpha=.5))
    plt.ylabel("% of time used", fontsize=18)
    plt.xticks(rotation=45, fontsize=18)
    plt.yticks(fontsize=18)

    fig.tight_layout()
    plt.savefig("perf_test.png", dpi=150)
