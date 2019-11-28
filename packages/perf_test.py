
import time as t
import pickle
import sys
import numpy as np
import matplotlib.pyplot as plt

from packages.best_fit.bf_common import varPars
from packages.synth_clust import zaWAverage
from packages.synth_clust import move_isochrone
from packages.synth_clust import cut_max_mag
from packages.synth_clust import mass_distribution
from packages.synth_clust import mass_interp
from packages.synth_clust import binarity
from packages.synth_clust import completeness_rm
from packages.synth_clust import add_errors
from packages.best_fit import likelihood


def main(model, lkl_method, obs_clust, fundam_params, synthcl_args):
    """
    1. This file needs to be in the top folder
    2. The 'perf_test.pickle' is created by 'best_fit_synth_cl'
    3. The input 'perf_test.pickle' should also be present in the top folder
    """

    varIdxs, ndim, ranges = varPars(fundam_params)

    theor_tracks, e_max, err_lst, completeness, max_mag_syn, st_dist_mass,\
        R_V, ext_coefs, N_fc, m_ini, cmpl_rnd, err_rnd = synthcl_args

    # Average a new isochrone
    s = t.time()
    isochrone, model_proper = zaWAverage.main(
        theor_tracks, fundam_params, varIdxs, model)
    t0 = t.time() - s

    # Generate synthetic cluster.
    # synth_clust = synth_cluster.main(
    #     isochrone, model_proper, *synthcl_args[1:])
    e, d, M_total, bin_frac = model_proper
    s = t.time()
    isoch_moved = move_isochrone.main(isochrone, e, d, R_V, ext_coefs, N_fc)
    t1 = t.time() - s
    s = t.time()
    isoch_cut = cut_max_mag.main(isoch_moved, max_mag_syn)
    t2 = t.time() - s
    synth_clust = []
    if isoch_cut.any():
        s = t.time()
        mass_dist = mass_distribution.main(st_dist_mass, M_total)
        t3 = t.time() - s
        s = t.time()
        isoch_mass = mass_interp.main(isoch_cut, mass_dist, m_ini)
        t4 = t.time() - s
        if isoch_mass.any():
            s = t.time()
            isoch_binar = binarity.main(isoch_mass, bin_frac, m_ini, N_fc)
            t5 = t.time() - s
            s = t.time()
            isoch_compl = completeness_rm.main(
                isoch_binar, completeness, cmpl_rnd)
            t6 = t.time() - s
            if isoch_compl.any():
                s = t.time()
                synth_clust = add_errors.main(
                    isoch_compl, err_lst, e_max, m_ini, err_rnd)
                t7 = t.time() - s

    s = t.time()
    likelihood.main(lkl_method, synth_clust, obs_clust)
    t8 = t.time() - s

    return t0, t1, t2, t3, t4, t5, t6, t7, t8


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

    Nruns = 10000
    print("Running")
    times_all = []

    s = t.time()
    for i in range(Nruns):
        if (i + 1) % 100:
            continue

        print(i)
        model = []
        for p in fundam_params:
            model.append(np.random.uniform(min(p), max(p)))

        times_all.append(main(
            model, lkl_method, obs_clust, fundam_params, synthcl_args))

    t_total = t.time() - s
    times_norm = 100. * np.sum(times_all, 0) / t_total
    cols = [
        'zaWAvrg', 'move', 'cut', 'M_dst', 'M_intrp', 'binar',
        'complete', 'errors', lkl_method]
    plt.title("N={}".format(Nruns))
    plt.grid(zorder=0)
    plt.bar(cols, times_norm, zorder=4)
    plt.ylabel("% of time used")
    plt.xticks(rotation=45)
    plt.show()
