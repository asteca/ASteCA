
import time as t
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt

from packages.inp import get_tracks
from packages.best_fit import obs_clust_prepare
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
2. It uses the default 'CLUSTER.dat' file in the 'defvals' folder
3. It reads Gaia EDR3 isochrones from the 'isochrones'
"""


def main(max_time=100, Mmin=100, Mmax=1000):
    """
    """
    np.random.seed(12345)

    lkl_method, obs_clust, completeness, err_lst,\
        max_mag_syn, fundam_params, ext_coefs, binar_flag, mean_bin_mr, N_fc,\
        m_ini_idx, st_dist_mass, theor_tracks, err_norm_rand, binar_probs,\
        ext_unif_rand = inParams(Mmin, Mmax)

    # Pack synthetic cluster arguments.
    synthcl_args = [
        fundam_params, completeness, err_lst, max_mag_syn, ext_coefs,
        binar_flag, mean_bin_mr, N_fc, m_ini_idx, st_dist_mass, theor_tracks,
        err_norm_rand, binar_probs, ext_unif_rand]

    print("Running")

    times_all = []
    elapsed, start, N_tot, models_sec, tstep = 0., t.time(), 0, [], 10
    while elapsed < max_time:

        model = []
        for p in fundam_params:
            model.append(np.random.uniform(min(p), max(p)))

        data = synthClust(model, lkl_method, obs_clust, synthcl_args)

        now = t.time()
        models_sec.append(1. / (now - start))

        N_tot += 1
        elapsed += now - start
        if elapsed >= tstep:
            print(("{:.1f} | {:.0f}").format(elapsed, models_sec[-1]))
            tstep += 10

        times_all.append(data)
        start = t.time()
        if elapsed >= max_time:
            break

    plot(Mmin, Mmax, lkl_method, fundam_params, elapsed, times_all,
         models_sec, N_tot)


def inParams(Mmin, Mmax):
    """
    This module is built to work with the default 'CLUSTER.dat' file and
    Gaia EDR3 Parsec+No isochrones
    """
    lkl_method, lkl_binning, lkl_manual_bins = 'tremmel', 'knuth', None

    zmin, zmax, amin, amax = 0.01, 0.02, 7, 9.5
    fundam_params_all = {'CLUSTER': [
        [zmin, zmax], [amin, amax], [0.0, 1.0], [8.0, 15.0],
        [Mmin, Mmax], [0.0, 0.1], [3.1]]}
    priors_mcee_all = {'CLUSTER': [['u'], ['u'], ['u'], ['u'], ['u'], ['u']]}

    #
    cldata = ascii.read('packages/defvals/input/CLUSTER.dat')
    cl_max_mag0 = list(zip(*[
        cldata['EDR3Name'], cldata['_x'], cldata['_x'], cldata['Gmag'],
        cldata['e_Gmag'], cldata['BP-RP'], cldata['e_BP-RP']]))
    cl_max_mag = []
    for st in cl_max_mag0:
        cl_max_mag.append(
            list(st[:3]) + [[st[3]]] + [[st[4]]] + [[st[5]]] + [[st[6]]]
            + list(np.array([np.ones(4) * np.nan]))
            + list(np.array([np.ones(4) * np.nan])) + [1])
    max_mag_syn = max(cldata['Gmag'])
    obs_clust = obs_clust_prepare.main(
        cl_max_mag, lkl_method, lkl_binning, lkl_manual_bins)

    completeness = [
        np.array([
            11.0578, 15.569595, 16.2036, 16.61846, 16.8903, 17.08015,
            17.24488, 17.407075, 17.56784, 17.7021, 17.81175, 17.924355,
            18.02364, 18.115365, 18.21082, 18.313775, 18.41914, 18.545955,
            18.67507, 18.80564, 18.9959]),
        np.array([
            1., 0.00716116, 0., 0., 0.07651589,
            0., 0., 0., 0.00714035, 0.02800014,
            0., 0.0454411, 0.00454196, 0., 0.05118521,
            0.33643902, 0.39760438, 0.53978769, 0.59136226, 0.6437164,
            0.78133205]), 0.0]
    err_lst = [
        np.array([9.97217805e-09, 6.68941895e-01, 2.74837659e-02]),
        np.array([4.56490040e-09, 6.75067801e-01, 1.37429298e-02])]

    pd = {'best_fit_algor': 'y',
          'fundam_params_all': fundam_params_all,
          'priors_mcee_all': priors_mcee_all,
          'all_syst_filters': [('gaiaedr3', 'G_RPmag', 'Gmag', 'G_BPmag')],
          'evol_track': 'PAR12+No', 'iso_paths': ['./isochrones/gaiaedr3'],
          'CMD_extra_pars': (
              'Mini', 'int_IMF', 'Mass', 'logL', 'logTe', 'logg', 'label',
              'mbolmag'), 'synth_rand_seed': None,
          'cmd_systs': {'gaiaedr3': (
              ('Gmag', 'G_BPmag', 'G_RPmag'),
              (6422.01, 5335.42, 7739.17))}, 'filters': [('gaiaedr3', 'Gmag')],
          'colors': [('gaiaedr3', 'G_BPmag,G_RPmag')],
          'IMF_name': 'kroupa_2002', 'N_interp': 500, 'min_bmass_ratio': .7,
          'lkl_method': lkl_method}
    td = get_tracks.main(pd, 'CLUSTER')

    return lkl_method, obs_clust, completeness, err_lst, max_mag_syn,\
        td['fundam_params'], td['ext_coefs'], td['binar_flag'],\
        td['mean_bin_mr'], td['N_fc'], td['m_ini_idx'],\
        td['st_dist_mass'], td['theor_tracks'], td['err_norm_rand'],\
        td['binar_probs'], td['ext_unif_rand']


def synthClust(model, lkl_method, obs_clust, synthcl_args):
    """
    """
    fundam_params, completeness, err_lst,\
        max_mag_syn, ext_coefs, binar_flag, mean_bin_mr, N_fc, m_ini_idx,\
        st_dist_mass, theor_tracks, err_norm_rand, binar_probs, ext_unif_rand,\
        = synthcl_args

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
    e, d, M_total, bin_frac, R_V = model_proper

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
                    isoch_compl, err_lst, err_norm_rand[ml])
                t8 = t.time() - s

    s = t.time()
    # Transposing is necessary for np.histogramdd() in the likelihood
    synth_clust = synth_clust[:sum(N_fc)].T
    _ = likelihood.main(lkl_method, synth_clust, obs_clust)
    t9 = t.time() - s

    return t0, t1, t2, t3, t4, t5, t6, t7, t8, t9


def plot(Mmin, Mmax, lkl_method, fundam_params, elapsed, times_all,
         models_sec, N_tot):
    """
    """
    times_norm = 100. * np.sum(times_all, 0) / elapsed
    cols = [
        'propModel', 'zaWAvrg', 'move', 'cut', 'M_dst', 'M_interp', 'binar',
        'complete', 'errors', lkl_method]

    for i, c in enumerate(cols):
        print("{}: {}".format(c, times_norm[i]))

    # fig, axs = plt.subplots(2, 1, figsize=(20, 20))
    fig, axs = plt.subplots(1, 1, figsize=(20, 10))

    # M_min, M_max = fundam_params[-2]
    # axs[0].set_title(
    #     "Exclude init 5% of runs | M=[{:.0f}, {:.0f}]".format(M_min, M_max),
    #     fontsize=18)
    # axs[0].plot(models_sec[-int(.95 * len(models_sec)):])
    # axs[0].grid(True, zorder=0)
    # axs[0].xaxis.set_tick_params(labelsize=18)
    # axs[0].yaxis.set_tick_params(labelsize=18)

    _16, _50, _84 = np.percentile(models_sec, (16, 50, 84))
    axs.set_title(
        ("M=[{:.0f}, {:.0f}] | N={}, t={:.0f} --> [{:.0f} +/-({:.0f}, "
         + "{:.0f}) m/s]").format(
            Mmin, Mmax, N_tot, elapsed, _50, _16, _84), fontsize=18)
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


if __name__ == '__main__':
    main()
