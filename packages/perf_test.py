
import time as t
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt

from pyinstrument import Profiler
import memray

from packages.inp import read_met_files
from packages.best_fit import prep_obs_params
from packages.best_fit import prep_synth_params
from packages.best_fit.bf_common import getSynthClust
from packages.best_fit.bf_common import varPars
from packages.best_fit import likelihood

"""
1. **This file needs to be in the top folder of the repo**
2. It uses the default 'CLUSTER.dat' file in the 'defvals' folder
3. It reads Gaia EDR3 isochrones from the 'isochrones' folder
"""

# Define the range used for metallicity and age (isochrones must contain these
# ranges)
zmin, zmax, amin, amax = 0.01, 0.02, 7.5, 8.5
emin, emax, drmin, drmax = 0, 2, 0, .5
dmmin, dmmax, bmin, bmax = 5, 20, 0., .5
# Mass range
Max_mass = 1000
# Maximum running time (in seconds)
max_time = 20

np.random.seed(12345)


def profilerCall():
    """
    """
    # This makes 'mass_distribution' appear at the expense of taking a little
    # longer to process the results.
    profiler = Profiler(interval=0.0001)
    # profiler = Profiler()

    pd, clp, td = inParams(Max_mass)

    print("Running")
    profiler.start()
    Nmodels = main(pd, clp, td)
    profiler.stop()

    print("Extracting times")
    times_all = extractTimes(profiler)
    print("Plotting")
    plot(Max_mass, Nmodels, times_all)

    print("Preparing profiler report...")
    profiler.open_in_browser()


def memrayCall():
    """
    Usage:

    memray summary output_file.bin
    memray flamegraph output_file.bin
    """
    pd, clp, td = inParams(Max_mass)
    with memray.Tracker("output_file.bin"):
        main(pd, clp, td)
    print("Finished")


def main(pd, clp, td):
    """
    """
    model_proper = np.array([_[0] for _ in td['fundam_params']])
    # Pack common args for the 'synth_cluster.py' function.
    syntClustArgs = (
        pd['DR_dist'], pd['alpha'], model_proper,
        clp['varIdxs'], clp['completeness'], clp['err_lst'],
        clp['max_mag_syn'], clp['N_obs_stars'], td['fundam_params'],
        td['ext_coefs'], td['N_fc'], td['m_ini_idx'], td['st_dist_mass'],
        td['theor_tracks'], td['rand_norm_vals'], td['rand_unif_vals'])

    elapsed, start, Nmodels, tstep = 0., t.time(), 0, 10
    while elapsed < max_time:

        model = []
        for p in td['fundam_params']:
            model.append(np.random.uniform(min(p), max(p)))
        model = np.array(model)[clp['varIdxs']]

        # Call synthetic cluster generation
        synth_clust = getSynthClust(model, True, syntClustArgs)[0]
        # Call likelihood
        _ = likelihood.main(pd['lkl_method'], synth_clust, clp['obs_clust'])

        Nmodels += 1
        elapsed += t.time() - start
        if elapsed >= tstep:
            print(("{:.1f} | {:.0f}").format(elapsed, Nmodels / elapsed))
            tstep += 10

        start = t.time()
        if elapsed >= max_time:
            break

    return Nmodels


def inParams(Max_mass):
    """
    This module is built to work with the default 'CLUSTER.dat' file and
    Gaia EDR3 Parsec+No isochrones
    """

    fundam_params_all = {'CLUSTER': [
        [zmin, zmax], [amin, amax], [emin, emax], [drmin, drmax],
        [dmmin, dmmax], [bmin, bmax], [3.1]]}
    priors_mcee_all = {'CLUSTER': [['u'], ['u'], ['u'], ['u'], ['u'], ['u']]}

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
          'IMF_name': 'kroupa_2002', 'Max_mass': Max_mass, 'DR_dist': 'normal',
          'alpha': 0.275, 'gamma': 'D&K', 'lkl_method': 'tremmel',
          'lkl_binning': 'knuth',
          'lkl_manual_bins': None}

    # Read "observed" cluster
    cldata = ascii.read('packages/defvals/input/CLUSTER.dat')
    cl_reg_fit0 = list(zip(*[
        cldata['EDR3Name'], cldata['_x'], cldata['_x'], cldata['Gmag'],
        cldata['e_Gmag'], cldata['BP-RP'], cldata['e_BP-RP']]))
    cl_reg_fit = []
    for st in cl_reg_fit0:
        cl_reg_fit.append(
            list(st[:3]) + [[st[3]]] + [[st[4]]] + [[st[5]]] + [[st[6]]]
            + list(np.array([np.ones(4) * np.nan]))
            + list(np.array([np.ones(4) * np.nan])) + [1])
    # max_mag_syn = max(cldata['Gmag'])
    clp = {'cl_reg_fit': cl_reg_fit}

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
    clp['completeness'] = [completeness[0], completeness[1], False]

    err_lst = [
        np.array([9.97217805e-09, 6.68941895e-01, 2.74837659e-02]),
        np.array([4.56490040e-09, 6.75067801e-01, 1.37429298e-02])]
    clp['err_lst'] = err_lst

    clp['varIdxs'], clp['ndim'], clp['ranges'] = varPars(
        fundam_params_all)

    td = read_met_files.main(pd, 'CLUSTER')
    clp = prep_obs_params.main(clp, **pd)
    td = prep_synth_params.main(pd, clp, td)

    return pd, clp, td


def extractTimes(profiler):
    """
    """

    def getTime(line):
        l_s = line.split()
        try:
            i = l_s.index('`-')
        except ValueError:
            i = l_s.index('|-')
        return float(l_s[i + 1])

    t0, t1, t2, t3, t4, t5, t6, t7, t8, t9 = [0] * 10
    for line in profiler.output_text().split('\n'):
        if 'properModel' in line:
            t0 = getTime(line)
        elif 'main' in line and 'zaWAverage' in line:
            t1 = getTime(line)
        elif 'main' in line and 'move_isochrone' in line:
            t2 = getTime(line)
        elif 'main' in line and 'cut_max_mag' in line:
            t3 = getTime(line)
        elif 'main' in line and 'mass_distribution' in line:
            t4 = getTime(line)
        elif 'main' in line and 'mass_interp' in line:
            t5 = getTime(line)
        elif 'main' in line and 'binarity' in line:
            t6 = getTime(line)
        elif 'main' in line and 'completeness_rm' in line:
            t7 = getTime(line)
        elif 'main' in line and 'add_errors' in line:
            t8 = getTime(line)
        elif 'main' in line and 'likelihood' in line:
            t9 = getTime(line)

    return np.array([t0, t1, t2, t3, t4, t5, t6, t7, t8, t9])


def plot(Max_mass, Nmodels, times_all):
    """
    """
    times_norm = 100. * times_all / times_all.sum()
    cols = [
        'propModel', 'zaWAvrg', 'move', 'cut', 'M_dst', 'M_interp', 'binar',
        'complete', 'errors', 'likelihood']
    models_sec = Nmodels / max_time

    for i, c in enumerate(cols):
        print("{}: {:.2f}".format(c, times_norm[i]))

    fig, axs = plt.subplots(1, 1, figsize=(20, 10))

    axs.set_title(
        ("M_max={:.0f} | N={}, t={:.0f} --> [{:.0f} m/s]").format(
            Max_mass, Nmodels, max_time, models_sec), fontsize=18)
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

    profilerCall()

    # memrayCall()
